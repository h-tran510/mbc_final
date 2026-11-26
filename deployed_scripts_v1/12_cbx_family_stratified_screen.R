library(survival)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

source("deployed_scripts_v1/00_utils.R")

paths <- make_run_dir(tag = "cbx_variance_filter_bypass", base = "02_analysis")
set.seed(42)

meth_se       <- readRDS("00_raw_data/TCGA_BRCA_Methylation450_SE.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

stopifnot("meth_clinical has no sample_barcode column -- rerun 07" =
            "sample_barcode" %in% names(meth_clinical))
keep_barcodes <- meth_clinical[["sample_barcode"]]

stopifnot("meth_clinical has no meth_sample_barcode column -- rerun 07" =
            "meth_sample_barcode" %in% names(meth_clinical))
assert_one_per_patient(meth_clinical$meth_sample_barcode, "meth_clinical rows", "rerun 07")
match_idx <- match(meth_clinical$meth_sample_barcode, colnames(meth_se))
stopifnot("Some meth_clinical methylation barcodes are not in the raw SE" =
            !anyNA(match_idx))
meth_se_clean <- meth_se[, match_idx]
colnames(meth_se_clean) <- keep_barcodes
cat("Raw SE subset to the finalized cohort:", ncol(meth_se_clean), "samples",
    "(matched on the exact methylation barcode)\n")

assay_names <- names(assays(meth_se_clean))
beta_assays <- intersect(c("Beta_value", "Beta_value_mean", "beta", "Beta"), assay_names)
if (length(beta_assays) > 0) {
  beta_mat <- assay(meth_se_clean, beta_assays[1])
} else {
  beta_mat <- assay(meth_se_clean, 1)
}
cat("Raw beta matrix:", nrow(beta_mat), "probes x", ncol(beta_mat), "samples\n")

cbx_genes <- c("CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8")
hp1_paralogs <- c("CBX1", "CBX3", "CBX5")
polycomb_cbx <- c("CBX2", "CBX4", "CBX6", "CBX7", "CBX8")

ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_raw <- ann_df[match(rownames(beta_mat), rownames(ann_df)), ]

cbx_pairs <- promoter_gene_pairs(ann_raw$UCSC_RefGene_Name, ann_raw$UCSC_RefGene_Group)
cbx_pairs <- cbx_pairs[cbx_pairs$gene %in% cbx_genes, ]

cat("\nCBX-family promoter probes found on the raw array (pre-QC):",
    length(unique(cbx_pairs$probe)), "probes,", nrow(cbx_pairs), "probe-gene pairs\n")

beta_cbx <- beta_mat[cbx_pairs$probe, , drop = FALSE]
ann_cbx  <- ann_raw[cbx_pairs$probe, ]

cat("Chromosomes represented:", paste(unique(as.character(ann_cbx$chr)), collapse = ", "), "\n")

gene_labels_cbx <- cbx_pairs$gene

cat("\nPer-gene raw promoter probe counts:\n")
print(table(gene_labels_cbx))

miss_frac <- rowMeans(is.na(beta_cbx))
keep_miss <- miss_frac <= 0.2
cat("\nProbes dropped for >20% missingness:", sum(!keep_miss), "of", length(keep_miss), "\n")
beta_cbx <- beta_cbx[keep_miss, , drop = FALSE]
gene_labels_cbx <- gene_labels_cbx[keep_miss]

na_rows <- which(rowSums(is.na(beta_cbx)) > 0)
if (length(na_rows) > 0) {
  row_medians <- apply(beta_cbx[na_rows, , drop = FALSE], 1, median, na.rm = TRUE)
  for (i in seq_along(na_rows)) {
    probe_row <- na_rows[i]
    missing_samples <- is.na(beta_cbx[probe_row, ])
    beta_cbx[probe_row, missing_samples] <- row_medians[i]
  }
  cat("Probes with imputed missing values:", length(na_rows), "\n")
}

cat("\nCBX promoter probes surviving to M-value transform:", nrow(beta_cbx), "\n")
cat("Per-gene counts after QC (excluding the variance filter):\n")
print(table(gene_labels_cbx))

beta_clipped <- pmin(pmax(beta_cbx, 0.001), 0.999)
m_cbx <- log2(beta_clipped / (1 - beta_clipped))
stopifnot(identical(colnames(m_cbx), keep_barcodes))

gene_sums   <- rowsum(m_cbx, group = gene_labels_cbx)
gene_counts <- table(gene_labels_cbx)

stopifnot(
  "rowsum/table group ordering diverged -- CBX promoter means would be misaligned" =
    identical(rownames(gene_sums), names(gene_counts))
)
gene_mat_cbx <- gene_sums / as.numeric(gene_counts)

cat("\nRecovered CBX gene_mat:", nrow(gene_mat_cbx), "genes x",
    ncol(gene_mat_cbx), "samples\n")
cat("Genes recovered:", paste(rownames(gene_mat_cbx), collapse = ", "), "\n")
cat("Still missing entirely:", paste(setdiff(cbx_genes, rownames(gene_mat_cbx)), collapse = ", "), "\n")

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

meth_clinical <- load_subtype_purity(meth_clinical)

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)
meth_c        <- meth_clinical[complete, ]
gene_mat_cbx_c <- gene_mat_cbx[, complete, drop = FALSE]
surv_c        <- surv_obj[complete]

cat("\nAnalysis cohort:", sum(complete), "patients,", sum(meth_c$event), "events\n")

cat("\n=== Individual genes, variance-filter bypassed, subtype+purity-stratified ===\n")

recovered_genes <- rownames(gene_mat_cbx_c)
res <- t(sapply(recovered_genes, function(g) {
  stratified_screen_one(gene_mat_cbx_c[g, ], surv_c, meth_c)
}))
res_df <- data.frame(
  gene = recovered_genes, coef = res[, "coef"], p = res[, "p"],
  n_probes_used = as.vector(gene_counts[recovered_genes]),
  status = ifelse(recovered_genes == "CBX5",
                   "CONFIRMATORY (discovered in this data via script 10)",
                   "BLIND (never screened before this run)"),
  stringsAsFactors = FALSE
)
blind <- recovered_genes != "CBX5"
res_df$padj <- NA_real_
res_df$padj[blind] <- p.adjust(res_df$p[blind], method = "BH")
cat("BH padj covers the", sum(blind), "BLIND genes only; CBX5 is reported",
    "unadjusted (padj = NA) because it is confirmatory (discovered in script 10).\n")
print(res_df[order(res_df$p), ])

hp1_present <- intersect(hp1_paralogs, recovered_genes)
if (length(hp1_present) >= 2) {
  cat("\n=== Combined HP1-paralog score ===\n")
  hp1_score <- rowMeans(scale(t(gene_mat_cbx_c[hp1_present, , drop = FALSE])))
  stopifnot("hp1_score must have one value per patient" =
              length(hp1_score) == ncol(gene_mat_cbx_c))
  r <- stratified_screen_one(hp1_score, surv_c, meth_c)
  cat(sprintf("HP1 combined (%s): coef=%.3f, p=%.3g\n",
              paste(hp1_present, collapse = "+"), r["coef"], r["p"]))

  blind_only <- setdiff(hp1_present, "CBX5")
  if (length(blind_only) >= 1) {
    blind_score <- if (length(blind_only) == 1) {
      scale(gene_mat_cbx_c[blind_only, ])[, 1]
    } else {
      rowMeans(scale(t(gene_mat_cbx_c[blind_only, , drop = FALSE])))
    }
    stopifnot("blind_score must have one value per patient" =
                length(blind_score) == ncol(gene_mat_cbx_c))
    r_blind <- stratified_screen_one(blind_score, surv_c, meth_c)
    cat(sprintf("CBX5-EXCLUDED, fully blind (%s): coef=%.3f, p=%.3g  <- the genuinely independent number\n",
                paste(blind_only, collapse = "+"), r_blind["coef"], r_blind["p"]))
  }
}

polycomb_present <- intersect(polycomb_cbx, recovered_genes)
if (length(polycomb_present) >= 2) {
  cat("\n=== Combined Polycomb-CBX score ===\n")
  poly_score <- rowMeans(scale(t(gene_mat_cbx_c[polycomb_present, , drop = FALSE])))
  stopifnot("poly_score must have one value per patient" =
              length(poly_score) == ncol(gene_mat_cbx_c))
  r_poly <- stratified_screen_one(poly_score, surv_c, meth_c)
  cat(sprintf("Polycomb combined (%s): coef=%.3f, p=%.3g\n",
              paste(polycomb_present, collapse = "+"), r_poly["coef"], r_poly["p"]))
}

write.csv(res_df, file.path(paths$tables, "cbx_recovered_individual_results.csv"),
          row.names = FALSE)
saveRDS(gene_mat_cbx, file.path(paths$models, "cbx_gene_mat_variance_filter_bypassed.rds"))

write_manifest(paths$run, list(
  script = "12_cbx_family_stratified_screen.R",
  genes_recovered = paste(recovered_genes, collapse = ";"),
  genes_still_missing = paste(setdiff(cbx_genes, recovered_genes), collapse = ";"),
  variance_filter = "BYPASSED for this panel only",
  seed = 42
), inputs = c("00_raw_data/TCGA_BRCA_Methylation450_SE.rds",
              "01_processed_data/meth_clinical.rds",
              unname(SUBTYPE_PURITY_INPUTS)))

cat("\n=== Done. Results in", paths$run, "===\n")
cat("Remember: this bypasses the main pipeline's variance filter for CBX genes",
    "only. Report it as a targeted candidate-panel analysis, not as a",
    "correction to the discovery screen. No coherent family-wide signature",
    "was found -- CBX5 and CBX2 are two individual findings, not one.\n")
