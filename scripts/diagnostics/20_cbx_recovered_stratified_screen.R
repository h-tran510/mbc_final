# =============================================================================
# 20_cbx_recovered_stratified_screen.R
#
# Scripts 18-19 established that CBX1/CBX3/CBX4/CBX6/CBX8 are completely
# absent from meth_m_matrix.rds, and CBX2/CBX5/CBX7 are down to 1-5 of their
# original 9-28 promoter probes. The cause (confirmed by reading
# 07_methylation_prepare.R step 6d) is NOT array design or cross-reactivity --
# it is an unsupervised, genome-wide "keep the top 10% most variable probes"
# filter applied before any gene-level analysis happens. CBX-family promoter
# CpGs simply didn't rank in the top 10% by variance across the 755-patient
# cohort.
#
# This script rebuilds the CBX-family promoter matrix from the RAW SE,
# replicating every QC step in 07_methylation_prepare.R EXCEPT the variance
# filter, for these 8 genes only. Skipping a genome-wide unsupervised filter
# for a small, pre-specified candidate panel is standard practice for
# targeted/hypothesis-driven analysis -- discovery-screen filters are built to
# manage a 450,000-test multiple-comparison problem that a 5-8 gene panel
# doesn't have. Report it explicitly as a departure from the main pipeline's
# QC standard, not as "the real result the main pipeline missed."
#
# CAUTION ON INTERPRETATION: low baseline variance mechanically limits the
# power of a between-patient regression to detect anything, real or not. A
# null result for a recovered gene here does not mean "not biologically
# relevant" -- it may just mean there wasn't enough between-patient spread to
# find out either way. Recall also the circularity point from script 17:
# CBX5 was discovered in this exact cohort; its recovered result here is
# still confirmatory, not independent.
#
# Run from the project root, after scripts 14/15/17.
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "cbx_variance_filter_bypass", base = "02_analysis")
set.seed(42)


# - - - 1. Load raw SE and the already-finalized clinical cohort - - -
# meth_clinical.rds already encodes every sample-level decision from
# 07_methylation_prepare.R (clinical matching, PFI validity, stage/age
# completeness) -- reuse it rather than re-deriving the cohort, so the
# patient set here is IDENTICAL to scripts 13-17's, just with different
# probes.

meth_se       <- readRDS("00_raw_data/TCGA_BRCA_Methylation450_SE.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

barcode_col <- if ("sample_barcode" %in% names(meth_clinical)) {
  "sample_barcode"
} else if ("meth_sample_barcode" %in% names(meth_clinical)) {
  "meth_sample_barcode"
} else {
  stop("Neither sample_barcode nor meth_sample_barcode found in meth_clinical -- inspect names(meth_clinical)")
}
cat("Using clinical barcode column:", barcode_col, "\n")

keep_barcodes <- meth_clinical[[barcode_col]]

# Exact-string matching against colnames(meth_se) failed -- the same
# barcode-suffix issue this project has hit before (vial/portion/analyte
# codes can differ between the raw SE and a saved clinical table even when
# both trace back to the same samples). Match on a truncated ID instead,
# the same approach run_alignment_checks() already uses elsewhere in this
# project, and fail loudly if truncation doesn't resolve it either.
se_id15_all <- substr(colnames(meth_se), 1, 15)

# The raw SE has repeat aliquots/portions for some samples -- multiple full
# barcodes collapsing to the same 15-char sample ID. Same situation
# 07_methylation_prepare.R resolves at the patient level ("keep the
# lexicographically first barcode"); apply the identical rule here at the
# sample level, so exactly one raw-SE column survives per 15-char ID.
n_dup <- sum(duplicated(se_id15_all))
if (n_dup > 0) {
  cat("Raw SE has", n_dup, "repeat aliquot(s)/portion(s) collapsing to an",
      "already-seen 15-char sample ID -- keeping the lexicographically",
      "first barcode for each, matching 07_methylation_prepare.R's",
      "convention.\n")
}
keep_order   <- order(colnames(meth_se))
se_dedup_idx <- keep_order[!duplicated(se_id15_all[keep_order])]
meth_se      <- meth_se[, se_dedup_idx]
se_id15      <- substr(colnames(meth_se), 1, 15)
stopifnot("Deduplication left duplicate truncated IDs in the raw SE" =
            sum(duplicated(se_id15)) == 0)

clin_id15 <- substr(keep_barcodes, 1, 15)
stopifnot(
  "Duplicate truncated IDs in meth_clinical" = sum(duplicated(clin_id15)) == 0,
  "Some clinical patients still not found in the raw SE after truncation" =
    all(clin_id15 %in% se_id15)
)

match_idx <- match(clin_id15, se_id15)
meth_se_clean <- meth_se[, match_idx]
# Relabel columns with the clinical table's own barcode so every downstream
# join in this script (subtype, purity) uses one consistent barcode source.
colnames(meth_se_clean) <- keep_barcodes
cat("Raw SE subset to the finalized cohort:", ncol(meth_se_clean), "samples",
    "(matched via 15-char truncated barcode)\n")


# - - - 2. Extract beta values (same assay-detection logic as 07) - - -

assay_names <- names(assays(meth_se_clean))
if ("Beta_value" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "Beta_value")
} else if ("Beta_value_mean" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "Beta_value_mean")
} else if ("beta" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "beta")
} else if ("Beta" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "Beta")
} else {
  beta_mat <- assay(meth_se_clean, 1)
}
cat("Raw beta matrix:", nrow(beta_mat), "probes x", ncol(beta_mat), "samples\n")


# - - - 3. Restrict to CBX-family promoter probes, BEFORE any filtering - - -
# This is done first specifically so the missing-value and imputation steps
# below only ever touch this small probe set -- no risk of the variance
# filter sneaking back in through some other code path.

cbx_genes <- c("CBX1", "CBX2", "CBX3", "CBX4", "CBX5", "CBX6", "CBX7", "CBX8")
hp1_paralogs <- c("CBX1", "CBX3", "CBX5")
polycomb_cbx <- c("CBX2", "CBX4", "CBX6", "CBX7", "CBX8")

ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_raw <- ann_df[match(rownames(beta_mat), rownames(ann_df)), ]

is_cbx_gene <- sapply(strsplit(ann_raw$UCSC_RefGene_Name, ";"), function(x) any(x %in% cbx_genes))
is_promoter <- grepl("TSS200|TSS1500", ann_raw$UCSC_RefGene_Group)
cbx_probe_mask <- !is.na(is_cbx_gene) & is_cbx_gene & is_promoter

cat("\nCBX-family promoter probes found on the raw array (pre-QC):",
    sum(cbx_probe_mask), "\n")

beta_cbx <- beta_mat[cbx_probe_mask, , drop = FALSE]
ann_cbx  <- ann_raw[cbx_probe_mask, ]

# Chromosome sanity check -- confirms these are all autosomal, so the sex-
# chromosome removal step in 07 (6c) would never have touched them anyway.
cat("Chromosomes represented:", paste(unique(as.character(ann_cbx$chr)), collapse = ", "), "\n")

gene_labels_cbx <- sapply(strsplit(ann_cbx$UCSC_RefGene_Name, ";"), function(x) {
  hit <- intersect(x, cbx_genes)
  if (length(hit) > 0) hit[1] else NA
})

cat("\nPer-gene raw promoter probe counts:\n")
print(table(gene_labels_cbx))


# - - - 4. Replicate 07's QC steps 6a/6b, SKIP 6c (n/a) and 6d (the bypass) - - -

# (a) Missing-value filter: >20% NA removed, same threshold as the main pipeline
miss_frac <- rowMeans(is.na(beta_cbx))
keep_miss <- miss_frac <= 0.2
cat("\nProbes dropped for >20% missingness:", sum(!keep_miss), "of", length(keep_miss), "\n")
beta_cbx <- beta_cbx[keep_miss, , drop = FALSE]
gene_labels_cbx <- gene_labels_cbx[keep_miss]

# (b) Median imputation for any remaining NAs
na_rows <- which(rowSums(is.na(beta_cbx)) > 0)
if (length(na_rows) > 0) {
  row_medians <- apply(beta_cbx[na_rows, , drop = FALSE], 1, median, na.rm = TRUE)
  for (i in seq_along(na_rows)) {
    idx <- na_rows[i]
    beta_cbx[idx, is.na(beta_cbx[idx, ])] <- row_medians[i]
  }
  cat("Probes with imputed missing values:", length(na_rows), "\n")
}

# (c) Sex chromosome removal -- not applicable, confirmed autosomal above.
# (d) DELIBERATELY SKIPPED. This is the entire point of this script.

cat("\nCBX promoter probes surviving to M-value transform:", nrow(beta_cbx), "\n")
cat("Per-gene counts after QC (excluding the variance filter):\n")
print(table(gene_labels_cbx))


# - - - 5. M-value transform, identical to 07 step 7 - - -

beta_clipped <- pmin(pmax(beta_cbx, 0.001), 0.999)
m_cbx <- log2(beta_clipped / (1 - beta_clipped))
stopifnot(identical(colnames(m_cbx), keep_barcodes))


# - - - 6. Aggregate to gene-level promoter values - - -
# Same convention as scripts 13/14/15/17: rowsum by gene, divide by probe
# count per gene -- an unweighted mean of M-values across each gene's
# surviving promoter probes.

gene_sums   <- rowsum(m_cbx, group = gene_labels_cbx)
gene_counts <- table(gene_labels_cbx)
gene_mat_cbx <- gene_sums / as.vector(gene_counts)

cat("\nRecovered CBX gene_mat:", nrow(gene_mat_cbx), "genes x",
    ncol(gene_mat_cbx), "samples\n")
cat("Genes recovered:", paste(rownames(gene_mat_cbx), collapse = ", "), "\n")
cat("Still missing entirely:", paste(setdiff(cbx_genes, rownames(gene_mat_cbx)), collapse = ", "), "\n")


# - - - 7. Build the analysis cohort (complete-case on covariates) - - -

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

pam50_path <- "01_processed_data/PAM50_classification.txt"
subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)
names(subtype_tab) <- c("barcode", "subtype")
subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)
subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]
meth_patient <- substr(meth_clinical[[barcode_col]], 1, 12)
meth_clinical$subtype <- factor(subtype_tab$subtype[match(meth_patient, subtype_tab$patient)])

library(readxl)
purity_tab <- as.data.frame(read_excel("01_processed_data/aran_purity_tcga.xlsx",
                                       sheet = "Supp Data 1", skip = 2))
names(purity_tab)[names(purity_tab) == "Sample ID"] <- "sample_barcode"
purity_tab$CPE <- as.numeric(purity_tab$CPE)
purity_tab <- purity_tab[purity_tab$`Cancer type` == "BRCA" &
                           !is.na(purity_tab$`Cancer type`), ]
purity_tab$sample_id_15 <- substr(purity_tab$sample_barcode, 1, 15)
meth_clinical$purity <- purity_tab$CPE[
  match(substr(meth_clinical[[barcode_col]], 1, 15), purity_tab$sample_id_15)
]

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)
meth_c        <- meth_clinical[complete, ]
gene_mat_cbx_c <- gene_mat_cbx[, complete, drop = FALSE]
surv_c        <- surv_obj[complete]

cat("\nAnalysis cohort:", sum(complete), "patients,", sum(meth_c$event), "events\n")


# - - - 8. Screen (same function as scripts 14/15/17) - - -

stratified_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity) + strata(meth_x$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

cat("\n=== Individual genes, variance-filter bypassed, subtype+purity-stratified ===\n")

recovered_genes <- rownames(gene_mat_cbx_c)
res <- t(sapply(recovered_genes, function(g) {
  stratified_screen_one(gene_mat_cbx_c[g, ], surv_c, meth_c)
}))
res_df <- data.frame(
  gene = recovered_genes, coef = res[, "coef"], p = res[, "p"],
  n_probes_used = as.vector(gene_counts[recovered_genes]),
  status = ifelse(recovered_genes == "CBX5",
                   "CONFIRMATORY (discovered in this data via script 14)",
                   "BLIND (never screened before this run)"),
  stringsAsFactors = FALSE
)
res_df$padj <- p.adjust(res_df$p, method = "BH")
print(res_df[order(res_df$p), ])


# - - - 9. Combined scores, now buildable with full HP1 coverage - - -

hp1_present <- intersect(hp1_paralogs, recovered_genes)
if (length(hp1_present) >= 2) {
  cat("\n=== Combined HP1-paralog score (recovered coverage) ===\n")
  # rowMeans, not colMeans: t(...) is patients x genes, scale() z-scores each
  # gene (column); averaging ACROSS genes for each patient needs rowMeans.
  # colMeans here would average down each column instead, collapsing to one
  # number per GENE rather than one score per patient -- the same bug already
  # caught once in script 13's B1 proliferation section.
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
  cat("\n=== Combined Polycomb-CBX score (recovered coverage) ===\n")
  poly_score <- rowMeans(scale(t(gene_mat_cbx_c[polycomb_present, , drop = FALSE])))
  stopifnot("poly_score must have one value per patient" =
              length(poly_score) == ncol(gene_mat_cbx_c))
  r_poly <- stratified_screen_one(poly_score, surv_c, meth_c)
  cat(sprintf("Polycomb combined (%s): coef=%.3f, p=%.3g\n",
              paste(polycomb_present, collapse = "+"), r_poly["coef"], r_poly["p"]))
}


# - - - 10. Save outputs - - -

write.csv(res_df, file.path(paths$tables, "cbx_recovered_individual_results.csv"),
          row.names = FALSE)
saveRDS(gene_mat_cbx, file.path(paths$models, "cbx_gene_mat_variance_filter_bypassed.rds"))

write_manifest(paths$run, list(
  script = "20_cbx_recovered_stratified_screen.R",
  genes_recovered = paste(recovered_genes, collapse = ";"),
  genes_still_missing = paste(setdiff(cbx_genes, recovered_genes), collapse = ";"),
  variance_filter = "BYPASSED for this panel only",
  seed = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
cat("Remember: this bypasses the main pipeline's variance filter for CBX genes",
    "only. Report it as a targeted candidate-panel analysis, not as a",
    "correction to the discovery screen. Low pre-filter variance may still",
    "limit power for genes that remain non-significant here.\n")
