# =============================================================================
# 21_full_variance_filter_bypass_screen.R
#
# Scripts 14/15 only ever screened the ~4,141 gene promoters whose CpG
# probes survived 07_methylation_prepare.R step 6(d): an unsupervised,
# genome-wide "keep the top 10% most variable probes" filter, applied
# before any gene identity is considered. Any gene whose promoter probes
# didn't rank in that top 10% was never in the screen at all -- not tested
# and found null, just silently absent. Script 20 showed this can matter
# even for genes that PARTIALLY survived: CBX2 had only 1 of 9 promoter
# probes in the original gene_mat and did not reach significance; with full
# 9-probe coverage recovered, it did (p=0.003). The same effect could be
# hiding elsewhere across the ~90% of the array the variance filter
# discarded.
#
# This script extends script 20's per-gene bypass to the WHOLE array:
# rebuild the promoter matrix from the raw SE, replicating every QC step in
# 07_methylation_prepare.R EXCEPT the variance filter, for ALL genes (not
# a curated panel), then rerun the naive-vs-stratified screen from
# 14_subtype_purity_stratified_screen.R on the result.
#
# THIS SCRIPT IS EXPLORATORY, NOT PART OF THE FINALIZED PIPELINE -- kept
# under scripts/diagnostics/ (not deployed_scripts_v1/) for exactly that
# reason. Any candidate it finds needs the same permutation + Her2-
# exclusion validation CBX5 got in script 15 before being trusted or moved
# into the deployed pipeline.
#
# SCOPE, READ BEFORE RUNNING:
#   - This is a MUCH bigger multiple-testing problem than script 14's --
#     likely ~15,000-20,000 genes instead of 4,141 (the 450K array design
#     covers ~99% of RefSeq genes). BH correction will be more conservative
#     accordingly, and that is correct, not a flaw to work around.
#   - Memory: building beta/M-value matrices at ~450,000 probes x 755
#     samples uses several GB of RAM per matrix, with multiple intermediate
#     copies. If this fails or thrashes on your machine, that is a hardware
#     limit, not a bug -- consider processing in probe-count chunks.
#   - Runtime: the screen itself (~20,000 Cox fits x 2 screens) should take
#     on the order of 15-30 minutes sequentially, similar order of
#     magnitude to script 14's "~4,100 fits in a few minutes" scaled up.
#     PERMUTATION VALIDATION IS DELIBERATELY NOT RUN HERE -- at this test
#     count it would take hours. Any new candidate found below needs its
#     own targeted permutation + Her2-exclusion check (same pattern as
#     script 15) before being trusted, exactly like CBX5 and CBX2 were.
#   - This is now a genuine, undirected discovery screen again, with the
#     multiple-testing burden that implies -- report it as such, not as a
#     confirmation exercise.
#
# Run from the project root.
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(readxl)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "full_variance_filter_bypass_screen", base = "02_analysis")
set.seed(42)


# - - - 1. Load raw SE and the already-finalized clinical cohort - - -
# Same barcode-matching logic as script 20: exact-string matching against
# colnames(meth_se) can fail on vial/portion/analyte suffix differences,
# and the raw SE can have repeat aliquots collapsing to the same 15-char
# sample ID.

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
se_id15_all   <- substr(colnames(meth_se), 1, 15)

n_dup <- sum(duplicated(se_id15_all))
if (n_dup > 0) {
  cat("Raw SE has", n_dup, "repeat aliquot(s)/portion(s) -- keeping the",
      "lexicographically first barcode per 15-char ID, matching",
      "07_methylation_prepare.R's convention.\n")
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
colnames(meth_se_clean) <- keep_barcodes
rm(meth_se); gc()
cat("Raw SE subset to the finalized cohort:", ncol(meth_se_clean), "samples\n")


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


# - - - 3. Replicate 07's QC steps 6a/6b/6c, SKIP 6d (the bypass) - - -

# (a) Missing-value filter: >20% NA removed, same threshold as the main pipeline
miss_frac <- rowMeans(is.na(beta_mat))
beta_mat  <- beta_mat[miss_frac <= 0.2, ]
cat("Probes after missing-value filter:", nrow(beta_mat), "\n")

# (b) Median imputation for any remaining NAs
na_rows <- which(rowSums(is.na(beta_mat)) > 0)
if (length(na_rows) > 0) {
  row_medians <- apply(beta_mat[na_rows, , drop = FALSE], 1, median, na.rm = TRUE)
  for (i in seq_along(na_rows)) {
    idx <- na_rows[i]
    beta_mat[idx, is.na(beta_mat[idx, ])] <- row_medians[i]
  }
  cat("Probes with imputed missing values:", length(na_rows), "\n")
}

# (c) Sex chromosome removal -- DOES matter here (full array, not just
# autosomal CBX genes), so unlike script 20 this step is not skipped.
probe_chr  <- as.character(seqnames(rowRanges(meth_se_clean)[rownames(beta_mat)]))
sex_probes <- probe_chr %in% c("chrX", "chrY")
beta_mat   <- beta_mat[!sex_probes, ]
cat("Probes after sex chromosome removal:", nrow(beta_mat), "\n")

# (d) DELIBERATELY SKIPPED. This is the entire point of this script.
cat("Variance filter (step 6d) skipped for this run.\n")

rm(meth_se_clean); gc()


# - - - 4. M-value transform, identical to 07 step 7 - - -

beta_clipped <- pmin(pmax(beta_mat, 0.001), 0.999)
rm(beta_mat); gc()
m_mat_full <- log2(beta_clipped / (1 - beta_clipped))
rm(beta_clipped); gc()
stopifnot(identical(colnames(m_mat_full), keep_barcodes))

cat("Full (unfiltered) M-value matrix:", nrow(m_mat_full), "probes x",
    ncol(m_mat_full), "samples\n")


# - - - 5. Promoter aggregation across the FULL probe set - - -

ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann_matched_full <- ann_df[match(rownames(m_mat_full), rownames(ann_df)), ]

is_promoter_full <- grepl("TSS200|TSS1500", ann_matched_full$UCSC_RefGene_Group)
m_mat_promoters_full <- m_mat_full[is_promoter_full, ]
gene_labels_full <- ann_matched_full$UCSC_RefGene_Name[is_promoter_full]
gene_labels_full_clean <- sapply(strsplit(gene_labels_full, split = ";"), function(x) x[1])

gene_sums_full   <- rowsum(m_mat_promoters_full, group = gene_labels_full_clean)
gene_counts_full <- table(gene_labels_full_clean)
gene_mat_full    <- gene_sums_full / as.vector(gene_counts_full)

cat("\nFull (variance-filter-bypassed) gene_mat:", nrow(gene_mat_full), "genes\n")

# Quantify the blind spot: compare against the variance-filtered gene_mat
# scripts 13/14/15 actually screened.
m_mat_orig <- readRDS("01_processed_data/meth_m_matrix.rds")
ann_matched_orig <- ann_df[match(rownames(m_mat_orig), rownames(ann_df)), ]
is_promoter_orig <- grepl("TSS200|TSS1500", ann_matched_orig$UCSC_RefGene_Group)
gene_labels_orig <- sapply(strsplit(ann_matched_orig$UCSC_RefGene_Name[is_promoter_orig], ";"),
                            function(x) x[1])
genes_originally_screened <- unique(gene_labels_orig)

newly_testable <- setdiff(rownames(gene_mat_full), genes_originally_screened)
cat(sprintf("\nGenes in the original variance-filtered screen: %d\n",
            length(genes_originally_screened)))
cat(sprintf("Genes now testable with full coverage: %d\n", nrow(gene_mat_full)))
cat(sprintf("Genes with a promoter that were COMPLETELY ABSENT from the",
            "original screen: %d (%.1f%% of the full set)\n",
            length(newly_testable), 100 * length(newly_testable) / nrow(gene_mat_full)))
rm(m_mat_orig, m_mat_full, m_mat_promoters_full, gene_sums_full); gc()


# - - - 6. Analysis cohort (identical construction to scripts 14/15) - - -

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

pam50_path <- "01_processed_data/PAM50_classification.txt"
subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)
names(subtype_tab) <- c("barcode", "subtype")
subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)
subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]
meth_patient <- substr(meth_clinical[[barcode_col]], 1, 12)
meth_clinical$subtype <- factor(subtype_tab$subtype[match(meth_patient, subtype_tab$patient)])

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
meth_c          <- meth_clinical[complete, ]
gene_mat_full_c <- gene_mat_full[, complete, drop = FALSE]
surv_c          <- surv_obj[complete]

cat("\nAnalysis cohort:", sum(complete), "patients,", sum(meth_c$event), "events\n")


# - - - 7. Run the naive and stratified screens across ALL genes - - -
# Identical model specification to script 14, just a much larger gene set.

naive_screen_one <- function(g) {
  fit <- try(coxph(surv_c ~ scale(g) + meth_c$age_years + meth_c$stage +
                      scale(meth_c$purity)), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}
stratified_screen_one <- function(g) {
  fit <- try(coxph(surv_c ~ scale(g) + meth_c$age_years + meth_c$stage +
                      scale(meth_c$purity) + strata(meth_c$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

n_genes <- nrow(gene_mat_full_c)
cat("\nRunning naive screen on", n_genes, "genes (this will take a while)...\n")
t0 <- Sys.time()
naive_res <- t(apply(gene_mat_full_c, 1, naive_screen_one))
cat("Naive screen:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

cat("\nRunning subtype+purity-stratified screen on", n_genes, "genes...\n")
t0 <- Sys.time()
strat_res <- t(apply(gene_mat_full_c, 1, stratified_screen_one))
cat("Stratified screen:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

results_full <- data.frame(
  gene            = rownames(gene_mat_full_c),
  naive_coef      = naive_res[, "coef"],
  naive_p         = naive_res[, "p"],
  stratified_coef = strat_res[, "coef"],
  stratified_p    = strat_res[, "p"],
  previously_screened = rownames(gene_mat_full_c) %in% genes_originally_screened,
  stringsAsFactors = FALSE
)
results_full$stratified_padj <- p.adjust(results_full$stratified_p, method = "BH")

cat(sprintf("\nBH correction now applied across %d tests (vs %d in the original screen) --",
            n_genes, length(genes_originally_screened)))
cat(" the significance bar is correspondingly stricter, correctly.\n")


# - - - 8. Candidates, split by whether they were previously testable - - -

candidates_full <- results_full[!is.na(results_full$stratified_padj) &
                                   results_full$stratified_padj < 0.10, ]
candidates_full <- candidates_full[order(candidates_full$stratified_padj), ]

cat("\n=== Candidates at BH q<0.10 across the FULL gene set ===\n")
print(candidates_full)

new_candidates <- candidates_full[!candidates_full$previously_screened, ]
recovered_candidates <- candidates_full[candidates_full$previously_screened &
                                           !(candidates_full$gene %in% c("C1orf86", "PLA2G3", "CBX5", "KLK9", "RAB6B")), ]

cat("\n=== Genuinely NEW candidates (had ZERO coverage in the original screen) ===\n")
print(new_candidates)
cat("\n=== Candidates that WERE in the original screen but not in its top-5",
    "list (may have gained power from fuller probe coverage here) ===\n")
print(recovered_candidates)

cat("\nNOTE: none of these are validated yet. Before trusting any of them,",
    "run the same permutation-null and Her2-exclusion checks used for",
    "CBX5 in 15_stratified_screen_validation.R, restricted to this",
    "candidate list -- do not report them on the strength of this screen",
    "alone, exactly as was true for the original 5 candidates. Do not move",
    "any result from this script into deployed_scripts_v1/ until that",
    "validation is done.\n")


# - - - 9. Save outputs - - -

write.csv(results_full, file.path(paths$tables, "full_bypass_screen_results.csv"),
          row.names = FALSE)
write.csv(candidates_full, file.path(paths$tables, "full_bypass_screen_candidates_fdr10.csv"),
          row.names = FALSE)
saveRDS(gene_mat_full, file.path(paths$models, "gene_mat_full_variance_filter_bypassed.rds"))

write_manifest(paths$run, list(
  script = "21_full_variance_filter_bypass_screen.R",
  n_genes_tested = n_genes,
  n_genes_original_screen = length(genes_originally_screened),
  n_genes_never_before_tested = length(newly_testable),
  n_candidates_fdr10 = nrow(candidates_full),
  n_new_candidates = nrow(new_candidates),
  seed = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
