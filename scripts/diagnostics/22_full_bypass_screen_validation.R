# =============================================================================
# 22_full_bypass_screen_validation.R
#
# Script 21's full-array variance-filter bypass screen found 89 gene
# promoters at BH q<0.10 out of 18,131 tested (vs. 4,141 in the original
# screen) -- a much larger, more standard genome-wide-style candidate list,
# with the multiple-testing burden that implies. That is still an
# exploratory hit list, not a validated one, exactly like the original 5
# candidates were before scripts 15/20 tested them. This script runs the
# same two checks on the new list:
#
#   TEST A (permutation / negative control): shuffle which patient's
#   methylation profile gets compared against which patient's survival/
#   covariates -- breaking any real promoter-survival relationship while
#   leaving age/stage/purity/subtype correctly paired with true outcomes --
#   and rerun the IDENTICAL screen-and-BH-correct procedure across all
#   18,131 genes. If 89 hits is unremarkable under this null, the
#   candidate list does not clear the bar.
#
#   TEST B (Her2-exclusion stability): refit each candidate with Her2
#   patients dropped. A real, broadly-supported effect should survive
#   losing 33 patients / 6 events. A candidate whose significance depends
#   on the Her2 stratum specifically is not a stable finding.
#
# THIS SCRIPT IS EXPLORATORY, NOT PART OF THE FINALIZED PIPELINE -- kept
# under scripts/diagnostics/ (not deployed_scripts_v1/), same as script 21.
# Any candidate that passes both tests here still needs the same honest
# treatment CBX5 got: literature grounding, and skepticism proportional to
# how many genes were screened to find it (89 out of 18,131, on an 85-event
# cohort -- far more opportunity for a strong-looking but spurious hit than
# the original 5-candidate list had).
#
# RUNTIME WARNING: the stratified screen across 18,131 genes took about 1
# minute in script 21. TEST A reruns that screen once per permutation.
# N_PERM=100 could take 30-90+ minutes even parallelized, several times
# longer than script 15's original version (4,141 genes). Drop N_PERM to
# 20-30 for a quick sanity check before committing to a full run.
#
# Run from the project root, after 21_full_variance_filter_bypass_screen.R.
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(SummarizedExperiment)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(readxl)
library(doParallel)
library(foreach)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "full_bypass_screen_validation", base = "02_analysis")
set.seed(42)

N_PERM <- 100   # see runtime warning above -- drop to 20-30 for a first pass


# - - - 1. Get gene_mat_full -- reuse script 21's saved output if present - - -
# Rebuilding gene_mat_full from the raw SE (beta extraction, QC, M-value
# transform across ~450,000 probes) is the expensive part of script 21 and
# produces a DETERMINISTIC result given the same inputs, so there is no
# correctness reason to redo it here. Reuse the most recent run's saved
# matrix if one exists; rebuild from scratch only as a fallback.

find_latest_bypass_run <- function() {
  if (!dir.exists("02_analysis")) return(NULL)
  run_dirs <- list.dirs("02_analysis", recursive = FALSE)
  run_dirs <- run_dirs[grepl("full_variance_filter_bypass_screen$", basename(run_dirs))]
  if (length(run_dirs) == 0) return(NULL)
  run_dirs[order(basename(run_dirs), decreasing = TRUE)][1]
}

latest_run <- find_latest_bypass_run()
gene_mat_path <- if (!is.null(latest_run)) {
  file.path(latest_run, "models", "gene_mat_full_variance_filter_bypassed.rds")
} else NULL

meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")
barcode_col <- if ("sample_barcode" %in% names(meth_clinical)) {
  "sample_barcode"
} else if ("meth_sample_barcode" %in% names(meth_clinical)) {
  "meth_sample_barcode"
} else {
  stop("Neither sample_barcode nor meth_sample_barcode found in meth_clinical -- inspect names(meth_clinical)")
}

if (!is.null(gene_mat_path) && file.exists(gene_mat_path)) {
  cat("Reusing gene_mat_full from:", gene_mat_path, "\n")
  gene_mat_full <- readRDS(gene_mat_path)
} else {
  cat("No cached gene_mat_full found -- rebuilding from the raw SE",
      "(this replicates script 21's steps 1-5, expect several minutes)...\n")

  keep_barcodes <- meth_clinical[[barcode_col]]
  meth_se <- readRDS("00_raw_data/TCGA_BRCA_Methylation450_SE.rds")
  se_id15_all <- substr(colnames(meth_se), 1, 15)
  keep_order   <- order(colnames(meth_se))
  se_dedup_idx <- keep_order[!duplicated(se_id15_all[keep_order])]
  meth_se      <- meth_se[, se_dedup_idx]
  se_id15      <- substr(colnames(meth_se), 1, 15)
  clin_id15    <- substr(keep_barcodes, 1, 15)
  stopifnot(all(clin_id15 %in% se_id15))
  match_idx <- match(clin_id15, se_id15)
  meth_se_clean <- meth_se[, match_idx]
  colnames(meth_se_clean) <- keep_barcodes
  rm(meth_se); gc()

  assay_names <- names(assays(meth_se_clean))
  beta_mat <- if ("Beta_value" %in% assay_names) assay(meth_se_clean, "Beta_value")
    else if ("Beta_value_mean" %in% assay_names) assay(meth_se_clean, "Beta_value_mean")
    else if ("beta" %in% assay_names) assay(meth_se_clean, "beta")
    else if ("Beta" %in% assay_names) assay(meth_se_clean, "Beta")
    else assay(meth_se_clean, 1)

  miss_frac <- rowMeans(is.na(beta_mat))
  beta_mat  <- beta_mat[miss_frac <= 0.2, ]
  na_rows <- which(rowSums(is.na(beta_mat)) > 0)
  if (length(na_rows) > 0) {
    row_medians <- apply(beta_mat[na_rows, , drop = FALSE], 1, median, na.rm = TRUE)
    for (i in seq_along(na_rows)) {
      idx <- na_rows[i]
      beta_mat[idx, is.na(beta_mat[idx, ])] <- row_medians[i]
    }
  }
  probe_chr  <- as.character(seqnames(rowRanges(meth_se_clean)[rownames(beta_mat)]))
  beta_mat   <- beta_mat[!(probe_chr %in% c("chrX", "chrY")), ]
  rm(meth_se_clean); gc()

  beta_clipped <- pmin(pmax(beta_mat, 0.001), 0.999)
  rm(beta_mat); gc()
  m_mat_full <- log2(beta_clipped / (1 - beta_clipped))
  rm(beta_clipped); gc()

  ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann_matched_full <- ann_df[match(rownames(m_mat_full), rownames(ann_df)), ]
  is_promoter_full <- grepl("TSS200|TSS1500", ann_matched_full$UCSC_RefGene_Group)
  m_mat_promoters_full <- m_mat_full[is_promoter_full, ]
  gene_labels_full_clean <- sapply(
    strsplit(ann_matched_full$UCSC_RefGene_Name[is_promoter_full], ";"),
    function(x) x[1]
  )
  gene_sums_full   <- rowsum(m_mat_promoters_full, group = gene_labels_full_clean)
  gene_counts_full <- table(gene_labels_full_clean)
  gene_mat_full    <- gene_sums_full / as.vector(gene_counts_full)
  rm(m_mat_full, m_mat_promoters_full, gene_sums_full); gc()
}

cat("gene_mat_full:", nrow(gene_mat_full), "genes x", ncol(gene_mat_full), "samples\n")


# - - - 2. Analysis cohort (identical construction to scripts 14/15/21) - - -

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

cat("Analysis cohort:", sum(complete), "patients,", sum(meth_c$event), "events\n")


# - - - 3. Screening functions (explicit args, so TEST A can permute freely) - - -

stratified_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity) + strata(meth_x$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}
run_full_screen <- function(gm, surv_x, meth_x) {
  res <- t(apply(gm, 1, stratified_screen_one, surv_x = surv_x, meth_x = meth_x))
  colnames(res) <- c("coef", "z", "p")
  res
}


# - - - 4. Recompute the observed (real) screen - - -
# Recomputed fresh so this validation is testing the exact same code path
# and candidate list as script 21, not a hand-copied gene list.

cat("\nRunning observed stratified screen on", nrow(gene_mat_full_c), "genes...\n")
t0 <- Sys.time()
obs_res <- run_full_screen(gene_mat_full_c, surv_c, meth_c)
cat("Done in", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

obs_padj   <- p.adjust(obs_res[, "p"], method = "BH")
obs_n_hits <- sum(obs_padj < 0.10, na.rm = TRUE)
obs_min_p  <- min(obs_res[, "p"], na.rm = TRUE)
candidates <- rownames(gene_mat_full_c)[!is.na(obs_padj) & obs_padj < 0.10]

cat("Observed:", obs_n_hits, "genes at BH q<0.10\n")


# =============================================================================
# TEST A: Permutation null
# =============================================================================

cat("\n=== TEST A: permutation null (", N_PERM, "permutations) ===\n")

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Running on", n_cores, "cores...\n")

t0 <- Sys.time()
perm_results <- foreach(
  r = seq_len(N_PERM),
  .packages = c("survival"),
  .combine  = rbind
) %dopar% {
  set.seed(1000 + r)
  perm_idx <- sample(ncol(gene_mat_full_c))
  gm_perm  <- gene_mat_full_c[, perm_idx, drop = FALSE]
  res_p  <- run_full_screen(gm_perm, surv_c, meth_c)
  padj_p <- p.adjust(res_p[, "p"], method = "BH")
  c(n_hits = sum(padj_p < 0.10, na.rm = TRUE), min_p = min(res_p[, "p"], na.rm = TRUE))
}
stopCluster(cl)
cat("Permutation run:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

perm_results <- as.data.frame(perm_results)

perm_p_hitcount <- (sum(perm_results$n_hits >= obs_n_hits) + 1) / (N_PERM + 1)
perm_p_minp     <- (sum(perm_results$min_p <= obs_min_p) + 1) / (N_PERM + 1)

cat("\nNull distribution of hit counts (BH q<0.10) across", N_PERM, "permutations:\n")
print(summary(perm_results$n_hits))
cat(sprintf("\nObserved: %d hits. Permutation p-value: %.3f\n",
            obs_n_hits, perm_p_hitcount))
cat(sprintf("Observed smallest p-value: %.3g. Permutation p-value: %.3f\n",
            obs_min_p, perm_p_minp))

pdf(file.path(paths$diagnostics, "full_bypass_permutation_null.pdf"), width = 7, height = 5)
hist(perm_results$n_hits, breaks = max(10, N_PERM / 5), col = "grey80",
     border = "white", xlab = "Genes at BH q<0.10 (null permutations)",
     main = sprintf("Permutation null (n=%d) vs observed (perm p=%.3f)",
                     N_PERM, perm_p_hitcount),
     xlim = range(c(perm_results$n_hits, obs_n_hits)))
abline(v = obs_n_hits, col = "tomato", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Observed = %d", obs_n_hits),
       col = "tomato", lty = 2, lwd = 2, bty = "n")
dev.off()


# =============================================================================
# TEST B: Her2-exclusion stability
# =============================================================================

cat("\n=== TEST B: stability with Her2 excluded (", length(candidates), "candidates) ===\n")

no_her2 <- meth_c$subtype != "Her2"
meth_no_her2     <- meth_c[no_her2, ]
meth_no_her2$subtype <- droplevels(meth_no_her2$subtype)
gene_mat_no_her2 <- gene_mat_full_c[, no_her2, drop = FALSE]
surv_no_her2     <- surv_c[no_her2]

cat("Cohort with Her2 excluded:", sum(no_her2), "patients,",
    sum(meth_no_her2$event), "events\n")

stability <- data.frame(
  gene              = candidates,
  coef_with_her2    = obs_res[candidates, "coef"],
  p_with_her2       = obs_res[candidates, "p"],
  padj_with_her2    = obs_padj[candidates],
  coef_without_her2 = NA_real_,
  p_without_her2    = NA_real_,
  stringsAsFactors  = FALSE
)

for (g in candidates) {
  r <- stratified_screen_one(gene_mat_no_her2[g, ], surv_no_her2, meth_no_her2)
  stability[stability$gene == g, "coef_without_her2"] <- r["coef"]
  stability[stability$gene == g, "p_without_her2"]    <- r["p"]
}

stability$same_direction <- sign(stability$coef_with_her2) == sign(stability$coef_without_her2)
stability$still_p05      <- stability$p_without_her2 < 0.05
stability$stable_both    <- stability$same_direction & stability$still_p05

stability <- stability[order(stability$padj_with_her2), ]
print(stability)

cat(sprintf("\n%d of %d candidates (%.0f%%) are stable to Her2 exclusion",
            sum(stability$stable_both), nrow(stability),
            100 * mean(stability$stable_both)))
cat(" (same direction AND still nominally p<0.05 without Her2).\n")


# - - - 5. Save outputs - - -

write.csv(data.frame(permutation = seq_len(N_PERM), perm_results),
          file.path(paths$tables, "full_bypass_permutation_null_results.csv"),
          row.names = FALSE)
write.csv(stability, file.path(paths$tables, "full_bypass_her2_exclusion_stability.csv"),
          row.names = FALSE)

write_manifest(paths$run, list(
  script                 = "22_full_bypass_screen_validation.R",
  n_genes_tested         = nrow(gene_mat_full_c),
  n_permutations         = N_PERM,
  observed_hits_fdr10    = obs_n_hits,
  permutation_p_hitcount = round(perm_p_hitcount, 4),
  permutation_p_minp     = round(perm_p_minp, 4),
  n_candidates_tested    = length(candidates),
  n_stable_both_checks   = sum(stability$stable_both),
  seed                   = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
cat("Reminder: this validates the 89-candidate list COLLECTIVELY (TEST A)",
    "and each candidate's own stability (TEST B) -- it does not by itself",
    "make any single gene trustworthy the way CBX5 required literature",
    "grounding and the Her2 check on top of a much smaller, less",
    "multiple-testing-burdened original list. Still keep this out of",
    "deployed_scripts_v1/ until you've decided how to report a list this",
    "size responsibly.\n")
