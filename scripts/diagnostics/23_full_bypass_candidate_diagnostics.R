# =============================================================================
# 23_full_bypass_candidate_diagnostics.R
#
# Script 22 validated the 91-candidate list from the full variance-filter
# bypass screen (script 21) with a permutation null and Her2-exclusion
# check. Both nominally passed, but two things in those results are red
# flags worth chasing before trusting the list:
#
#   1. The permutation null (script 22, TEST A) is wildly overdispersed --
#      min 0, median 0, mean 43, MAX 4118 hits across 100 permutations.
#      That is the signature of highly correlated (non-independent) tests,
#      not a clean null. When one shuffle happens to pair an unusual
#      patient's outcome with their profile, many correlated promoters can
#      cross the significance line together -- so BH's FDR guarantee is
#      less trustworthy here than it was on the original, smaller screen.
#
#   2. 88 of 91 candidates (96.7%) share the same coefficient sign
#      (negative). That is far more extreme than the original screen's
#      already-unresolved ~74% negative sign bias, and is the classic
#      signature of one shared underlying axis driving many "independent"
#      hits at once, rather than 91 distinct biological findings.
#
# Both point the same direction: bypassing the variance filter for the
# whole array likely reintroduced a large block of promoters that load on
# the same residual PC1/global-methylation axis B7 (09_pc1_biological_
# characterization.R) found still partially unexplained by subtype+purity
# (R^2 = 0.60, 40% left over). This script tests that directly:
#
#   - Spearman correlation of each candidate against the already-
#     characterized PC1 score.
#   - Refit each candidate's stratified model WITH PC1 added as an extra
#     covariate (same logic as B7, applied gene-by-gene): does the
#     candidate survive once PC1 is controlled for on top of subtype and
#     purity?
#   - Formal sign-bias test (binom.test) on the candidate list, for the
#     record.
#
# THIS SCRIPT IS EXPLORATORY, kept under scripts/diagnostics/, same as
# scripts 21-22.
#
# Run from the project root, after 21_full_variance_filter_bypass_screen.R
# and 22_full_bypass_screen_validation.R.
# =============================================================================


# - - - 0. Setup - - -

library(survival)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "full_bypass_candidate_diagnostics", base = "02_analysis")


# - - - 1. Locate the most recent bypass-screen run and PC1 characterization run - - -

find_latest_run <- function(pattern) {
  if (!dir.exists("02_analysis")) return(NULL)
  run_dirs <- list.dirs("02_analysis", recursive = FALSE)
  run_dirs <- run_dirs[grepl(pattern, basename(run_dirs))]
  if (length(run_dirs) == 0) return(NULL)
  run_dirs[order(basename(run_dirs), decreasing = TRUE)][1]
}

bypass_run <- find_latest_run("full_variance_filter_bypass_screen$")
pc1_run    <- find_latest_run("pc1")

stopifnot("No full_variance_filter_bypass_screen run found -- run script 21 first" =
            !is.null(bypass_run))
stopifnot("No PC1 characterization run found -- run 09/13_pc1_biological_characterization.R first" =
            !is.null(pc1_run))

cat("Using bypass-screen run:", bypass_run, "\n")
cat("Using PC1 run:", pc1_run, "\n")

gene_mat_full <- readRDS(file.path(bypass_run, "models", "gene_mat_full_variance_filter_bypassed.rds"))
candidates_df <- read.csv(file.path(bypass_run, "tables", "full_bypass_screen_candidates_fdr10.csv"),
                           stringsAsFactors = FALSE)

pc1_clinical_path <- file.path(pc1_run, "meth_clinical_with_pc1.rds")
stopifnot("meth_clinical_with_pc1.rds not found in the PC1 run directory" =
            file.exists(pc1_clinical_path))
meth_clinical_pc1 <- readRDS(pc1_clinical_path)

cat("Candidates loaded:", nrow(candidates_df), "\n")
cat("PC1 available for:", sum(!is.na(meth_clinical_pc1$pc1)), "patients\n")


# - - - 2. Align PC1 onto the current cohort - - -

barcode_col <- if ("sample_barcode" %in% names(meth_clinical_pc1)) {
  "sample_barcode"
} else if ("meth_sample_barcode" %in% names(meth_clinical_pc1)) {
  "meth_sample_barcode"
} else {
  stop("Neither sample_barcode nor meth_sample_barcode found -- inspect names(meth_clinical_pc1)")
}

meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")
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

# Join PC1 in via 15-char truncated barcode (same convention used throughout
# this project for cross-file joins).
pc1_id15  <- substr(meth_clinical_pc1[[barcode_col]], 1, 15)
meth_id15 <- substr(meth_clinical[[barcode_col]], 1, 15)
meth_clinical$pc1 <- meth_clinical_pc1$pc1[match(meth_id15, pc1_id15)]

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage) &
  !is.na(meth_clinical$pc1)
meth_c          <- meth_clinical[complete, ]
gene_mat_full_c <- gene_mat_full[, complete, drop = FALSE]
surv_c          <- surv_obj[complete]

cat("\nAnalysis cohort (with PC1 available):", sum(complete), "patients,",
    sum(meth_c$event), "events\n")


# - - - 3. Sign-bias test on the candidate list - - -

cat("\n=== Sign bias among the", nrow(candidates_df), "candidates ===\n")
n_neg <- sum(candidates_df$stratified_coef < 0)
n_tot <- nrow(candidates_df)
bt <- binom.test(n_neg, n_tot, p = 0.5)
cat(sprintf("%d of %d (%.1f%%) negative. binom.test p = %.3g\n",
            n_neg, n_tot, 100 * n_neg / n_tot, bt$p.value))
cat("A p-value this extreme means near-unanimous sign is essentially never",
    "expected by chance -- consistent with one shared axis driving most of",
    "these hits rather than independent gene-specific effects.\n")


# - - - 4. Per-candidate correlation with PC1 - - -

cat("\n=== Per-candidate correlation with PC1 ===\n")

pc1_z <- meth_c$pc1

cor_results <- t(sapply(candidates_df$gene, function(g) {
  if (!(g %in% rownames(gene_mat_full_c))) return(c(rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(gene_mat_full_c[g, ], pc1_z, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}))
candidates_df$pc1_cor_rho <- cor_results[, "rho"]
candidates_df$pc1_cor_p   <- cor_results[, "p"]

cat(sprintf("Median |rho| with PC1 across candidates: %.2f\n",
            median(abs(candidates_df$pc1_cor_rho), na.rm = TRUE)))
cat(sprintf("Candidates with |rho| > 0.5 (strongly PC1-linked): %d of %d\n",
            sum(abs(candidates_df$pc1_cor_rho) > 0.5, na.rm = TRUE), nrow(candidates_df)))


# - - - 5. Does each candidate survive adjusting for PC1, on top of subtype+purity? - - -
# Same logic as B7 (09_pc1_biological_characterization.R), applied per gene:
# if PC1 absorbs the candidate's apparent effect, its coefficient should
# shrink toward zero and/or lose significance once PC1 is in the model.

cat("\n=== Refitting each candidate with PC1 added as a covariate ===\n")

refit_with_pc1 <- function(g) {
  fit <- try(coxph(surv_c ~ scale(gene_mat_full_c[g, ]) + meth_c$age_years +
                      meth_c$stage + scale(meth_c$purity) + scale(meth_c$pc1) +
                      strata(meth_c$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), p = unname(s[5]))
}

pc1_adj_res <- t(sapply(candidates_df$gene, function(g) {
  if (!(g %in% rownames(gene_mat_full_c))) return(c(coef = NA, p = NA))
  refit_with_pc1(g)
}))
candidates_df$coef_pc1_adjusted <- pc1_adj_res[, "coef"]
candidates_df$p_pc1_adjusted    <- pc1_adj_res[, "p"]
candidates_df$survives_pc1_adjustment <- !is.na(candidates_df$p_pc1_adjusted) &
  candidates_df$p_pc1_adjusted < 0.05 &
  sign(candidates_df$coef_pc1_adjusted) == sign(candidates_df$stratified_coef)

n_survive <- sum(candidates_df$survives_pc1_adjustment, na.rm = TRUE)
cat(sprintf("\n%d of %d candidates (%.0f%%) survive PC1 adjustment",
            n_survive, nrow(candidates_df), 100 * n_survive / nrow(candidates_df)))
cat(" (same direction AND still nominally p<0.05 with PC1 in the model).\n")
cat("A LOW survival rate here is the decisive evidence that most of this",
    "candidate list is the residual PC1/global-methylation confound",
    "resurfacing at high multiplicity, not independent findings -- exactly",
    "what the overdispersed permutation null and the extreme sign bias",
    "both already suggested.\n")

candidates_df <- candidates_df[order(candidates_df$stratified_padj), ]
print(candidates_df[, c("gene", "stratified_coef", "stratified_padj",
                         "pc1_cor_rho", "coef_pc1_adjusted", "p_pc1_adjusted",
                         "survives_pc1_adjustment")])


# - - - 6. Save outputs - - -

write.csv(candidates_df, file.path(paths$tables, "full_bypass_candidates_pc1_diagnostics.csv"),
          row.names = FALSE)

pdf(file.path(paths$diagnostics, "candidate_pc1_correlation_histogram.pdf"), width = 7, height = 5)
hist(candidates_df$pc1_cor_rho, breaks = 20, col = "steelblue", border = "white",
     xlab = "Spearman rho with PC1", main = "Candidate promoters' correlation with PC1")
abline(v = 0, lty = 2, col = "grey40")
dev.off()

write_manifest(paths$run, list(
  script = "23_full_bypass_candidate_diagnostics.R",
  n_candidates = nrow(candidates_df),
  sign_bias_n_negative = n_neg,
  sign_bias_p = signif(bt$p.value, 4),
  median_abs_pc1_cor = round(median(abs(candidates_df$pc1_cor_rho), na.rm = TRUE), 3),
  n_survive_pc1_adjustment = n_survive
))

cat("\n=== Done. Results in", paths$run, "===\n")
