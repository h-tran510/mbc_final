# - - - 0. Set up - - -

# modeling packages
library(survival)
library(glmnet)

# data related packages
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(dplyr)

# metric packages
library(timeROC)
library(riskRegression)
library(prodlim)
library(survminer)

# call helper functions
source("scripts/00_utils.R")

K <- 4
n_repeats <- 10
k_inner <- 4

# call helper function make_run_dir with tag and base
paths <- make_run_dir(tag = "null_and_directional", base = "02_analysis")

# write manifest about session
params <- c(K = K, n_repeats = n_repeats, K_inner = k_inner,
            prescreen = 500, endpoint = "PFI", seed = 42)
write_manifest(paths$run, params)

# - - - 1. Load processed methylation data - - -
m_mat        <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

cat("M-value matrix:", nrow(m_mat), "probes x", ncol(m_mat), "samples\n")
cat("Clinical rows:", nrow(meth_clinical), "\n")
cat("Events (deaths):", sum(meth_clinical$event), "\n")

# Truncate IDs
meth_sample_id <- substr(colnames(m_mat), 1, 15)
clin_sample_id <- substr(meth_clinical$sample_barcode, 1, 15)

# Build match
match_idx <- match(clin_sample_id, meth_sample_id)

# Reorder m_mat's columns
m_mat <- m_mat[, match_idx]

# safety check in case missing matches in the future
stopifnot(all(!is.na(match_idx)))

stopifnot(identical(substr(colnames(m_mat), 1, 15),
                     substr(meth_clinical$sample_barcode, 1, 15)))

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

# - - - 1b. Promoter aggregation - - -
# annotate probes, identify promoters, average beta values

# load the annotation object
ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# filter m_mat to promoter associated probes
ann_probe_idx <- match(rownames(m_mat), rownames(ann_df))

ann_matched <- ann_df[ann_probe_idx, ]
is_promoter <- grepl("TSS200|TSS1500", ann_matched$UCSC_RefGene_Group)

m_mat_promoters <- m_mat[is_promoter, ]
gene_labels <- ann_matched$UCSC_RefGene_Name[is_promoter]

# note: approximately 9096 promoter-associated probes

# aggregate (mean) probes per gene and produces gene_mat
# sum(grepl(";", gene_labels)) # note: 3532 labels split with a semicolon

multi_labels <- strsplit(gene_labels, split = ";")
gene_labels_clean <- sapply(multi_labels, function(x) x[1])
gene_sums <- rowsum(m_mat_promoters, group = gene_labels_clean)

gene_counts <- table(gene_labels_clean)
gene_mat <- gene_sums / as.vector(gene_counts)

# - - - 2. Build design matrix - - -
# Pre-screen to top 5000 probes by variance
# Clinical covariates placed first and marked unpenalized (penalty.factor = 0)

#probe_var is now using gene_mat
probe_var     <- apply(gene_mat, 1, var)
n_prescreened <- min(500, nrow(gene_mat))
top_var_idx   <- order(probe_var, decreasing = TRUE)[seq_len(n_prescreened)]
m_lasso       <- gene_mat[top_var_idx, ]

cat("Promoters entering LASSO:", nrow(m_lasso), "\n")

clin_mat      <- model.matrix(~ age_years + stage, data = meth_clinical)[, -1, drop = FALSE]
probe_mat_t <- t(m_lasso)
X             <- cbind(clin_mat, probe_mat_t)
clin_colnames <- colnames(clin_mat)
gene_col_idx <- which(!colnames(X) %in% clin_colnames)

penalty_factors <- c(rep(0, ncol(clin_mat)), rep(1, ncol(probe_mat_t)))

# - - - 2b. Running tests - - -
 
# confirm if event count is enough
cat("\n D2. Endpoint / power \n")
 
n_events <- sum(meth_clinical$event)
n_total  <- nrow(meth_clinical)
 
cat("N samples:", n_total, "\n")
cat("N events:", n_events, sprintf("(%.1f%%)\n", 100 * n_events / n_total))
cat("Median follow-up (days), censored only:",
    median(meth_clinical$time[meth_clinical$event == 0]), "\n")
cat("Median time to event (days):",
    median(meth_clinical$time[meth_clinical$event == 1]), "\n")
 
# Rule of thumb: a Cox model supports roughly n_events / 10 free parameters.
cat("Supportable covariates (events/10):", floor(n_events / 10), "\n")
cat("Features entering LASSO:", ncol(X), "\n")
cat("Ratio (features per event):", round(ncol(X) / n_events, 1), "\n")
 
# The clinical model already uses some of that budget
cat("Clinical df used:", length(clin_colnames), "\n")
 
# Follow-up distribution -- administrative censoring at study end vs
# genuine loss to follow-up behave very differently
hist(meth_clinical$time[meth_clinical$event == 0],
     breaks = 40, main = "Censoring time distribution", xlab = "Days")

# run univariate model and exmaine p-value distribution
cat("\n D4. Univariate Cox screen (all promoters) \n")
 
uni_screen <- function(g) {
  fit <- try(coxph(surv_obj ~ scale(g) + meth_clinical$age_years + meth_clinical$stage),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}
 
uni_res <- as.data.frame(t(apply(gene_mat, 1, uni_screen)))
uni_res$gene <- rownames(gene_mat)
uni_res$p_adj <- p.adjust(uni_res$p, method = "BH")
 
cat("Promoters tested:", nrow(uni_res), "\n")
cat("Raw p < 0.05:", sum(uni_res$p < 0.05, na.rm = TRUE),
    sprintf("(expect ~%.0f by chance)\n", 0.05 * nrow(uni_res)))
cat("BH-adjusted q < 0.10:", sum(uni_res$p_adj < 0.10, na.rm = TRUE), "\n")
cat("BH-adjusted q < 0.05:", sum(uni_res$p_adj < 0.05, na.rm = TRUE), "\n")
 
hist(uni_res$p, breaks = 50, main = "Univariate Cox p-values (all promoters)",
     xlab = "p-value")
abline(h = nrow(uni_res) / 50, col = "red", lty = 2)  # uniform expectation

pdf(file.path(paths$diagnostics, "d4_univariate_p_histogram.pdf"), width = 6, height = 5)
hist(uni_res$p, breaks = 50, main = "Univariate Cox p-values (all promoters)",
     xlab = "p-value")
abline(h = nrow(uni_res) / 50, col = "red", lty = 2)
dev.off()

top_uni <- head(uni_res[order(uni_res$p), ], 30)
print(top_uni[, c("gene", "coef", "p", "p_adj")])
 
cat("\nOverlap between top-500 by VARIANCE and top-500 by univariate p:\n")
top500_var <- rownames(gene_mat)[order(apply(gene_mat, 1, var), decreasing = TRUE)][1:500]
top500_cox <- uni_res$gene[order(uni_res$p)][1:500]
cat("  Overlap:", length(intersect(top500_var, top500_cox)), "of 500\n")
cat("  (Expected by chance:", round(500 * 500 / nrow(gene_mat)), ")\n")
 
# see if promoter aggregation is working
cat("\n D6. Aggregation audit \n")
 
cat("\n=== D6. Aggregation audit ===\n")
 
gene_sums_check   <- rowsum(m_mat_promoters, group = gene_labels_clean)
gene_counts_check <- table(gene_labels_clean)
cat("Names identical:",
    identical(rownames(gene_sums_check), names(gene_counts_check)), "\n")
cat("NA gene labels:", sum(is.na(gene_labels_clean)), "\n")
 
grp_split  <- strsplit(ann_matched$UCSC_RefGene_Group[is_promoter], ";")
name_split <- strsplit(gene_labels, ";")
first_is_promoter <- mapply(function(g) {
  length(g) > 0 && grepl("TSS200|TSS1500", g[1])
}, grp_split)
cat("Probes where position 1 is NOT the TSS entry:",
    sum(!first_is_promoter), "of", length(first_is_promoter), "\n")
 
tss200_idx  <- grepl("TSS200",  ann_matched$UCSC_RefGene_Group[is_promoter])
tss1500_idx <- grepl("TSS1500", ann_matched$UCSC_RefGene_Group[is_promoter])
genes_with_both <- intersect(
  unique(gene_labels_clean[tss200_idx]),
  unique(gene_labels_clean[tss1500_idx])
)
cat("Genes with both TSS200 and TSS1500 probes:", length(genes_with_both), "\n")
 
within_gene_cor <- sapply(head(genes_with_both, 300), function(g) {
  a <- m_mat_promoters[tss200_idx  & gene_labels_clean == g, , drop = FALSE]
  b <- m_mat_promoters[tss1500_idx & gene_labels_clean == g, , drop = FALSE]
  if (nrow(a) == 0 || nrow(b) == 0) return(NA)
  cor(colMeans(a), colMeans(b))
})
cat("Median TSS200-vs-TSS1500 correlation within gene:",
    round(median(within_gene_cor, na.rm = TRUE), 3), "\n")
cat("Proportion NEGATIVELY correlated:",
    round(mean(within_gene_cor < 0, na.rm = TRUE), 3), "\n")
 
pdf(file.path(paths$diagnostics, "d6_tss200_vs_tss1500_correlation.pdf"),
    width = 6, height = 5)
hist(within_gene_cor, breaks = 40,
     main = "Within-gene correlation: TSS200 vs TSS1500",
     xlab = "Pearson r")
abline(v = 0, col = "red", lty = 2)
dev.off()
 
set.seed(42)
probe_subset <- sample(nrow(m_mat_promoters), min(2000, nrow(m_mat_promoters)))
probe_p <- apply(m_mat_promoters[probe_subset, ], 1, function(g) {
  fit <- try(coxph(surv_obj ~ scale(g) + meth_clinical$age_years + meth_clinical$stage),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NA)
  unname(summary(fit)$coefficients[1, 5])
})
cat("Best gene-level p:", signif(min(uni_res$p, na.rm = TRUE), 3), "\n")
cat("Best probe-level p (2000 sampled):", signif(min(probe_p, na.rm = TRUE), 3), "\n")
 
# see if variance filter is selecting anything prognostically insignificant
cat("\n D7. Composition of the variance-filtered set \n")
 
gene_chr <- tapply(ann_matched$chr[is_promoter], gene_labels_clean, function(x) x[1])
 
cat("Chromosome distribution, top 500 by variance:\n")
print(sort(table(gene_chr[top500_var]), decreasing = TRUE)[1:10])
cat("Sex-chromosome genes in top 500 by variance:",
    sum(gene_chr[top500_var] %in% c("chrX", "chrY"), na.rm = TRUE), "\n")
cat("Sex-chromosome genes in top 500 by univariate p:",
    sum(gene_chr[top500_cox] %in% c("chrX", "chrY"), na.rm = TRUE), "\n")
 
pdf(file.path(paths$diagnostics, "d7_variance_vs_pvalue.pdf"),
    width = 6, height = 5)
plot(log10(apply(gene_mat, 1, var)),
     -log10(uni_res[rownames(gene_mat), "p"]),
     pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.3),
     xlab = "log10(variance)", ylab = "-log10(univariate Cox p)",
     main = "Does variance predict outcome association?")
abline(v = log10(sort(apply(gene_mat, 1, var), decreasing = TRUE)[500]),
       col = "red", lty = 2)
legend("topright", "top-500 variance cutoff", lty = 2, col = "red", bty = "n")
dev.off()

cat("\nAll diagnostic plots saved to:", paths$diagnostics, "\n")

# D10. Directional bias in the significant hits
# -----------------------------------------------------------------------------
# Under independent biology, prognostic promoters should split roughly evenly
# between protective (negative coef) and risk-increasing (positive coef).
# A strong skew means something global is pushing many promoters the same way.
 
cat("\n=== D10. Directional bias ===\n")
 
sig <- uni_res[!is.na(uni_res$p) & uni_res$p < 0.05, ]
n_neg <- sum(sig$coef < 0)
n_pos <- sum(sig$coef > 0)
 
cat("Significant at p < 0.05:", nrow(sig), "\n")
cat("  Negative coef (higher methylation -> LOWER hazard):", n_neg, "\n")
cat("  Positive coef (higher methylation -> HIGHER hazard):", n_pos, "\n")
cat(sprintf("  Proportion negative: %.3f\n", n_neg / nrow(sig)))
 
# Formal test: is the sign split different from 50/50?
print(binom.test(n_neg, nrow(sig), p = 0.5))
 
# Same check on ALL promoters, not just significant ones -- a global shift
# should show up as a shifted coefficient distribution overall
cat("\nMean coefficient across ALL promoters:",
    round(mean(uni_res$coef, na.rm = TRUE), 4), "\n")
cat("(Expect ~0 if no global factor)\n")
 
pdf(file.path(paths$diagnostics, "d10_coefficient_sign_distribution.pdf"),
    width = 6, height = 5)
hist(uni_res$coef, breaks = 60,
     main = "Univariate Cox coefficients (all promoters)",
     xlab = "coef (log hazard ratio per SD methylation)")
abline(v = 0, col = "red", lty = 2, lwd = 2)
dev.off()
 
 
# -----------------------------------------------------------------------------
# D11. Estimate what fraction of promoters are truly non-null
# -----------------------------------------------------------------------------
# Storey-style: the flat right-hand portion of the p-value histogram is almost
# entirely null promoters. Its height tells you what proportion (pi0) are null.
# If the implied non-null fraction is implausibly large for independent
# biology, a shared factor is the parsimonious explanation.
 
cat("\n=== D11. Non-null fraction ===\n")
 
lambda_cut <- 0.5
pi0 <- mean(uni_res$p > lambda_cut, na.rm = TRUE) / (1 - lambda_cut)
pi0 <- min(pi0, 1)
 
cat(sprintf("Estimated pi0 (proportion null):     %.3f\n", pi0))
cat(sprintf("Estimated non-null proportion:       %.3f\n", 1 - pi0))
cat(sprintf("Implied number of non-null promoters: %.0f of %d\n",
            (1 - pi0) * nrow(uni_res), nrow(uni_res)))
 
# Interpretation: a few dozen genuinely prognostic promoters is biologically
# ordinary. Many hundreds, all leaning the same direction, is the signature
# of a sample-level covariate rather than promoter-specific biology.
 
# - - -  3. Full-data LASSO Cox (for probe selection and deduplication) - - -
set.seed(42)
cv_fit <- cv.glmnet(
  x              = X,
  y              = surv_obj,
  family         = "cox",
  penalty.factor = penalty_factors,
  alpha          = 1,
  nfolds         = 10,
  standardize    = TRUE
)

cat("\nOptimal lambda (1se):", cv_fit$lambda.1se, "\n")
cat("Optimal lambda (min):", cv_fit$lambda.min, "\n")

extract_probes <- function(fit, s, clin_colnames) {
  coefs        <- as.matrix(coef(fit, s = s))
  selected_all <- rownames(coefs)[coefs[, 1] != 0]
  selected_all[!selected_all %in% clin_colnames]
}

selected_1se <- extract_probes(cv_fit, "lambda.1se", clin_colnames)
selected_min <- extract_probes(cv_fit, "lambda.min", clin_colnames)

cat("\nLASSO-selected probes at lambda.1se:", length(selected_1se), "\n")
if (length(selected_1se) > 0) print(selected_1se)

cat("\nLASSO-selected probes at lambda.min:", length(selected_min), "\n")
if (length(selected_min) > 0) print(head(selected_min, 20))

selected_probes <- if (length(selected_1se) > 0) selected_1se else selected_min
lambda_used     <- if (length(selected_1se) > 0) "lambda.1se" else "lambda.min"

# -----------------------------------------------------------------------------
# D12. Find the shared factor
# -----------------------------------------------------------------------------
# Three candidates, in order of likelihood for TCGA methylation:
#   1. Global methylation level (mean M-value per sample)
#   2. Tumor purity / normal-cell contamination
#   3. Batch (plate / processing date, from the barcode)
 
cat("\n=== D12. Candidate shared factors ===\n")
 
# --- D12a. Global mean methylation per sample ---
global_meth <- colMeans(gene_mat)
 
cat("Global methylation vs survival:\n")
print(summary(coxph(surv_obj ~ scale(global_meth) +
                      meth_clinical$age_years + meth_clinical$stage)))
 
# How much of each promoter's association is explained by this one number?
cat("\nCorrelation of each promoter with global methylation:\n")
cor_with_global <- apply(gene_mat, 1, function(g) cor(g, global_meth))
cat("  Median |r|:", round(median(abs(cor_with_global)), 3), "\n")
cat("  Proportion with |r| > 0.5:", round(mean(abs(cor_with_global) > 0.5), 3), "\n")
 
# Are the SIGNIFICANT promoters more correlated with global methylation
# than the non-significant ones? If yes, that's the tell.
sig_genes <- sig$gene
cat("\n  Median |r| among p<0.05 promoters:",
    round(median(abs(cor_with_global[sig_genes]), na.rm = TRUE), 3), "\n")
cat("  Median |r| among p>=0.05 promoters:",
    round(median(abs(cor_with_global[setdiff(uni_res$gene, sig_genes)]),
                 na.rm = TRUE), 3), "\n")
 
pdf(file.path(paths$diagnostics, "d12_global_methylation_vs_significance.pdf"),
    width = 6, height = 5)
plot(abs(cor_with_global[uni_res$gene]), -log10(uni_res$p),
     pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.3),
     xlab = "|correlation with global methylation|",
     ylab = "-log10(univariate Cox p)",
     main = "Is significance driven by a global factor?")
dev.off()
 
# --- D12b. Principal components of the promoter matrix ---
# PC1 captures the dominant shared axis of variation. If PC1 is strongly
# prognostic and the top univariate hits all load heavily on it, you have
# one factor, not hundreds of independent signals.
pca <- prcomp(t(gene_mat), scale. = TRUE)
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
 
cat("\nVariance explained by first 5 PCs:\n")
print(round(var_explained[1:5], 3))
 
cat("\nPC1-PC5 association with survival (adjusted for age + stage):\n")
for (i in 1:5) {
  f <- coxph(surv_obj ~ scale(pca$x[, i]) +
               meth_clinical$age_years + meth_clinical$stage)
  s <- summary(f)$coefficients[1, ]
  cat(sprintf("  PC%d: coef = %+.3f, p = %.3g (explains %.1f%% of variance)\n",
              i, unname(s[1]), unname(s[5]), 100 * var_explained[i]))
}
 
cat("\nCorrelation of PC1 with global methylation:",
    round(cor(pca$x[, 1], global_meth), 3), "\n")
 
pdf(file.path(paths$diagnostics, "d12_pca_variance_explained.pdf"),
    width = 6, height = 5)
barplot(var_explained[1:20], names.arg = 1:20,
        xlab = "PC", ylab = "Proportion of variance explained",
        main = "Scree plot: promoter methylation")
dev.off()
 
# --- D12c. Batch / plate effect from the TCGA barcode ---
# Barcode field 6 is the plate ID. If plate is prognostic, that is a
# processing artifact, not biology.
plate <- sapply(strsplit(meth_clinical$sample_barcode, "-"), `[`, 6)
cat("\nDistinct plates:", length(unique(plate)), "\n")
print(head(sort(table(plate), decreasing = TRUE), 10))
 
if (length(unique(plate)) > 1 && length(unique(plate)) < nrow(meth_clinical) / 5) {
  cat("\nPlate association with survival:\n")
  print(anova(coxph(surv_obj ~ factor(plate))))
 
  cat("\nPlate association with global methylation (ANOVA):\n")
  print(summary(aov(global_meth ~ factor(plate))))
}
 
# --- D12d. Tumor purity, if you have it ---
# TCGA purity estimates (ABSOLUTE / ESTIMATE) are not in the standard clinical
# file. If you have a purity column, uncomment:
# cat("\nPurity vs survival:\n")
# print(summary(coxph(surv_obj ~ scale(meth_clinical$purity) +
#                       meth_clinical$age_years + meth_clinical$stage)))
# cat("Purity vs global methylation r:",
#     round(cor(meth_clinical$purity, global_meth, use = "complete.obs"), 3), "\n")

# - - - 3b. Elastic net alpha search - - -
 
# alpha_grid <- seq(0, 1, by = 0.1)
# alpha_cvm  <- numeric(length(alpha_grid))
# alpha_fits <- vector("list", length(alpha_grid))
 
# cat("\nSearching alpha for elastic net...\n")
 
# for (i in seq_along(alpha_grid)) {
#   set.seed(42)
#   fit_i <- cv.glmnet(
#     x              = X,
#     y              = surv_obj,
#     family         = "cox",
#     penalty.factor = penalty_factors,
#     alpha          = alpha_grid[i],
#     nfolds         = 10,
#     standardize    = TRUE
#   )
#   alpha_fits[[i]] <- fit_i
#   alpha_cvm[i]    <- min(fit_i$cvm) 
 
#   cat(sprintf("  alpha = %.1f   cvm(lambda.min) = %.4f\n", alpha_grid[i], alpha_cvm[i]))
# }
 
# best_idx    <- which.min(alpha_cvm)
# best_alpha  <- alpha_grid[best_idx]
# enet_cv_fit <- alpha_fits[[best_idx]]
 
# cat(sprintf("\nBest alpha by CV deviance: %.1f\n", best_alpha))
# cat("(alpha = 1 is pure LASSO -- already fit in step 3 as cv_fit; alpha = 0 is ridge)\n")
 
# selected_enet_1se <- extract_probes(enet_cv_fit, "lambda.1se", clin_colnames)
# selected_enet_min <- extract_probes(enet_cv_fit, "lambda.min", clin_colnames)
 
# cat("\nElastic net-selected probes at lambda.1se:", length(selected_enet_1se), "\n")
# if (length(selected_enet_1se) > 0) print(selected_enet_1se)
 
# cat("\nProbes in elastic net set but dropped by LASSO (alpha = 1):\n")
# print(setdiff(selected_enet_1se, selected_1se))
 
# cat("\nProbes kept by LASSO but dropped by elastic net:\n")
# print(setdiff(selected_1se, selected_enet_1se))
 
# # Save the alpha search itself for reference/reporting
# write.csv(
#   data.frame(alpha = alpha_grid, cvm_lambda_min = alpha_cvm),
#   file.path(paths$tables, "methylation_enet_alpha_search.csv"),
#   row.names = FALSE
# )

# - - - 4. Post-LASSO correlation deduplication (for probe reporting) - - -
# Clusters selected probes at |r| > 0.7 and retains the one with the largest
# |LASSO coefficient| per cluster. Used for biological interpretation only —
# the risk score in step 5 uses all per-fold coefficients without deduplication.

if (length(selected_probes) > 1) {
  cor_mat  <- cor(t(gene_mat[selected_probes, , drop = FALSE]), method = "pearson")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hc       <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.3)

  lasso_coefs_raw <- as.matrix(coef(cv_fit, s = lambda_used))

  dedup_probes <- vapply(unique(clusters), function(cl) {
    probes_in_cl <- names(clusters)[clusters == cl]
    probes_in_cl <- intersect(probes_in_cl, rownames(lasso_coefs_raw))
    if (length(probes_in_cl) == 0) return(NA_character_)
    coef_vals <- abs(lasso_coefs_raw[probes_in_cl, 1, drop = TRUE])
    if (all(is.na(coef_vals))) return(NA_character_)
    probes_in_cl[which.max(coef_vals)]
  }, character(1))

  selected_probes <- unique(dedup_probes[!is.na(dedup_probes)])

  cat(
    "\nRepresentative probes after correlation deduplication (|r| > 0.7):",
    length(selected_probes), "\n"
  )
}

cat("\nRepresentative LASSO probe set:\n")
print(selected_probes)

# - - 5. Repeated K-fold CV for stable risk scores and C-index distribution - -
# 20 x 10-fold nested CV. In each outer fold, LASSO is refitted on the training
# partition (inner 5-fold CV for lambda), and the probe-only linear predictor is
# scored on the held-out fold — no data leakage into the final model.
# Averaging risk scores across 20 repeats reduces fold-assignment variance.

all_risk_scores <- matrix(NA_real_, nrow = nrow(meth_clinical), ncol = n_repeats)
c_index_reps    <- numeric(n_repeats)

cat("\nRunning", n_repeats, "x", K, "-fold CV...\n")

for (r in seq_len(n_repeats)) {
  set.seed(r)
  folds_ids <- integer(nrow(meth_clinical))

  # Build stratified outer fold
  events_idx <- which(meth_clinical$event == 1)
  cen_idx <- which(meth_clinical$event == 0)
  folds_ids[events_idx] <- sample(rep(1:K, length.out = length(events_idx)))
  folds_ids[cen_idx] <- sample(rep(1:K, length.out = length(cen_idx)))

  cv_risk_rep <- numeric(nrow(meth_clinical))

  for (k in seq_len(K)) {
    train_idx <- folds_ids != k
    test_idx  <- folds_ids == k

    # Build stratified inner fold
    events_in <- which(meth_clinical$event[train_idx] == 1)
    cen_in <- which(meth_clinical$event[train_idx] == 0)
    inner_foldid <- numeric(sum(train_idx))
    inner_foldid[events_in] <- sample(rep(1:k_inner, length.out = length(events_in)))
    inner_foldid[cen_in] <- sample(rep(1:k_inner, length.out = length(cen_in)))

    inner_cv <- cv.glmnet(
      x              = X[train_idx, ],
      y              = surv_obj[train_idx],
      family         = "cox",
      penalty.factor = penalty_factors,
      foldid = inner_foldid,
      standardize    = TRUE
    )

    inner_coefs   <- as.matrix(coef(inner_cv, s = "lambda.1se"))
    probe_coefs_k <- inner_coefs[colnames(X)[gene_col_idx], 1]

    if (all(probe_coefs_k == 0)) {
      inner_coefs   <- as.matrix(coef(inner_cv, s = "lambda.min"))
      probe_coefs_k <- inner_coefs[colnames(X)[gene_col_idx], 1]
    }

    cv_risk_rep[test_idx] <- as.numeric(
      X[test_idx, gene_col_idx, drop = FALSE] %*% probe_coefs_k
    )
  }

  all_risk_scores[, r] <- cv_risk_rep
  c_index_reps[r] <- concordance(surv_obj ~ cv_risk_rep, reverse = TRUE)$concordance

  cat(sprintf("  Repeat %2d  C-index: %.3f\n", r, c_index_reps[r]))
}

mean_cv_risk <- rowMeans(all_risk_scores)

# - - - 5b. Direction sanity checks for methylation risk score - - -
# Cox/glmnet linear predictors are usually "risk" scores:
# higher score = higher hazard = shorter survival.
# Therefore, the correct C-index may require reverse = TRUE
# or using the negative risk score.

# cat("\n── Direction sanity checks ──\n")

# c_raw <- concordance(surv_obj ~ mean_cv_risk)$concordance
# c_rev <- concordance(surv_obj ~ mean_cv_risk, reverse = TRUE)$concordance
# c_neg <- concordance(surv_obj ~ I(-mean_cv_risk))$concordance

# cat(sprintf("C-index, raw mean_cv_risk:          %.3f\n", c_raw))
# cat(sprintf("C-index, reverse = TRUE:            %.3f\n", c_rev))
# cat(sprintf("C-index, using -mean_cv_risk:       %.3f\n", c_neg))

# # Higher risk score should correspond to higher event rate / shorter survival.
# risk_group_check <- factor(
#   ifelse(mean_cv_risk >= median(mean_cv_risk), "High score", "Low score"),
#   levels = c("Low score", "High score")
# )

# cat("\nEvent counts by risk-score group:\n")
# print(table(risk_group_check, meth_clinical$event))

# cat("\nEvent rates by risk-score group:\n")
# print(prop.table(table(risk_group_check, meth_clinical$event), margin = 1))

# km_check <- survdiff(surv_obj ~ risk_group_check)
# km_check_p <- 1 - pchisq(km_check$chisq, df = 1)

# cat(sprintf("\nLog-rank p-value for direction check: %.4f\n", km_check_p))

# cox_check <- coxph(surv_obj ~ mean_cv_risk)
# cat("\nUnivariable Cox model for methylation risk score:\n")
# print(summary(cox_check))

# - - - 6. C-index: stability summary and clinical comparison - - -
clinical_model <- coxph(Surv(time, event) ~ age_years + stage, data = meth_clinical)
c_clinical     <- summary(clinical_model)$concordance["C"]
c_mean_score <- concordance(surv_obj ~ mean_cv_risk, reverse = TRUE)$concordance

cat("\n── C-index summary ──\n")
cat(sprintf("  Clinical only (age + stage):             %.3f\n", c_clinical))
cat(sprintf("  CV methylation score (mean repeats):     %.3f\n", c_mean_score))
cat(sprintf("  Per-repeat mean  (SD):                   %.3f  (%.3f)\n",
            mean(c_index_reps), sd(c_index_reps)))
cat(sprintf("  Per-repeat range:                        %.3f - %.3f\n",
            min(c_index_reps), max(c_index_reps)))
cat(sprintf("  95%% empirical CI:                        %.3f - %.3f\n",
            quantile(c_index_reps, 0.025), quantile(c_index_reps, 0.975)))

# - - - 7. Final Cox model on averaged CV methylation risk scores - - -
model_df <- cbind(
  meth_clinical[, c("time", "event", "age_years", "stage")],
  cv_meth_risk_score = mean_cv_risk
)

final_model <- coxph(
  Surv(time, event) ~ cv_meth_risk_score + age_years + stage,
  data = model_df
)

cat("\n── Final Cox model (CV methylation risk score + clinical covariates) ──\n")
print(summary(final_model))

c_final <- summary(final_model)$concordance["C"]
cat(sprintf("\nFinal model C-index: %.3f  (SE: %.3f)\n",
            c_final, summary(final_model)$concordance["se(C)"]))

# - - - 7b. Determine IBS + time-dependent AUC - - -

# model with timeROC()
timeroc_result <- timeROC(T = meth_clinical$time,
  delta = meth_clinical$event,
  marker = mean_cv_risk,
  cause = 1,
  weighting = "marginal",
  times = c(365, 730, 1825, 3652),
  ROC = TRUE,
  iid = FALSE
)

# model integrated brier score
brier_score <- Score(
  object = list("Methylation" = final_model, "Clinical" = clinical_model),
  formula = Surv(time, event) ~ 1,
  data = model_df,
  metrics = "brier",
  summary = "IPA",
  cause = 1
)

cat(sprintf("AUC at 1yr: %.3f\n", timeroc_result$AUC[1]))
print(brier_score)

# - - - 8. Kaplan-Meier split on median risk score - - -
risk_group <- factor(
  ifelse(mean_cv_risk >= median(mean_cv_risk), "High", "Low"),
  levels = c("Low", "High")
)

km_fit  <- survfit(surv_obj ~ risk_group)
lr_test <- survdiff(surv_obj ~ risk_group)
lr_pval <- 1 - pchisq(lr_test$chisq, df = 1)

cat(sprintf("\nLog-rank test p-value (High vs Low risk): %.4f\n", lr_pval))
cat("Group sizes:\n"); print(table(risk_group))

pdf(file.path(paths$diagnostics, "methylation_risk_groups.pdf"),
    width = 7, height = 6)
plot(
  km_fit,
  col       = c("steelblue", "tomato"),
  lty       = 1,
  lwd       = 2,
  xlab      = "Time (days)",
  ylab      = "Overall Survival Probability",
  main      = "Kaplan-Meier: Methylation Risk Score Groups (Median Split)",
  mark.time = TRUE
)
legend(
  "topright",
  legend = c(
    paste0("Low risk  (n=", sum(risk_group == "Low"),  ")"),
    paste0("High risk (n=", sum(risk_group == "High"), ")")
  ),
  col = c("steelblue", "tomato"),
  lty = 1, lwd = 2, bty = "n"
)
mtext(paste0("Log-rank p = ", format.pval(lr_pval, digits = 3)), side = 3, line = 0.3)
dev.off()

pdf(file.path(paths$diagnostics, "c_index_stability_methylation.pdf"),
    width = 6, height = 5)
hist(
  c_index_reps,
  breaks  = 12,
  col     = "steelblue",
  border  = "white",
  xlab    = "C-index (out-of-fold)",
  main    = paste0("Methylation CV C-index Distribution (", n_repeats, " x ", K, "-fold)"),
  xlim    = c(
    min(c_index_reps, c_clinical) - 0.03,
    max(c_index_reps, c_clinical) + 0.03
  )
)
abline(v = mean(c_index_reps), col = "tomato", lwd = 2, lty = 2)
abline(v = c_clinical,         col = "gray40", lwd = 2, lty = 3)
legend(
  "topright",
  legend = c(
    paste0("Methylation mean: ", round(mean(c_index_reps), 3)),
    paste0("Clinical only:    ", round(c_clinical, 3))
  ),
  col = c("tomato", "gray40"), lty = c(2, 3), lwd = 2, bty = "n"
)
dev.off()

cat("\nPlots saved to diagnostics/\n")

# plot ROC curves
plot(timeroc_result, time = 365,  col = "blue",      title = TRUE, lwd = 2)  # 1yr, no add
plot(timeroc_result, time = 730,  col = "red",        add = TRUE, lwd = 2)   # 2yr
plot(timeroc_result, time = 1825, col = "darkgreen",  add = TRUE, lwd = 2)   # 5yr
plot(timeroc_result, time = 3652, col = "purple",     add = TRUE, lwd = 2)   # 10yr

legend("bottomright",
       legend = c("1yr", "2yr", "5yr", "10yr"),
       col = c("blue", "red", "darkgreen", "purple"), lwd = 2, bty = "n")

# plot calibration plot
plotCalibration(
  x = brier_score,
  models = c("Methylation", "Clinical"),
  times = 3652,
  method = "quantile",
  q = 5
)

# plot AUC(t), which is AUC over time
pdf(file.path(paths$diagnostics, "auc_over_time.pdf"), width = 6, height = 5)
plot(
  x = c(365, 730, 1825, 3652),
  y = timeroc_result$AUC,
  type = "b", pch = 19, col = "steelblue", lwd = 2,
  xlab = "Time (days)", ylab = "AUC(t)",
  ylim = c(0.4, 1),
  main = "Time-dependent AUC: Methylation Risk Score"
)
abline(h = 0.5, lty = 2, col = "gray50")
dev.off()

# forest plot of cox hazard values
pdf(file.path(paths$diagnostics, "forest_plot_final_model.pdf"), width = 7, height = 5)
print(ggforest(final_model, data = model_df))
dev.off()

# - - - 9. Save outputs - - -
# save models
saveRDS(cv_fit, file.path(paths$models, "methylation_lasso_cv_fit.rds"))
saveRDS(final_model, file.path(paths$models, "methylation_final_model.rds"))

# save diagnostics (selected genes, risk scores, and c index per repeat)
write.csv(
  data.frame(probe = selected_probes),
  file.path(paths$tables, "methylation_selected_probes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_barcode     = meth_clinical$sample_barcode,
    patient_id         = meth_clinical$patient_id,
    mean_cv_risk_score = mean_cv_risk,
    risk_group         = as.character(risk_group)
  ),
  file.path(paths$tables, "methylation_cv_risk_scores.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(repeat_id = seq_len(n_repeats), c_index = c_index_reps),
  file.path(paths$tables, "methylation_c_index_per_repeat.csv"),
  row.names = FALSE
)

cat("\nOutputs saved to 02_analysis/\n")
