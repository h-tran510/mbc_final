# - - - 0. Set up - - -

library(survival)
library(glmnet)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(dplyr)
library(timeROC)
library(riskRegression)
library(prodlim)
library(survminer)
library(doParallel)
library(foreach)
library(readxl)

source("scripts/00_utils.R")

K          <- 4
n_repeats  <- 10
k_inner    <- 4
n_screen   <- 200   # promoters kept per fold in the supervised prescreen
run_batch_correction <- FALSE   # see Part 1b -- review the ComBat caution first

paths <- make_run_dir(tag = "v3_tumor_purity", base = "02_analysis")

# - - - 1. Load processed methylation data - - -

m_mat        <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

cat("M-value matrix:", nrow(m_mat), "probes x", ncol(m_mat), "samples\n")
cat("Clinical rows:", nrow(meth_clinical), "\n")
cat("Events (deaths):", sum(meth_clinical$event), "\n")

# Mandatory alignment gate -- fails loudly instead of letting a silent
# mismatch propagate (see the barcode-suffix issue diagnosed earlier)
match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)


# - - - 1a. Tumor purity (Aran et al. 2015 consensus purity estimate) - - -

cat("\n=== Purity integration ===\n")

purity_path <- "01_processed_data/aran_purity_tcga.xlsx"

if (file.exists(purity_path)) {
  # First 2 rows are a title/blank row before the real header
  purity_tab <- read_excel(purity_path, sheet = "Supp Data 1", skip = 2)
  purity_tab <- as.data.frame(purity_tab)
  names(purity_tab)[names(purity_tab) == "Sample ID"] <- "sample_barcode"
  purity_tab$CPE <- as.numeric(purity_tab$CPE)  # read_excel can import this as character
  purity_tab <- purity_tab[purity_tab$`Cancer type` == "BRCA" & !is.na(purity_tab$`Cancer type`), ]
  purity_tab$sample_id_15 <- substr(purity_tab$sample_barcode, 1, 15)  # drop aliquot suffix

  clin_id_15   <- substr(meth_clinical$sample_barcode, 1, 15)
  purity_match <- match(clin_id_15, purity_tab$sample_id_15)
  meth_clinical$purity <- purity_tab$CPE[purity_match]

  cat("Patients with purity estimate:",
      sum(!is.na(meth_clinical$purity)), "of", nrow(meth_clinical), "\n")

  if (sum(!is.na(meth_clinical$purity)) > 50) {
    cat("\nPurity vs survival:\n")
    print(summary(coxph(Surv(time, event) ~ scale(purity) + age_years + stage,
                         data = meth_clinical)))
  }
} else {
  cat("No purity file found at", purity_path, "-- proceeding without purity.\n")
  cat("Download Aran et al. 2015 (Nature Communications 6:8971) Supplementary",
      "Data 1 and save it there to enable purity adjustment.\n")
  meth_clinical$purity <- NA_real_
}

purity_available <- !any(is.na(meth_clinical$purity))
cat("Purity available for all samples:", purity_available, "\n")

# - - - 1b. Batch correction (optional -- off by default) - - -

cat("\n=== Batch correction ===\n")

if (run_batch_correction) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Install sva (BiocManager::install('sva')) to run batch correction.")
  }
  library(sva)

  plate <- sapply(strsplit(meth_clinical$sample_barcode, "-"), `[`, 6)

  if (!purity_available) {
    cat("Purity missing for some samples -- omitting purity from the ComBat",
        "protective model (keeping all samples) rather than dropping them.\n")
    mod_matrix <- model.matrix(~ age_years + stage + event, data = meth_clinical)
  } else {
    mod_matrix <- model.matrix(~ age_years + stage + event + purity,
                                data = meth_clinical)
  }
  stopifnot(nrow(mod_matrix) == ncol(m_mat))  # guards against silent row drop

  m_mat_corrected <- ComBat(dat = m_mat, batch = plate, mod = mod_matrix,
                             par.prior = TRUE)
  saveRDS(m_mat_corrected, file.path(paths$run, "meth_m_matrix_combat.rds"))

  m_mat <- m_mat_corrected   # downstream code uses m_mat transparently
  cat("Batch-corrected matrix in use for all downstream steps.\n")
} else {
  cat("Batch correction skipped (run_batch_correction = FALSE). Recall: plate",
      "does not predict survival directly (p = 0.43), only global",
      "methylation (p = 0.012) -- this is optional cleanup, not required to",
      "trust the survival results.\n")
}


# - - - 1c. Promoter aggregation - - -
# annotate probes, identify promoters, average M-values per gene.
# Only the annotation lookup is cached -- gene_mat itself depends on m_mat
# (which changes if batch correction is enabled above), so it is NOT cached
# across runs to avoid silently serving stale values from a different run.

ann_cache_path <- "01_processed_data/ann_matched_cache.rds"
ann_matched <- cache_or_compute(ann_cache_path, function() {
  ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann_probe_idx <- match(rownames(m_mat), rownames(ann_df))
  ann_df[ann_probe_idx, ]
})

is_promoter <- grepl("TSS200|TSS1500", ann_matched$UCSC_RefGene_Group)
m_mat_promoters <- m_mat[is_promoter, ]
gene_labels <- ann_matched$UCSC_RefGene_Name[is_promoter]

multi_labels <- strsplit(gene_labels, split = ";")
gene_labels_clean <- sapply(multi_labels, function(x) x[1])
gene_sums <- rowsum(m_mat_promoters, group = gene_labels_clean)

gene_counts <- table(gene_labels_clean)
gene_mat <- gene_sums / as.vector(gene_counts)   # ALL promoters, unfiltered

cat("\nTotal promoters (unfiltered):", nrow(gene_mat), "\n")


# - - - 2. Legacy design matrix (variance-prescreened top 500) - - -
# Kept as a fast, cheap reference point -- NOT the source of the final model.

probe_var     <- apply(gene_mat, 1, var)
n_prescreened <- min(500, nrow(gene_mat))
top_var_idx   <- order(probe_var, decreasing = TRUE)[seq_len(n_prescreened)]
m_lasso       <- gene_mat[top_var_idx, ]

clin_mat      <- model.matrix(~ age_years + stage, data = meth_clinical)[, -1, drop = FALSE]
probe_mat_t   <- t(m_lasso)
X             <- cbind(clin_mat, probe_mat_t)
clin_colnames <- colnames(clin_mat)
gene_col_idx  <- which(!colnames(X) %in% clin_colnames)

penalty_factors <- c(rep(0, ncol(clin_mat)), rep(1, ncol(probe_mat_t)))


# - - - 3. Legacy full-data LASSO Cox (fast reference, not cross-validated) - - -

set.seed(42)
cv_fit <- cv.glmnet(
  x = X, y = surv_obj, family = "cox", penalty.factor = penalty_factors,
  alpha = 1, nfolds = 10, standardize = TRUE
)

selected_1se <- extract_probes(cv_fit, "lambda.1se", clin_colnames)
selected_min <- extract_probes(cv_fit, "lambda.min", clin_colnames)
selected_probes <- if (length(selected_1se) > 0) selected_1se else selected_min
lambda_used     <- if (length(selected_1se) > 0) "lambda.1se" else "lambda.min"

cat("\nLegacy variance-filtered LASSO probes at", lambda_used, ":",
    length(selected_probes), "\n")
if (length(selected_probes) > 0) print(selected_probes)


# - - - 4. Post-LASSO correlation deduplication (legacy set, for reporting) - - -

if (length(selected_probes) > 1) {
  cor_mat  <- cor(t(gene_mat[selected_probes, , drop = FALSE]), method = "pearson")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hc       <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.3)

  lasso_coefs_raw <- as.matrix(coef(cv_fit, s = lambda_used))
  dedup_probes <- vapply(unique(clusters), function(cl) {
    probes_in_cl <- intersect(names(clusters)[clusters == cl],
                               rownames(lasso_coefs_raw))
    if (length(probes_in_cl) == 0) return(NA_character_)
    coef_vals <- abs(lasso_coefs_raw[probes_in_cl, 1, drop = TRUE])
    if (all(is.na(coef_vals))) return(NA_character_)
    probes_in_cl[which.max(coef_vals)]
  }, character(1))
  selected_probes <- unique(dedup_probes[!is.na(dedup_probes)])
}

cat("\nRepresentative legacy LASSO probe set:\n")
print(selected_probes)


# - - - 5. Parallelized, purity-aware supervised prescreen (nested CV) - - -
# Replaces the old sequential top-500-variance nested CV. Screens ALL
# promoters by univariate Cox p-value, re-selected inside each training fold
# to avoid leakage, across gene_mat (not the variance-filtered m_lasso).

screen_train <- function(gm, sv, clin, idx, n_keep, purity_available) {
  pv <- apply(gm[, idx, drop = FALSE], 1, function(g) {
    res <- if (purity_available) {
      uni_screen(g, sv[idx], clin$age_years[idx], clin$stage[idx],
                 purity = clin$purity[idx])
    } else {
      uni_screen(g, sv[idx], clin$age_years[idx], clin$stage[idx])
    }
    res["p"]
  })
  names(sort(pv))[seq_len(n_keep)]
}

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("\nRunning", n_repeats, "x", K, "-fold supervised screen across",
    n_cores, "cores...\n")

results <- foreach(
  r = seq_len(n_repeats),
  .packages = c("survival", "glmnet")
) %dopar% {

  set.seed(r)
  folds_ids  <- integer(nrow(meth_clinical))
  events_idx <- which(meth_clinical$event == 1)
  cen_idx    <- which(meth_clinical$event == 0)
  folds_ids[events_idx] <- sample(rep(1:K, length.out = length(events_idx)))
  folds_ids[cen_idx]    <- sample(rep(1:K, length.out = length(cen_idx)))

  cv_risk_rep   <- numeric(nrow(meth_clinical))
  fold_selected <- character(0)

  for (k in seq_len(K)) {
    tr_idx <- which(folds_ids != k)
    te_idx <- which(folds_ids == k)

    keep_genes <- screen_train(gene_mat, surv_obj, meth_clinical, tr_idx,
                                n_screen, purity_available)
    fold_selected <- c(fold_selected, keep_genes)

    clin_formula <- if (purity_available) ~ age_years + stage + purity
                                      else ~ age_years + stage
    clin_k <- model.matrix(clin_formula, data = meth_clinical)[, -1, drop = FALSE]
    stopifnot(nrow(clin_k) == nrow(meth_clinical))

    X_k  <- cbind(clin_k, t(gene_mat[keep_genes, , drop = FALSE]))
    pf_k <- c(rep(0, ncol(clin_k)), rep(1, length(keep_genes)))
    gene_cols_k <- (ncol(clin_k) + 1):ncol(X_k)

    fit_k <- cv.glmnet(X_k[tr_idx, ], surv_obj[tr_idx], family = "cox",
                        penalty.factor = pf_k, nfolds = k_inner,
                        standardize = TRUE)

    cf <- as.matrix(coef(fit_k, s = "lambda.1se"))[colnames(X_k)[gene_cols_k], 1]
    if (all(cf == 0)) {
      cf <- as.matrix(coef(fit_k, s = "lambda.min"))[colnames(X_k)[gene_cols_k], 1]
    }
    cv_risk_rep[te_idx] <- as.numeric(X_k[te_idx, gene_cols_k, drop = FALSE] %*% cf)
  }

  list(risk = cv_risk_rep, selected = fold_selected,
       c_index = concordance(surv_obj ~ cv_risk_rep, reverse = TRUE)$concordance)
}
stopCluster(cl)

all_risk_scores <- sapply(results, function(x) x$risk)
c_index_reps    <- sapply(results, function(x) x$c_index)
selected_freq   <- table(unlist(lapply(results, function(x) x$selected)))
mean_cv_risk    <- rowMeans(all_risk_scores)

cat("\n-- Per-repeat results (supervised screen) --\n")
for (r in seq_len(n_repeats)) {
  cat(sprintf("  Repeat %2d  C-index: %.3f\n", r, c_index_reps[r]))
}

cat("\nMost consistently screened-in promoters (of", n_repeats * K, "folds):\n")
print(head(sort(selected_freq, decreasing = TRUE), 25))

cat("\nOverlap between legacy variance-LASSO set and top supervised-screen genes:\n")
top_supervised <- names(sort(selected_freq, decreasing = TRUE))[seq_len(min(25, length(selected_freq)))]
print(intersect(selected_probes, top_supervised))


# - - - 6. C-index: stability summary and clinical comparison - - -

clinical_formula <- if (purity_available) {
  Surv(time, event) ~ age_years + stage + purity
} else {
  Surv(time, event) ~ age_years + stage
}
clinical_model <- coxph(clinical_formula, data = meth_clinical, x = TRUE, y = TRUE)
c_clinical     <- summary(clinical_model)$concordance["C"]
c_mean_score   <- concordance(surv_obj ~ mean_cv_risk, reverse = TRUE)$concordance

cat("\n-- C-index summary --\n")
cat(sprintf("  Clinical only:                        %.3f\n", c_clinical))
cat(sprintf("  CV methylation score (mean repeats):  %.3f\n", c_mean_score))
cat(sprintf("  Per-repeat mean (SD):                 %.3f (%.3f)\n",
            mean(c_index_reps), sd(c_index_reps)))
cat(sprintf("  95%% empirical CI:                     %.3f - %.3f\n",
            quantile(c_index_reps, 0.025), quantile(c_index_reps, 0.975)))


# - - - 7. Final Cox model on averaged CV methylation risk scores - - -

model_df <- cbind(
  meth_clinical[, c("time", "event", "age_years", "stage")],
  cv_meth_risk_score = mean_cv_risk
)
if (purity_available) model_df$purity <- meth_clinical$purity

final_formula <- if (purity_available) {
  Surv(time, event) ~ cv_meth_risk_score + age_years + stage + purity
} else {
  Surv(time, event) ~ cv_meth_risk_score + age_years + stage
}
final_model <- coxph(final_formula, data = model_df, x = TRUE, y = TRUE)

cat("\n-- Final Cox model --\n")
print(summary(final_model))

c_final <- summary(final_model)$concordance["C"]
cat(sprintf("\nFinal model C-index: %.3f (SE: %.3f)\n",
            c_final, summary(final_model)$concordance["se(C)"]))


# - - - 7b. IBS + time-dependent AUC - - -

timeroc_result <- timeROC(T = meth_clinical$time, delta = meth_clinical$event,
                           marker = mean_cv_risk, cause = 1,
                           weighting = "marginal",
                           times = c(365, 730, 1825, 3652),
                           ROC = TRUE, iid = FALSE)

brier_score <- Score(
  object = list("Methylation" = final_model, "Clinical" = clinical_model),
  formula = Surv(time, event) ~ 1, data = model_df,
  metrics = "brier", summary = "IPA", cause = 1,
  times = c(365, 730, 1825, 3652),
  plots = "calibration"
)

cat(sprintf("AUC at 1yr: %.3f\n", timeroc_result$AUC[1]))
print(brier_score$Brier$score)


# - - - 8. Kaplan-Meier split on median risk score - - -

risk_group <- factor(ifelse(mean_cv_risk >= median(mean_cv_risk), "High", "Low"),
                      levels = c("Low", "High"))
km_fit  <- survfit(surv_obj ~ risk_group)
lr_test <- survdiff(surv_obj ~ risk_group)
lr_pval <- 1 - pchisq(lr_test$chisq, df = 1)

pdf(file.path(paths$diagnostics, "methylation_risk_groups.pdf"), width = 7, height = 6)
plot(km_fit, col = c("steelblue", "tomato"), lty = 1, lwd = 2,
     xlab = "Time (days)", ylab = "Overall Survival Probability",
     main = "Kaplan-Meier: Methylation Risk Score Groups (Median Split)",
     mark.time = TRUE)
legend("topright",
       legend = c(paste0("Low risk  (n=", sum(risk_group == "Low"), ")"),
                  paste0("High risk (n=", sum(risk_group == "High"), ")")),
       col = c("steelblue", "tomato"), lty = 1, lwd = 2, bty = "n")
mtext(paste0("Log-rank p = ", format.pval(lr_pval, digits = 3)), side = 3, line = 0.3)
dev.off()

pdf(file.path(paths$diagnostics, "c_index_stability_methylation.pdf"), width = 6, height = 5)
hist(c_index_reps, breaks = 12, col = "steelblue", border = "white",
     xlab = "C-index (out-of-fold)",
     main = paste0("Methylation CV C-index Distribution (", n_repeats, " x ", K, "-fold)"),
     xlim = c(min(c_index_reps, c_clinical) - 0.03, max(c_index_reps, c_clinical) + 0.03))
abline(v = mean(c_index_reps), col = "tomato", lwd = 2, lty = 2)
abline(v = c_clinical, col = "gray40", lwd = 2, lty = 3)
legend("topright",
       legend = c(paste0("Methylation mean: ", round(mean(c_index_reps), 3)),
                  paste0("Clinical only:    ", round(c_clinical, 3))),
       col = c("tomato", "gray40"), lty = c(2, 3), lwd = 2, bty = "n")
dev.off()

# ROC curves -- fixed: first call starts fresh (no add=TRUE), times now
# correctly cover 365/730/1825/3652 (1yr/2yr/5yr/10yr) with no duplicate
pdf(file.path(paths$diagnostics, "timeroc_curves.pdf"), width = 6, height = 6)
plot(timeroc_result, time = 365,  col = "blue",     title = TRUE, lwd = 2)
plot(timeroc_result, time = 730,  col = "red",       add = TRUE, lwd = 2)
plot(timeroc_result, time = 1825, col = "darkgreen", add = TRUE, lwd = 2)
plot(timeroc_result, time = 3652, col = "purple",    add = TRUE, lwd = 2)
legend("bottomright", legend = c("1yr", "2yr", "5yr", "10yr"),
       col = c("blue", "red", "darkgreen", "purple"), lwd = 2, bty = "n")
dev.off()

# Calibration plot -- fixed: "Clinical only" -> "Clinical" to match the name
# actually registered in Score() above
pdf(file.path(paths$diagnostics, "calibration_plot.pdf"), width = 6, height = 6)
plotCalibration(x = brier_score, models = c("Methylation", "Clinical"),
                 times = 3652, method = "quantile", q = 5)
dev.off()

pdf(file.path(paths$diagnostics, "auc_over_time.pdf"), width = 6, height = 5)
plot(x = c(365, 730, 1825, 3652), y = timeroc_result$AUC, type = "b", pch = 19,
     col = "steelblue", lwd = 2, xlab = "Time (days)", ylab = "AUC(t)",
     ylim = c(0.4, 1), main = "Time-dependent AUC: Methylation Risk Score")
abline(h = 0.5, lty = 2, col = "gray50")
dev.off()

pdf(file.path(paths$diagnostics, "forest_plot_final_model.pdf"), width = 7, height = 5)
print(ggforest(final_model, data = model_df))
dev.off()

cat("\nPlots saved to", paths$diagnostics, "\n")


# - - - 9. Save outputs - - -

params <- c(K = K, n_repeats = n_repeats, K_inner = k_inner, n_screen = n_screen,
            endpoint = "PFI", seed = 42,
            purity_used = purity_available, batch_corrected = run_batch_correction)
write_manifest(paths$run, params)

saveRDS(cv_fit, file.path(paths$models, "legacy_lasso_cv_fit.rds"))
saveRDS(final_model, file.path(paths$models, "methylation_final_model.rds"))

write.csv(data.frame(gene = names(selected_freq), n_folds = as.integer(selected_freq)),
          file.path(paths$tables, "supervised_screen_selection_frequency.csv"),
          row.names = FALSE)

write.csv(data.frame(probe = selected_probes),
          file.path(paths$tables, "legacy_variance_lasso_probes.csv"),
          row.names = FALSE)

write.csv(data.frame(sample_barcode = meth_clinical$sample_barcode,
                      patient_id = meth_clinical$patient_id,
                      mean_cv_risk_score = mean_cv_risk,
                      risk_group = as.character(risk_group)),
          file.path(paths$tables, "methylation_cv_risk_scores.csv"),
          row.names = FALSE)

write.csv(data.frame(repeat_id = seq_len(n_repeats), c_index = c_index_reps),
          file.path(paths$tables, "methylation_c_index_per_repeat.csv"),
          row.names = FALSE)

cat("\nOutputs saved to", paths$run, "\n")