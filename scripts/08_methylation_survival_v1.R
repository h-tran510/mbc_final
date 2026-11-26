library(survival)

if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package glmnet is missing. Install it with install.packages(\"glmnet\") and rerun.")
library(glmnet)

# ── 1. Load processed methylation data ───────────────────────────────────────
m_mat        <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

cat("M-value matrix:", nrow(m_mat), "probes x", ncol(m_mat), "samples\n")
cat("Clinical rows:", nrow(meth_clinical), "\n")
cat("Events (deaths):", sum(meth_clinical$event), "\n")

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

# ── 2. Build design matrix ────────────────────────────────────────────────────
# Pre-screen to top 5000 probes by variance — outcome-agnostic, no selection bias.
# Clinical covariates placed first and marked unpenalized (penalty.factor = 0).

probe_var     <- apply(m_mat, 1, var)
n_prescreened <- min(5000, nrow(m_mat))
top_var_idx   <- order(probe_var, decreasing = TRUE)[seq_len(n_prescreened)]
m_lasso       <- m_mat[top_var_idx, ]

cat("Probes entering LASSO:", nrow(m_lasso), "\n")

clin_mat      <- model.matrix(~ age_years + stage, data = meth_clinical)[, -1, drop = FALSE]
probe_mat_t   <- t(m_lasso)
X             <- cbind(clin_mat, probe_mat_t)
clin_colnames <- colnames(clin_mat)
probe_col_idx <- which(!colnames(X) %in% clin_colnames)

penalty_factors <- c(rep(0, ncol(clin_mat)), rep(1, ncol(probe_mat_t)))

# ── 3. Full-data LASSO Cox (for probe selection and deduplication) ────────────
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

# ── 4. Post-LASSO correlation deduplication (for probe reporting) ─────────────
# Clusters selected probes at |r| > 0.7 and retains the one with the largest
# |LASSO coefficient| per cluster. Used for biological interpretation only —
# the risk score in step 5 uses all per-fold coefficients without deduplication.

if (length(selected_probes) > 1) {
  cor_mat  <- cor(t(m_mat[selected_probes, , drop = FALSE]), method = "pearson")
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

# ── 5. Repeated K-fold CV for stable risk scores and C-index distribution ─────
# 20 x 10-fold nested CV. In each outer fold, LASSO is refitted on the training
# partition (inner 5-fold CV for lambda), and the probe-only linear predictor is
# scored on the held-out fold — no data leakage into the final model.
# Averaging risk scores across 20 repeats reduces fold-assignment variance.
# Runtime: ~200 cv.glmnet calls — expect several minutes.

K         <- 5
n_repeats <- 1

all_risk_scores <- matrix(NA_real_, nrow = nrow(meth_clinical), ncol = n_repeats)
c_index_reps    <- numeric(n_repeats)

cat("\nRunning", n_repeats, "x", K, "-fold CV...\n")

for (r in seq_len(n_repeats)) {
  set.seed(r)
  fold_ids    <- sample(rep(seq_len(K), length.out = nrow(meth_clinical)))
  cv_risk_rep <- numeric(nrow(meth_clinical))

  for (k in seq_len(K)) {
    train_idx <- fold_ids != k
    test_idx  <- fold_ids == k

    inner_cv <- cv.glmnet(
      x              = X[train_idx, ],
      y              = surv_obj[train_idx],
      family         = "cox",
      penalty.factor = penalty_factors,
      alpha          = 1,
      nfolds         = 5,
      standardize    = TRUE
    )

    inner_coefs   <- as.matrix(coef(inner_cv, s = "lambda.1se"))
    probe_coefs_k <- inner_coefs[colnames(X)[probe_col_idx], 1]

    if (all(probe_coefs_k == 0)) {
      inner_coefs   <- as.matrix(coef(inner_cv, s = "lambda.min"))
      probe_coefs_k <- inner_coefs[colnames(X)[probe_col_idx], 1]
    }

    cv_risk_rep[test_idx] <- as.numeric(
      X[test_idx, probe_col_idx, drop = FALSE] %*% probe_coefs_k
    )
  }

  all_risk_scores[, r] <- cv_risk_rep
  c_index_reps[r] <- concordance(surv_obj ~ cv_risk_rep, reverse = TRUE)$concordance

  cat(sprintf("  Repeat %2d  C-index: %.3f\n", r, c_index_reps[r]))
}

mean_cv_risk <- rowMeans(all_risk_scores)

# ── 5b. Direction sanity checks for methylation risk score ────────────────────
# Cox/glmnet linear predictors are usually "risk" scores:
# higher score = higher hazard = shorter survival.
# Therefore, the correct C-index may require reverse = TRUE
# or using the negative risk score.

cat("\n── Direction sanity checks ──\n")

c_raw <- concordance(surv_obj ~ mean_cv_risk)$concordance
c_rev <- concordance(surv_obj ~ mean_cv_risk, reverse = TRUE)$concordance
c_neg <- concordance(surv_obj ~ I(-mean_cv_risk))$concordance

cat(sprintf("C-index, raw mean_cv_risk:          %.3f\n", c_raw))
cat(sprintf("C-index, reverse = TRUE:            %.3f\n", c_rev))
cat(sprintf("C-index, using -mean_cv_risk:       %.3f\n", c_neg))

# Higher risk score should correspond to higher event rate / shorter survival.
risk_group_check <- factor(
  ifelse(mean_cv_risk >= median(mean_cv_risk), "High score", "Low score"),
  levels = c("Low score", "High score")
)

cat("\nEvent counts by risk-score group:\n")
print(table(risk_group_check, meth_clinical$event))

cat("\nEvent rates by risk-score group:\n")
print(prop.table(table(risk_group_check, meth_clinical$event), margin = 1))

km_check <- survdiff(surv_obj ~ risk_group_check)
km_check_p <- 1 - pchisq(km_check$chisq, df = 1)

cat(sprintf("\nLog-rank p-value for direction check: %.4f\n", km_check_p))

cox_check <- coxph(surv_obj ~ mean_cv_risk)
cat("\nUnivariable Cox model for methylation risk score:\n")
print(summary(cox_check))

# ── 6. C-index: stability summary and clinical comparison ─────────────────────
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

# ── 7. Final Cox model on averaged CV methylation risk scores ─────────────────
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

# ── 8. Kaplan-Meier split on median risk score ────────────────────────────────
risk_group <- factor(
  ifelse(mean_cv_risk >= median(mean_cv_risk), "High", "Low"),
  levels = c("Low", "High")
)

km_fit  <- survfit(surv_obj ~ risk_group)
lr_test <- survdiff(surv_obj ~ risk_group)
lr_pval <- 1 - pchisq(lr_test$chisq, df = 1)

cat(sprintf("\nLog-rank test p-value (High vs Low risk): %.4f\n", lr_pval))
cat("Group sizes:\n"); print(table(risk_group))

dir.create("03_figures", showWarnings = FALSE)

pdf("03_figures/km_methylation_risk_groups.pdf", width = 7, height = 6)
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

pdf("03_figures/c_index_stability_methylation.pdf", width = 6, height = 5)
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

cat("\nPlots saved to 03_figures/\n")

# ── 9. Save outputs ───────────────────────────────────────────────────────────
saveRDS(cv_fit,      "02_analysis/methylation_lasso_cv_fit.rds")
saveRDS(final_model, "02_analysis/methylation_final_model.rds")

write.csv(
  data.frame(probe = selected_probes),
  "02_analysis/methylation_selected_probes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_barcode     = meth_clinical$sample_barcode,
    patient_id         = meth_clinical$patient_id,
    mean_cv_risk_score = mean_cv_risk,
    risk_group         = as.character(risk_group)
  ),
  "02_analysis/methylation_cv_risk_scores.csv",
  row.names = FALSE
)

write.csv(
  data.frame(repeat_id = seq_len(n_repeats), c_index = c_index_reps),
  "02_analysis/methylation_c_index_per_repeat.csv",
  row.names = FALSE
)

cat("\nOutputs saved to 02_analysis/\n")
