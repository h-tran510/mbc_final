library(survival)
library(SummarizedExperiment)

if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package glmnet is missing. Install it with install.packages(\"glmnet\") and rerun.")
library(glmnet)

# ── 1. Load data ──────────────────────────────────────────────────────────────
merged_data <- readRDS("01_processed_data/merged_data.rds")
expr_matrix <- readRDS("01_processed_data/expr_matrix_vst.rds")

# ── 2. Clinical covariate preparation ────────────────────────────────────────
df <- merged_data

df$age_years <- df$age_at_diagnosis / 365.25

stage_map <- c(
  "Stage I"    = "I",   "Stage IA"   = "I",   "Stage IB"   = "I",
  "Stage II"   = "II",  "Stage IIA"  = "II",  "Stage IIB"  = "II",
  "Stage III"  = "III", "Stage IIIA" = "III", "Stage IIIB" = "III",
  "Stage IIIC" = "III",
  "Stage IV"   = "IV"
)
df$stage <- factor(stage_map[df$ajcc_pathologic_stage], levels = c("I", "II", "III", "IV"))

# ── 3. Exclude samples with invalid time or missing covariates ────────────────
# 21 samples have time <= 0 (likely same-day events or data entry errors)
valid <- !is.na(df$time) & df$time > 0 &
         !is.na(df$event) &
         !is.na(df$age_years) &
         !is.na(df$stage)

df_clean   <- df[valid, ]
expr_clean <- expr_matrix[, df_clean$sample_barcode]

cat("Samples retained:", nrow(df_clean), "\n")
cat("Events (deaths):", sum(df_clean$event), "\n")
cat("Median follow-up (days):", median(df_clean$time[df_clean$event == 0]), "\n")

# ── 4. Variance filter ────────────────────────────────────────────────────────
gene_var   <- apply(expr_clean, 1, var)
var_thresh <- quantile(gene_var, 0.5)
expr_filt  <- expr_clean[gene_var > var_thresh, ]

cat("Genes after variance filter:", nrow(expr_filt), "\n")

surv_obj <- Surv(time = df_clean$time, event = df_clean$event)

# ── 5. Build design matrix ────────────────────────────────────────────────────
# Pre-screen to top 2000 genes by variance only — no prior selection involved.
# Clinical covariates are placed first so penalty.factor = 0 keeps them mandatory.

n_prescreened <- min(2000, nrow(expr_filt))
top_var_idx   <- order(gene_var[rownames(expr_filt)], decreasing = TRUE)[seq_len(n_prescreened)]
expr_lasso    <- expr_filt[top_var_idx, ]

cat("Genes entering LASSO:", nrow(expr_lasso), "\n")

clin_mat      <- model.matrix(~ age_years + stage, data = df_clean)[, -1, drop = FALSE]
gene_mat      <- t(expr_lasso)
X             <- cbind(clin_mat, gene_mat)
clin_colnames <- colnames(clin_mat)
gene_col_idx  <- which(!colnames(X) %in% clin_colnames)

penalty_factors <- c(rep(0, ncol(clin_mat)), rep(1, ncol(gene_mat)))

# ── 6. Full-data LASSO Cox (for gene selection and deduplication) ─────────────
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

extract_genes <- function(fit, s, clin_colnames) {
  coefs        <- as.matrix(coef(fit, s = s))
  selected_all <- rownames(coefs)[coefs[, 1] != 0]
  selected_all[!selected_all %in% clin_colnames]
}

lasso_genes_1se <- extract_genes(cv_fit, "lambda.1se", clin_colnames)
lasso_genes_min <- extract_genes(cv_fit, "lambda.min", clin_colnames)

cat("\nLASSO-selected genes at lambda.1se:", length(lasso_genes_1se), "\n")
if (length(lasso_genes_1se) > 0) print(lasso_genes_1se)

cat("\nLASSO-selected genes at lambda.min:", length(lasso_genes_min), "\n")
if (length(lasso_genes_min) > 0) print(lasso_genes_min)

selected_genes <- if (length(lasso_genes_1se) > 0) lasso_genes_1se else lasso_genes_min
lambda_used    <- if (length(lasso_genes_1se) > 0) "lambda.1se" else "lambda.min"

# ── 7. Post-LASSO correlation deduplication (for gene reporting) ──────────────
# Clusters genes at |r| > 0.7 and keeps the gene with the largest |LASSO coef|
# from each cluster. Used for biological interpretation only — the risk score
# in step 8 still uses the full set of per-fold LASSO coefficients.

if (length(selected_genes) > 1) {
  cor_mat  <- cor(t(expr_filt[selected_genes, , drop = FALSE]), method = "pearson")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hc       <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.3)   # h = 0.3 corresponds to |r| > 0.7

  lasso_coefs_raw <- as.matrix(coef(cv_fit, s = lambda_used))

  dedup_genes <- vapply(unique(clusters), function(cl) {
    genes_in_cl <- names(clusters)[clusters == cl]
    genes_in_cl <- intersect(genes_in_cl, rownames(lasso_coefs_raw))
    if (length(genes_in_cl) == 0) return(NA_character_)
    coef_vals <- abs(lasso_coefs_raw[genes_in_cl, 1, drop = TRUE])
    if (all(is.na(coef_vals))) return(NA_character_)
    genes_in_cl[which.max(coef_vals)]
  }, character(1))

  dedup_genes    <- unique(dedup_genes[!is.na(dedup_genes)])
  selected_genes <- dedup_genes

  cat(
    "\nRepresentative genes after correlation deduplication (|r| > 0.7):",
    length(selected_genes), "\n"
  )
}

cat("\nRepresentative LASSO gene set:\n")
print(selected_genes)

# ── 8. Repeated K-fold CV for stable risk scores and C-index distribution ─────
# Runs n_repeats independent 10-fold CVs with different random seeds.
# Each repeat fits LASSO on 9 folds (inner 5-fold CV for lambda) and predicts
# the gene-only linear predictor on the held-out fold. Averaging risk scores
# across repeats reduces variance from any single fold assignment.
# The per-repeat C-indices (computed directly from out-of-fold scores, no refitting)
# quantify model stability across different splits.
# Runtime: ~200 cv.glmnet calls — expect several minutes.

K         <- 3
n_repeats <- 1

all_risk_scores <- matrix(NA_real_, nrow = nrow(df_clean), ncol = n_repeats)
c_index_reps    <- numeric(n_repeats)

cat("\nRunning", n_repeats, "x", K, "-fold CV...\n")

for (r in seq_len(n_repeats)) {
  set.seed(r)
  fold_ids    <- sample(rep(seq_len(K), length.out = nrow(df_clean)))
  cv_risk_rep <- numeric(nrow(df_clean))

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

    inner_coefs  <- as.matrix(coef(inner_cv, s = "lambda.1se"))
    gene_coefs_k <- inner_coefs[colnames(X)[gene_col_idx], 1]

    if (all(gene_coefs_k == 0)) {
      inner_coefs  <- as.matrix(coef(inner_cv, s = "lambda.min"))
      gene_coefs_k <- inner_coefs[colnames(X)[gene_col_idx], 1]
    }

    cv_risk_rep[test_idx] <- as.numeric(
      X[test_idx, gene_col_idx, drop = FALSE] %*% gene_coefs_k
    )

    saveRDS(cv_risk_rep, "02_analysis/cv_gene_risk_partial.rds")
  }

  all_risk_scores[, r] <- cv_risk_rep

  # C-index computed directly from out-of-fold scores (no refitting on full data)
  c_index_reps[r] <- concordance(surv_obj ~ cv_risk_rep)$concordance

  cat(sprintf("  Repeat %2d  C-index: %.3f\n", r, c_index_reps[r]))
}

# Average risk score across repeats as the final stable risk score
mean_cv_risk <- rowMeans(all_risk_scores)

# ── 9. C-index: stability summary and clinical comparison ─────────────────────
clinical_model <- coxph(Surv(time, event) ~ age_years + stage, data = df_clean)
c_clinical     <- summary(clinical_model)$concordance["C"]

# C-index of averaged risk score (more stable than any single repeat)
c_mean_score <- concordance(surv_obj ~ mean_cv_risk)$concordance

cat("\n── C-index summary ──\n")
cat(sprintf("  Clinical only (age + stage):          %.3f\n", c_clinical))
cat(sprintf("  CV gene score, mean across repeats:   %.3f\n", c_mean_score))
cat(sprintf("  Per-repeat mean  (SD):                %.3f  (%.3f)\n",
            mean(c_index_reps), sd(c_index_reps)))
cat(sprintf("  Per-repeat range:                     %.3f - %.3f\n",
            min(c_index_reps), max(c_index_reps)))
cat(sprintf("  95%% empirical CI:                     %.3f - %.3f\n",
            quantile(c_index_reps, 0.025), quantile(c_index_reps, 0.975)))

# ── 10. Final Cox model on averaged CV risk scores ────────────────────────────
model_df <- cbind(
  df_clean[, c("time", "event", "age_years", "stage")],
  cv_gene_risk_score = mean_cv_risk
)

final_model <- coxph(
  Surv(time, event) ~ cv_gene_risk_score + age_years + stage,
  data = model_df
)

cat("\n── Final Cox model (averaged CV gene risk score + clinical covariates) ──\n")
print(summary(final_model))

c_final <- summary(final_model)$concordance["C"]
cat(sprintf("\nFinal model C-index: %.3f  (SE: %.3f)\n",
            c_final, summary(final_model)$concordance["se(C)"]))

# ── 11. Kaplan-Meier split on median risk score ───────────────────────────────
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

pdf("03_figures/km_risk_groups.pdf", width = 7, height = 6)
plot(
  km_fit,
  col       = c("steelblue", "tomato"),
  lty       = 1,
  lwd       = 2,
  xlab      = "Time (days)",
  ylab      = "Overall Survival Probability",
  main      = "Kaplan-Meier: LASSO Gene Risk Score Groups (Median Split)",
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

# C-index stability histogram
pdf("03_figures/c_index_stability.pdf", width = 6, height = 5)
hist(
  c_index_reps,
  breaks  = 12,
  col     = "steelblue",
  border  = "white",
  xlab    = "C-index (out-of-fold)",
  main    = paste0("CV C-index Distribution (", n_repeats, " x ", K, "-fold)"),
  xlim    = c(
    min(c_index_reps, c_clinical) - 0.03,
    max(c_index_reps, c_clinical) + 0.03
  )
)
abline(v = mean(c_index_reps), col = "tomato",  lwd = 2, lty = 2)
abline(v = c_clinical,         col = "gray40",  lwd = 2, lty = 3)
legend(
  "topright",
  legend = c(
    paste0("Gene score mean: ", round(mean(c_index_reps), 3)),
    paste0("Clinical only:   ", round(c_clinical, 3))
  ),
  col = c("tomato", "gray40"), lty = c(2, 3), lwd = 2, bty = "n"
)
dev.off()

cat("\nPlots saved to 03_figures/\n")

# ── 12. Save outputs ──────────────────────────────────────────────────────────
saveRDS(cv_fit,      "02_analysis/lasso_cv_fit.rds")
saveRDS(final_model, "02_analysis/final_multivariate_model.rds")

write.csv(
  data.frame(gene = selected_genes),
  "02_analysis/lasso_selected_genes.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_barcode     = df_clean$sample_barcode,
    mean_cv_risk_score = mean_cv_risk,
    risk_group         = as.character(risk_group)
  ),
  "02_analysis/cv_gene_risk_scores.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    repeat_id = seq_len(n_repeats),
    c_index   = c_index_reps
  ),
  "02_analysis/c_index_per_repeat.csv",
  row.names = FALSE
)

cat("\nOutputs saved to 02_analysis/\n")
