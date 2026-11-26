library(survival)
library(SummarizedExperiment)

if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package glmnet is missing: install the pinned environment (repository root README.md, section Install) and rerun.")
library(glmnet)
library(doParallel)
library(foreach)

source("deployed_scripts_v1/00_utils.R")

merged_data <- readRDS("01_processed_data/merged_data.rds")
expr_matrix <- readRDS("01_processed_data/expr_matrix_vst.rds")
assert_one_per_patient(merged_data$sample_barcode, "merged_data")
assert_one_per_patient(colnames(expr_matrix), "VST expression matrix columns")

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

valid <- !is.na(df$time) & df$time > 0 &
         !is.na(df$event) &
         !is.na(df$age_years) &
         !is.na(df$stage)

df_clean   <- df[valid, ]
expr_clean <- expr_matrix[, df_clean$sample_barcode]

cat("Samples retained:", nrow(df_clean), "\n")
cat("Events (deaths):", sum(df_clean$event), "\n")
cat("Median follow-up (days):", median(df_clean$time[df_clean$event == 0]), "\n")

gene_var   <- row_vars(expr_clean)
var_thresh <- quantile(gene_var, 0.5)
expr_filt  <- expr_clean[gene_var > var_thresh, ]

cat("Genes after variance filter:", nrow(expr_filt), "\n")

surv_obj <- Surv(time = df_clean$time, event = df_clean$event)

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

lasso_genes_1se <- extract_probes(cv_fit, "lambda.1se", clin_colnames)
lasso_genes_min <- extract_probes(cv_fit, "lambda.min", clin_colnames)

cat("\nLASSO-selected genes at lambda.1se:", length(lasso_genes_1se), "\n")
if (length(lasso_genes_1se) > 0) print(lasso_genes_1se)

cat("\nLASSO-selected genes at lambda.min:", length(lasso_genes_min), "\n")
if (length(lasso_genes_min) > 0) print(lasso_genes_min)

if (length(lasso_genes_1se) > 0) {
  selected_genes <- lasso_genes_1se
  lambda_used <- "lambda.1se"
} else {
  selected_genes <- lasso_genes_min
  lambda_used <- "lambda.min"
}

if (length(selected_genes) > 1) {
  cor_mat  <- cor(t(expr_filt[selected_genes, , drop = FALSE]), method = "pearson")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hc       <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.3)

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

K         <- full_or_test(10, 3)
n_repeats <- full_or_test(10, 1)

run_repeated_cv <- function(X, surv_obj, penalty_factors, gene_col_idx, K,
                            n_repeats, n_cores, worker_log_dir) {
  if (n_cores > 1) {
    cl <- make_worker_cluster(n_cores, worker_log_dir)
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }

  foreach(r = seq_len(n_repeats), .packages = c("survival", "glmnet")) %dopar% {
    set.seed(r)
    fold_ids    <- sample(rep(seq_len(K), length.out = nrow(X)))
    cv_risk_rep <- numeric(nrow(X))

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
    }

    list(risk    = cv_risk_rep,
         c_index = concordance(surv_obj ~ cv_risk_rep, reverse = TRUE)$concordance)
  }
}

n_cores <- min(n_repeats, n_workers(test_cap = 4))
cat("\nRunning", n_repeats, "x", K, "-fold CV on", n_cores, "worker(s)...\n")

cv_results <- run_repeated_cv(X, surv_obj, penalty_factors, gene_col_idx, K,
                              n_repeats, n_cores, "02_analysis")
all_risk_scores <- vapply(cv_results, function(x) x$risk, numeric(nrow(df_clean)))
c_index_reps    <- vapply(cv_results, function(x) x$c_index, numeric(1))
for (r in seq_len(n_repeats)) {
  cat(sprintf("  Repeat %2d  C-index: %.3f\n", r, c_index_reps[r]))
}

mean_cv_risk <- rowMeans(all_risk_scores)

clinical_model <- coxph(Surv(time, event) ~ age_years + stage, data = df_clean)
c_clinical     <- summary(clinical_model)$concordance["C"]

c_mean_score <- concordance(surv_obj ~ mean_cv_risk, reverse = TRUE)$concordance

cat("\n── C-index summary ──\n")
cat(sprintf("  Clinical only (age + stage):          %.3f\n", c_clinical))
cat(sprintf("  CV gene score, mean across repeats:   %.3f\n", c_mean_score))
if (n_repeats > 1) {
  cat(sprintf("  Per-repeat mean  (SD):                %.3f  (%.3f)\n",
              mean(c_index_reps), sd(c_index_reps)))
  cat(sprintf("  Per-repeat range:                     %.3f - %.3f\n",
              min(c_index_reps), max(c_index_reps)))
  cat(sprintf("  95%% empirical CI:                     %.3f - %.3f\n",
              quantile(c_index_reps, 0.025), quantile(c_index_reps, 0.975)))
} else {
  cat(sprintf("  Per-repeat mean  (SD):                %.3f  (n/a (1 repeat))\n",
              c_index_reps[1]))
  cat("  Per-repeat range:                     n/a (1 repeat)\n")
  cat("  95% empirical CI:                     n/a (1 repeat)\n")
}

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

risk_group <- factor(
  ifelse(mean_cv_risk >= median(mean_cv_risk), "High", "Low"),
  levels = c("Low", "High")
)

km_fit  <- survfit(surv_obj ~ risk_group)
lr_test <- survdiff(surv_obj ~ risk_group)
lr_pval <- pchisq(lr_test$chisq, df = length(lr_test$n) - 1, lower.tail = FALSE)

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

pdf("03_figures/c_index_stability.pdf", width = 6, height = 5)
if (n_repeats > 1) {
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
} else {
  plot.new()
  title(main = sprintf("CV C-index distribution: n/a (1 repeat)\nC-index %.3f, clinical only %.3f",
                       c_index_reps[1], c_clinical))
  cat("C-index stability histogram: n/a (1 repeat)\n")
}
dev.off()

cat("\nPlots saved to 03_figures/\n")

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
