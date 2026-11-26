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

source("deployed_scripts_v1/00_utils.R")

K          <- full_or_test(4, 3)
n_repeats  <- full_or_test(10, 2)
k_inner    <- full_or_test(4, 3)
n_screen   <- full_or_test(200, 50)
run_batch_correction <- FALSE
use_subtype          <- TRUE

paths <- make_run_dir(tag = "v3_tumor_purity", base = "02_analysis")

m_mat        <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

cat("M-value matrix:", nrow(m_mat), "probes x", ncol(m_mat), "samples\n")
cat("Clinical rows:", nrow(meth_clinical), "\n")
cat("Events (PFI):", sum(meth_clinical$event), "\n")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

cat("\n=== Purity integration ===\n")

sp <- load_subtype_purity(meth_clinical, need_subtype = use_subtype,
                          need_purity = FALSE)
meth_clinical$purity <- sp$purity
if (use_subtype) meth_clinical$subtype <- sp$subtype

if (attr(sp, "purity_loaded")) {
  cat("Patients with purity estimate:",
      sum(!is.na(meth_clinical$purity)), "of", nrow(meth_clinical), "\n")

  if (sum(!is.na(meth_clinical$purity)) > 50) {
    cat("\nPurity vs survival:\n")
    print(summary(coxph(Surv(time, event) ~ scale(purity) + age_years + stage,
                         data = meth_clinical)))
  }
} else {
  cat("Proceeding without purity. Download Aran et al. 2015 (Nature",
      "Communications 6:8971) Supplementary Data 1 to",
      "01_processed_data/aran_purity_tcga.xlsx to enable purity adjustment.\n")
}
rm(sp)

purity_missing <- sum(is.na(meth_clinical$purity))
if (purity_missing > 0 && purity_missing < nrow(meth_clinical)) {
  warning(sprintf(
    "%d of %d samples lack a purity estimate; those rows are dropped from purity-adjusted models by complete-case deletion.",
    purity_missing, nrow(meth_clinical)), call. = FALSE)
}
purity_available <- mean(!is.na(meth_clinical$purity)) > 0.90

cat("Purity available for", nrow(meth_clinical) - purity_missing, "of",
    nrow(meth_clinical), "samples; purity used in models:", purity_available, "\n")

if (purity_available && purity_missing > 0) {
  keep_rows <- !is.na(meth_clinical$purity)
  cat("Restricting the cross-validated screen to", sum(keep_rows),
      "purity-complete samples so glmnet and coxph share one cohort.\n")
  meth_clinical <- meth_clinical[keep_rows, ]
  m_mat         <- m_mat[, keep_rows, drop = FALSE]
  surv_obj      <- Surv(time = meth_clinical$time, event = meth_clinical$event)
}

if (use_subtype) {
  cat("\n=== PAM50 subtype in the clinical baseline ===\n")
  keep_rows <- !is.na(meth_clinical$subtype)
  cat("Dropping", sum(!keep_rows), "of", nrow(meth_clinical),
      "samples without a PAM50 call; subtype-complete cohort:", sum(keep_rows), "\n")
  meth_clinical <- meth_clinical[keep_rows, ]
  m_mat         <- m_mat[, keep_rows, drop = FALSE]
  surv_obj      <- Surv(time = meth_clinical$time, event = meth_clinical$event)
  meth_clinical$subtype <- droplevels(meth_clinical$subtype)
  print(table(subtype = meth_clinical$subtype, event = meth_clinical$event))
}

clin_terms <- c("age_years", "stage",
                if (purity_available) "purity",
                if (use_subtype) "subtype")
cat("Clinical baseline covariates:", paste(clin_terms, collapse = " + "), "\n")

cat("\n=== Batch correction ===\n")

if (run_batch_correction) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("sva is missing: install the pinned environment (repository root README.md, section Install) to run batch correction.")
  }
  library(sva)

  plate <- vapply(strsplit(meth_clinical$meth_sample_barcode, "-"),
                  function(x) if (length(x) >= 6) x[6] else NA_character_,
                  character(1))
  stopifnot("Could not parse a plate ID from every barcode" = !anyNA(plate))

  if (!purity_available) {
    cat("Purity missing for some samples -- omitting purity from the ComBat",
        "protective model (keeping all samples) rather than dropping them.\n")
  }
  mod_terms  <- c("age_years", "stage", "event",
                  if (purity_available) "purity",
                  if (use_subtype) "subtype")
  cat("ComBat protected covariates:", paste(mod_terms, collapse = " + "), "\n")
  mod_matrix <- model.matrix(reformulate(mod_terms), data = meth_clinical)
  stopifnot(nrow(mod_matrix) == ncol(m_mat))

  m_mat_corrected <- ComBat(dat = m_mat, batch = plate, mod = mod_matrix,
                             par.prior = TRUE)
  saveRDS(m_mat_corrected, file.path(paths$run, "meth_m_matrix_combat.rds"))

  m_mat <- m_mat_corrected
  cat("Batch-corrected matrix in use for all downstream steps.\n")
} else {
  cat("Batch correction skipped (run_batch_correction = FALSE). Recall: plate",
      "does not predict survival directly (p = 0.43), only global",
      "methylation (p = 0.012) -- this is optional cleanup, not required to",
      "trust the survival results.\n")
}

gene_mat <- load_promoter_gene_mat(m_mat)

cat("\nTotal promoters (unfiltered):", nrow(gene_mat), "\n")

probe_var <- row_vars(gene_mat)
n_prescreened <- min(500, nrow(gene_mat))
top_var_idx   <- order(probe_var, decreasing = TRUE)[seq_len(n_prescreened)]
m_lasso       <- gene_mat[top_var_idx, ]

clin_mat      <- model.matrix(~ age_years + stage, data = meth_clinical)[, -1, drop = FALSE]
probe_mat_t   <- t(m_lasso)
X             <- cbind(clin_mat, probe_mat_t)
clin_colnames <- colnames(clin_mat)
gene_col_idx  <- which(!colnames(X) %in% clin_colnames)

penalty_factors <- c(rep(0, ncol(clin_mat)), rep(1, ncol(probe_mat_t)))

set.seed(42)
cv_fit <- cv.glmnet(
  x = X, y = surv_obj, family = "cox", penalty.factor = penalty_factors,
  alpha = 1, nfolds = 10, standardize = TRUE
)

selected_1se <- extract_probes(cv_fit, "lambda.1se", clin_colnames)
selected_min <- extract_probes(cv_fit, "lambda.min", clin_colnames)
if (length(selected_1se) > 0) {
  selected_probes <- selected_1se
  lambda_used <- "lambda.1se"
} else {
  selected_probes <- selected_min
  lambda_used <- "lambda.min"
}

cat("\nLegacy variance-filtered LASSO probes at", lambda_used, ":",
    length(selected_probes), "\n")
if (length(selected_probes) > 0) print(selected_probes)

if (length(selected_probes) > 1) {
  cor_mat  <- cor(t(gene_mat[selected_probes, , drop = FALSE]), method = "pearson")
  dist_mat <- as.dist(1 - abs(cor_mat))
  hc       <- hclust(dist_mat, method = "complete")
  clusters <- cutree(hc, h = 0.3)

  lasso_coefs_raw <- as.matrix(coef(cv_fit, s = lambda_used))
  best_probe <- vapply(unique(clusters), function(cluster) {
    probes <- names(clusters)[clusters == cluster]
    probes <- intersect(probes, rownames(lasso_coefs_raw))
    if (length(probes) == 0) return(NA_character_)

    abs_coefs <- abs(lasso_coefs_raw[probes, 1, drop = TRUE])
    if (all(is.na(abs_coefs))) return(NA_character_)
    probes[which.max(abs_coefs)]
  }, character(1))
  selected_probes <- unique(best_probe[!is.na(best_probe)])
}

cat("\nRepresentative legacy LASSO probe set:\n")
print(selected_probes)

screen_train <- function(gm, sv, clin, idx, n_keep, purity_available, use_subtype) {
  sv_i      <- sv[idx]
  age_i     <- clin$age_years[idx]
  stage_i   <- clin$stage[idx]
  purity_i  <- if (purity_available) clin$purity[idx] else NULL
  subtype_i <- if (use_subtype) droplevels(factor(clin$subtype[idx])) else NULL
  promoter_p <- apply(gm[, idx, drop = FALSE], 1, function(promoter) {
    fit <- uni_screen(promoter, sv_i, age_i, stage_i, purity = purity_i, subtype = subtype_i)
    fit["p"]
  })
  names(sort(promoter_p))[seq_len(n_keep)]
}

run_supervised_cv <- function(gene_mat, surv_obj, meth_clinical, K, n_repeats,
                              k_inner, n_screen, purity_available, use_subtype,
                              clin_terms, screen_train, n_cores, worker_log_dir) {
  cl <- make_worker_cluster(n_cores, worker_log_dir)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)

  foreach(
    r = seq_len(n_repeats),
    .packages = c("survival", "glmnet"),
    .export   = c("uni_screen")
  ) %dopar% {

    set.seed(r)
    folds_ids  <- integer(nrow(meth_clinical))
    events_idx <- which(meth_clinical$event == 1)
    cen_idx    <- which(meth_clinical$event == 0)
    folds_ids[events_idx] <- sample(rep(1:K, length.out = length(events_idx)))
    folds_ids[cen_idx]    <- sample(rep(1:K, length.out = length(cen_idx)))

    cv_risk_rep   <- numeric(nrow(meth_clinical))
    fold_selected <- character(0)
    fold_c        <- NULL
    fold_events   <- NULL

    for (k in seq_len(K)) {
      tr_idx <- which(folds_ids != k)
      te_idx <- which(folds_ids == k)
      if (use_subtype) {
        ev <- table(meth_clinical$subtype[tr_idx][meth_clinical$event[tr_idx] == 1])
        fold_events <- rbind(fold_events,
                             data.frame(repeat_id = r, fold = k, subtype = names(ev),
                                        train_events = as.integer(ev)))
      }

      keep_genes <- screen_train(gene_mat, surv_obj, meth_clinical, tr_idx,
                                  n_screen, purity_available, use_subtype)
      fold_selected <- c(fold_selected, keep_genes)

      clin_formula <- reformulate(clin_terms)
      clin_k <- model.matrix(clin_formula, data = meth_clinical)[, -1, drop = FALSE]
      stopifnot(nrow(clin_k) == nrow(meth_clinical))

      X_k  <- cbind(clin_k, t(gene_mat[keep_genes, , drop = FALSE]))
      pf_k <- c(rep(0, ncol(clin_k)), rep(1, length(keep_genes)))
      gene_cols_k <- (ncol(clin_k) + 1):ncol(X_k)

      fit_k <- cv.glmnet(X_k[tr_idx, ], surv_obj[tr_idx], family = "cox",
                          penalty.factor = pf_k, nfolds = k_inner,
                          standardize = TRUE)

      gene_names_k <- colnames(X_k)[gene_cols_k]
      coef_mat_1se <- as.matrix(coef(fit_k, s = "lambda.1se"))
      stopifnot("LASSO coefficient names do not cover the screened genes" =
                  all(gene_names_k %in% rownames(coef_mat_1se)))

      coef_used <- coef_mat_1se
      cf <- coef_used[gene_names_k, 1]
      if (all(cf == 0, na.rm = TRUE)) {
        coef_used <- as.matrix(coef(fit_k, s = "lambda.min"))
        stopifnot(all(gene_names_k %in% rownames(coef_used)))
        cf <- coef_used[gene_names_k, 1]
      }
      cv_risk_rep[te_idx] <- as.numeric(X_k[te_idx, gene_cols_k, drop = FALSE] %*% cf)

      stopifnot(all(colnames(X_k) %in% rownames(coef_used)))
      lp_mc_te  <- as.numeric(X_k[te_idx, , drop = FALSE] %*% coef_used[colnames(X_k), 1])
      clin_fit  <- coxph(reformulate(clin_terms, response = "Surv(time, event)"),
                         data = meth_clinical[tr_idx, ])
      lp_cl_te  <- predict(clin_fit, newdata = meth_clinical[te_idx, ], type = "lp")
      fold_c <- rbind(fold_c, c(
        repeat_id = r, fold = k, n_test = length(te_idx),
        events_test = sum(meth_clinical$event[te_idx]),
        c_clinical  = unname(concordance(surv_obj[te_idx] ~ lp_cl_te, reverse = TRUE)$concordance),
        c_meth_clin = unname(concordance(surv_obj[te_idx] ~ lp_mc_te, reverse = TRUE)$concordance)))
    }

    list(risk = cv_risk_rep, selected = fold_selected, fold_c = fold_c,
         fold_events = fold_events,
         c_index = concordance(surv_obj ~ cv_risk_rep, reverse = TRUE)$concordance)
  }
}

n_cores <- n_workers()
cat("\nRunning", n_repeats, "x", K, "-fold supervised screen across",
    n_cores, "cores...\n")

results <- run_supervised_cv(gene_mat, surv_obj, meth_clinical, K, n_repeats,
                             k_inner, n_screen, purity_available, use_subtype,
                             clin_terms, screen_train, n_cores, paths$run)

all_risk_scores <- vapply(results, function(x) x$risk,
                          numeric(nrow(meth_clinical)))
if (is.null(dim(all_risk_scores))) {
  all_risk_scores <- matrix(all_risk_scores, ncol = 1)
}
c_index_reps    <- vapply(results, function(x) x$c_index, numeric(1))
selected_freq   <- table(unlist(lapply(results, function(x) x$selected)))
mean_cv_risk    <- rowMeans(all_risk_scores)

if (use_subtype) {
  fold_events <- do.call(rbind, lapply(results, `[[`, "fold_events"))
  cat("\n-- Events per subtype per outer training fold --\n")
  print(xtabs(train_events ~ subtype + interaction(repeat_id, fold, sep = "/", lex.order = TRUE),
              data = fold_events))
  low <- fold_events[fold_events$train_events < 5, ]
  if (nrow(low) > 0) {
    warning(sprintf("%d training fold(s) have a subtype with fewer than 5 events (%s); its subtype coefficient is poorly estimated in those folds.",
                    nrow(low), paste(sprintf("repeat %d fold %d: %s %d", low$repeat_id, low$fold,
                                             low$subtype, low$train_events), collapse = "; ")),
            call. = FALSE)
  }
}

cat("\n-- Per-repeat results (supervised screen) --\n")
for (r in seq_len(n_repeats)) {
  cat(sprintf("  Repeat %2d  C-index: %.3f\n", r, c_index_reps[r]))
}

cat("\nMost consistently screened-in promoters (of", n_repeats * K, "folds):\n")
print(head(sort(selected_freq, decreasing = TRUE), 25))

cat("\nOverlap between legacy variance-LASSO set and top supervised-screen genes:\n")
top_supervised <- names(sort(selected_freq, decreasing = TRUE))[seq_len(min(25, length(selected_freq)))]
print(intersect(selected_probes, top_supervised))

heldout <- as.data.frame(do.call(rbind, lapply(results, `[[`, "fold_c")))
heldout$c_diff <- heldout$c_meth_clin - heldout$c_clinical
c_ho_clin <- mean(heldout$c_clinical)
c_ho_mc   <- mean(heldout$c_meth_clin)
c_ho_diff <- mean(heldout$c_diff)
cat("\n-- Held-out C-index (outer CV folds,", nrow(heldout), "folds) --\n")
cat(sprintf("  Clinical baseline (%s):         %.3f\n", paste(clin_terms, collapse = " + "), c_ho_clin))
cat(sprintf("  Methylation + clinical baseline (LASSO):  %.3f\n", c_ho_mc))
cat(sprintf("  Difference, mean (SD over folds):         %+.3f (%.3f); better in %d of %d folds\n",
            c_ho_diff, if (nrow(heldout) > 1) sd(heldout$c_diff) else NA_real_,
            sum(heldout$c_diff > 0), nrow(heldout)))

clinical_formula <- reformulate(clin_terms, response = "Surv(time, event)")
clinical_model <- coxph(clinical_formula, data = meth_clinical, x = TRUE, y = TRUE)
c_clinical     <- summary(clinical_model)$concordance["C"]
c_mean_score   <- concordance(surv_obj ~ mean_cv_risk, reverse = TRUE)$concordance

cat("\n-- C-index summary (in-sample unless stated) --\n")
cat(sprintf("  Clinical only, in-sample (%s): %.3f\n", paste(clin_terms, collapse = " + "), c_clinical))
cat(sprintf("  CV methylation score (mean repeats):  %.3f\n", c_mean_score))
cat(sprintf("  Per-repeat mean (SD):                 %.3f (%.3f)\n",
            mean(c_index_reps), sd(c_index_reps)))
cat(sprintf("  95%% empirical CI:                     %.3f - %.3f\n",
            quantile(c_index_reps, 0.025), quantile(c_index_reps, 0.975)))

model_df <- cbind(
  meth_clinical[, c("time", "event", "age_years", "stage")],
  cv_meth_risk_score = mean_cv_risk
)
if (purity_available) model_df$purity <- meth_clinical$purity
if (use_subtype) model_df$subtype <- meth_clinical$subtype

final_formula <- reformulate(c("cv_meth_risk_score", clin_terms),
                             response = "Surv(time, event)")
final_model <- coxph(final_formula, data = model_df, x = TRUE, y = TRUE)

cat("\n-- Final Cox model --\n")
print(summary(final_model))

c_final <- summary(final_model)$concordance["C"]
cat(sprintf("\nFinal model C-index, in-sample (Cox refit on the full cohort): %.3f (SE: %.3f)\n",
            c_final, summary(final_model)$concordance["se(C)"]))

cat(sprintf(paste0("\n-- Headline (subtype in baseline: %s) --\n",
                   "  Held-out (outer CV folds):\n",
                   "    Clinical baseline C-index:               %.3f\n",
                   "    Methylation + clinical baseline:         %.3f\n",
                   "    Delta (methylation over baseline):       %+.3f\n",
                   "  In-sample (full cohort, optimistic):\n",
                   "    Clinical baseline C-index:               %.3f\n",
                   "    Methylation score + clinical baseline:   %.3f\n",
                   "    Delta (methylation over baseline):       %+.3f\n"),
            use_subtype, c_ho_clin, c_ho_mc, c_ho_diff,
            c_clinical, c_final, c_final - c_clinical))

eval_times_requested <- c(365, 730, 1825, 3652)
max_event_t <- max(meth_clinical$time[meth_clinical$event == 1], na.rm = TRUE)
eval_times  <- eval_times_requested[eval_times_requested < max_event_t]

dropped_times <- setdiff(eval_times_requested, eval_times)
if (length(dropped_times) > 0) {
  warning(sprintf(
    "Dropping evaluation time(s) beyond the last event (%d days): %s",
    round(max_event_t), paste(dropped_times, collapse = ", ")), call. = FALSE)
}
stopifnot("No evaluation time falls within observed follow-up" =
            length(eval_times) > 0)

eval_labels <- c("365" = "1yr", "730" = "2yr", "1825" = "5yr", "3652" = "10yr")
eval_cols   <- c("365" = "blue", "730" = "red",
                 "1825" = "darkgreen", "3652" = "purple")

timeroc_result <- timeROC(T = meth_clinical$time, delta = meth_clinical$event,
                           marker = mean_cv_risk, cause = 1,
                           weighting = "marginal",
                           times = eval_times,
                           ROC = TRUE, iid = FALSE)

brier_score <- Score(
  object = list("Methylation" = final_model, "Clinical" = clinical_model),
  formula = Surv(time, event) ~ 1, data = model_df,
  metrics = "brier", summary = "IPA", cause = 1,
  times = eval_times,
  plots = "calibration"
)

cat(sprintf("AUC at %s: %.3f\n", eval_labels[as.character(eval_times[1])],
            timeroc_result$AUC[1]))
print(brier_score$Brier$score)

risk_group <- factor(ifelse(mean_cv_risk >= median(mean_cv_risk), "High", "Low"),
                      levels = c("Low", "High"))
km_fit  <- survfit(surv_obj ~ risk_group)
lr_test <- survdiff(surv_obj ~ risk_group)

lr_pval <- pchisq(lr_test$chisq, df = length(lr_test$n) - 1, lower.tail = FALSE)

with_pdf(file.path(paths$diagnostics, "methylation_risk_groups.pdf"), {
  plot(km_fit, col = c("steelblue", "tomato"), lty = 1, lwd = 2,
       xlab = "Time (days)", ylab = "Progression-free probability",
       main = "Kaplan-Meier: Methylation Risk Score Groups (Median Split)",
       mark.time = TRUE)
  legend("topright",
         legend = c(paste0("Low risk  (n=", sum(risk_group == "Low"), ")"),
                    paste0("High risk (n=", sum(risk_group == "High"), ")")),
         col = c("steelblue", "tomato"), lty = 1, lwd = 2, bty = "n")
  mtext(paste0("Log-rank p = ", format.pval(lr_pval, digits = 3)), side = 3, line = 0.3)
}, width = 7, height = 6)

with_pdf(file.path(paths$diagnostics, "c_index_stability_methylation.pdf"), {
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
}, width = 6, height = 5)

with_pdf(file.path(paths$diagnostics, "timeroc_curves.pdf"), {
  for (i in seq_along(eval_times)) {
    tt <- eval_times[i]
    plot(timeroc_result, time = tt, col = eval_cols[as.character(tt)],
         add = (i > 1), title = (i == 1), lwd = 2)
  }
  legend("bottomright", legend = eval_labels[as.character(eval_times)],
         col = eval_cols[as.character(eval_times)], lwd = 2, bty = "n")
}, width = 6, height = 6)

with_pdf(file.path(paths$diagnostics, "calibration_plot.pdf"), {
  plotCalibration(x = brier_score, models = c("Methylation", "Clinical"),
                   times = max(eval_times), method = "quantile", q = 5)
}, width = 6, height = 6)

with_pdf(file.path(paths$diagnostics, "auc_over_time.pdf"), {
  plot(x = eval_times, y = timeroc_result$AUC, type = "b", pch = 19,
       col = "steelblue", lwd = 2, xlab = "Time (days)", ylab = "AUC(t)",
       ylim = c(0.4, 1), main = "Time-dependent AUC: Methylation Risk Score")
  abline(h = 0.5, lty = 2, col = "gray50")
}, width = 6, height = 5)

with_pdf(file.path(paths$diagnostics, "forest_plot_final_model.pdf"), {
  print(ggforest(final_model, data = model_df))
}, width = 7, height = 5)

cat("\nPlots saved to", paths$diagnostics, "\n")

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

write.csv(heldout, file.path(paths$tables, "heldout_c_index_per_fold.csv"),
          row.names = FALSE)

params <- c(K = K, n_repeats = n_repeats, K_inner = k_inner, n_screen = n_screen,
            endpoint = "PFI", seed = 42,
            purity_used = purity_available, subtype_used = use_subtype,
            heldout_c_clinical  = round(c_ho_clin, 4),
            heldout_c_meth_clin = round(c_ho_mc, 4),
            heldout_c_diff      = round(c_ho_diff, 4),
            insample_c_clinical = round(unname(c_clinical), 4),
            insample_c_final    = round(unname(c_final), 4),
            batch_corrected = run_batch_correction,
            eval_times = paste(eval_times, collapse = ";"))
write_manifest(paths$run, params, inputs = METH_SCREEN_INPUTS)

cat("\nOutputs saved to", paths$run, "\n")
