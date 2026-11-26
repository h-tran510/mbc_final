library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(doParallel)
library(foreach)

source("deployed_scripts_v1/00_utils.R")

paths <- make_run_dir(tag = "sign_bias_diagnostics", base = "02_analysis")
fig_dir <- "03_figures/sign_bias"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

N_PERM  <- full_or_test(200, 5)
n_cores <- n_workers()
N_PROM_TEST  <- 400
STORE_PMAT   <- N_PERM <= 200

m_mat         <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

gene_mat <- load_promoter_gene_mat(m_mat)
rm(m_mat)

meth_clinical <- load_subtype_purity(meth_clinical)

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)

meth_c     <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c     <- surv_obj[complete]
meth_c$subtype <- droplevels(meth_c$subtype)
meth_c$stage   <- droplevels(meth_c$stage)

cat("Analysis cohort:", sum(complete), "patients,",
    sum(meth_c$event), "events,", nrow(gene_mat_c), "promoters\n")

stratified_pc1_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity) + scale(meth_x$pc1) +
                      strata(meth_x$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

frac_neg <- function(coef) mean(coef[!is.na(coef)] < 0)

pca <- prcomp(t(gene_mat), center = TRUE, scale. = TRUE, rank. = 3)
if (cor(pca$x[, 1], colMeans(gene_mat)) < 0) {
  pca$x[, 1]        <- -pca$x[, 1]
  pca$rotation[, 1] <- -pca$rotation[, 1]
  cat("PC1 sign flipped so that higher PC1 = higher global methylation.\n")
}
var_expl <- pca$sdev[1:3]^2 / sum(pca$sdev^2)
stopifnot("PCA score rows are not in meth_clinical sample order" =
            identical(rownames(pca$x), meth_clinical$meth_sample_barcode))
meth_c$pc1 <- pca$x[complete, 1]
cat(sprintf("\nVariance explained: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%\n",
            100 * var_expl[1], 100 * var_expl[2], 100 * var_expl[3]))

screen_rows <- if (TEST_MODE) {
  sample(rownames(gene_mat_c), min(N_PROM_TEST, nrow(gene_mat_c)))
} else rownames(gene_mat_c)
gm_s <- gene_mat_c[screen_rows, , drop = FALSE]
cat("Promoters screened:", nrow(gm_s),
    if (TEST_MODE) "(test-mode subset)" else "", "\n")

cat("\n=== (a) Naive vs stratified vs stratified + PC1 ===\n")
t0 <- Sys.time()
naive_res <- run_full_screen(gm_s, surv_c, meth_c, naive_screen_one)
strat_res <- run_full_screen(gm_s, surv_c, meth_c, stratified_screen_one)
pc1_res   <- run_full_screen(gm_s, surv_c, meth_c, stratified_pc1_screen_one)
cat("Three observed screens:", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

obs_frac <- c(naive = frac_neg(naive_res[, "coef"]),
              stratified = frac_neg(strat_res[, "coef"]),
              stratified_pc1 = frac_neg(pc1_res[, "coef"]))
cat(sprintf("Negative-coefficient fraction: naive %.1f%% | stratified %.1f%% | stratified+PC1 %.1f%%\n",
            100 * obs_frac[1], 100 * obs_frac[2], 100 * obs_frac[3]))

bh_hits <- function(p) sum(p.adjust(p, method = "BH") < 0.10, na.rm = TRUE)
cat(sprintf("BH q<0.10 hits: naive %d | stratified %d | stratified+PC1 %d\n",
            bh_hits(naive_res[, "p"]), bh_hits(strat_res[, "p"]), bh_hits(pc1_res[, "p"])))

candidates <- screen_candidates(strat_res[, "p"])
cat("Candidates (stratified screen, BH q<0.10):",
    if (length(candidates) > 0) paste(candidates, collapse = ", ") else "none", "\n")

screen_table <- data.frame(
  promoter   = screen_rows,
  pc1_loading = pca$rotation[screen_rows, 1],
  naive_coef = naive_res[, "coef"], naive_p = naive_res[, "p"],
  strat_coef = strat_res[, "coef"], strat_p = strat_res[, "p"],
  strat_pc1_coef = pc1_res[, "coef"], strat_pc1_p = pc1_res[, "p"],
  stringsAsFactors = FALSE
)

ok <- !is.na(screen_table$strat_coef)
cor_coef_load <- cor(screen_table$strat_coef[ok], screen_table$pc1_loading[ok])
cat(sprintf("cor(stratified coef, PC1 loading) = %.3f\n", cor_coef_load))

tercile <- cut(screen_table$pc1_loading,
               quantile(screen_table$pc1_loading, c(0, 1/3, 2/3, 1)),
               include.lowest = TRUE, labels = c("low", "mid", "high"))
by_tercile <- sapply(split(screen_table, tercile), function(d)
  c(n = nrow(d), strat_neg = frac_neg(d$strat_coef), strat_pc1_neg = frac_neg(d$strat_pc1_coef)))
cat("Negative fraction by PC1-loading tercile:\n")
print(round(by_tercile, 3))

if (length(candidates) > 0) {
  cat("\nCandidates in each screen:\n")
  print(screen_table[screen_table$promoter %in% candidates, ], row.names = FALSE)
} else {
  cat("\nNo candidates at BH q<0.10 -- candidate table and (d) skipped.\n")
}

cat("\n=== (b)+(c) Freedman-Lane vs column-permutation null (", N_PERM,
    "permutations) ===\n")

fl <- freedman_lane_prepare(gm_s, meth_c)

run_two_nulls <- function(gm_s, fl, surv_c, meth_c, N_PERM,
                          n_cores, store_pmat, run_full_screen,
                          stratified_screen_one, freedman_lane_perm,
                          worker_log_dir, fit_fail_max = PERM_FIT_FAIL_MAX) {
  force(fit_fail_max)
  cl <- make_worker_cluster(n_cores, worker_log_dir)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)

  foreach(
    r = seq_len(N_PERM),
    .packages = c("survival"),
    .export   = c("run_full_screen", "stratified_screen_one", "freedman_lane_perm",
                  "strata_perm")
  ) %dopar% {
    tryCatch({
      set.seed(1000 + r)
      perm_idx <- sample(ncol(gm_s))

      summarise <- function(res) {
        c(n_hits = sum(p.adjust(res[, "p"], method = "BH") < 0.10, na.rm = TRUE),
          min_p  = min(res[, "p"], na.rm = TRUE),
          frac_neg = mean(res[!is.na(res[, "coef"]), "coef"] < 0))
      }

      res_col <- run_full_screen(gm_s[, perm_idx, drop = FALSE], surv_c, meth_c,
                                 stratified_screen_one)

      gm_fl  <- freedman_lane_perm(fl, perm_idx)
      res_fl <- run_full_screen(gm_fl, surv_c, meth_c, stratified_screen_one)

      n_fit_failed <- max(sum(is.na(res_col[, "p"])), sum(is.na(res_fl[, "p"])))
      n_fits       <- nrow(res_col)
      if (n_fit_failed > fit_fail_max * n_fits) {
        message(sprintf("perm %d seed %d failed: %d of %d Cox fits returned NA",
                        r, 1000 + r, n_fit_failed, n_fits))
        list(permutation = r, failed = TRUE, n_fit_failed = n_fit_failed, n_fits = n_fits)
      } else {
        list(permutation = r, failed = FALSE,
             n_fit_failed = n_fit_failed, n_fits = n_fits,
             col = summarise(res_col), fl = summarise(res_fl),
             p_fl = if (store_pmat) res_fl[, "p"] else NULL)
      }
    }, error = function(e) {
      message(sprintf("perm %d seed %d failed: %s", r, 1000 + r, conditionMessage(e)))
      list(permutation = r, failed = TRUE, n_fit_failed = NA, n_fits = NA)
    })
  }
}

t0 <- Sys.time()
perm_out <- run_two_nulls(gm_s, fl, surv_c, meth_c, N_PERM,
                          n_cores, STORE_PMAT, run_full_screen,
                          stratified_screen_one, freedman_lane_perm, paths$run)
cat("Permutation run:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

perm_failed <- vapply(perm_out, function(x) isTRUE(x$failed), logical(1))
failed_fits <- vapply(perm_out, function(x) as.numeric(x$n_fit_failed), numeric(1))
total_fits <- vapply(perm_out, function(x) as.numeric(x$n_fits), numeric(1))
perm_fits <- perm_fit_failures(failed_fits, total_fits, perm_failed)

failed_ids <- vapply(perm_out[perm_failed], function(x) x$permutation, numeric(1))
n_failed <- check_perm_failures(failed_ids, N_PERM)
perm_out <- perm_out[!perm_failed]
perm_ids <- vapply(perm_out, function(x) x$permutation, numeric(1))
n_perm_ok <- length(perm_out)

null_col <- as.data.frame(do.call(rbind, lapply(perm_out, `[[`, "col")))
null_fl <- as.data.frame(do.call(rbind, lapply(perm_out, `[[`, "fl")))
p_mat_fl <- NULL
if (STORE_PMAT) {
  p_mat_fl <- do.call(cbind, lapply(perm_out, `[[`, "p_fl"))
}

obs_hits <- bh_hits(strat_res[, "p"])
obs_minp <- min(strat_res[, "p"], na.rm = TRUE)
obs_fn   <- obs_frac[["stratified"]]

perm_p_ge <- function(null, obs) (sum(null >= obs) + 1) / (length(null) + 1)
perm_p_le <- function(null, obs) (sum(null <= obs) + 1) / (length(null) + 1)
perm_p_2s <- function(null, obs) min(1, 2 * min(perm_p_ge(null, obs), perm_p_le(null, obs)))

null_compare <- data.frame(
  null = c("column_permutation", "freedman_lane_within_strata"),
  obs_hits_q10 = obs_hits,
  null_hits_median = c(median(null_col$n_hits), median(null_fl$n_hits)),
  null_hits_mean   = c(mean(null_col$n_hits), mean(null_fl$n_hits)),
  perm_p_hits  = c(perm_p_ge(null_col$n_hits, obs_hits), perm_p_ge(null_fl$n_hits, obs_hits)),
  obs_min_p    = obs_minp,
  perm_p_min_p = c(perm_p_le(null_col$min_p, obs_minp), perm_p_le(null_fl$min_p, obs_minp)),
  obs_frac_neg = obs_fn,
  null_frac_neg_median = c(median(null_col$frac_neg), median(null_fl$frac_neg)),
  null_frac_neg_q025 = c(quantile(null_col$frac_neg, 0.025), quantile(null_fl$frac_neg, 0.025)),
  null_frac_neg_q975 = c(quantile(null_col$frac_neg, 0.975), quantile(null_fl$frac_neg, 0.975)),
  obs_frac_neg_quantile = c(mean(null_col$frac_neg <= obs_fn), mean(null_fl$frac_neg <= obs_fn)),
  perm_p_frac_neg_2s = c(perm_p_2s(null_col$frac_neg, obs_fn), perm_p_2s(null_fl$frac_neg, obs_fn)),
  stringsAsFactors = FALSE
)

cat("\n(b) Sign fraction: observed stratified negative fraction vs each null\n")
cat("    (replaces script 10's binom.test, which assumed independent promoters)\n")
for (i in seq_len(nrow(null_compare))) {
  with(null_compare[i, ], cat(sprintf(
    "  %-27s observed %.3f | null median %.3f, 95%% [%.3f, %.3f] | observed at null quantile %.2f | 2-sided perm p %.3f\n",
    null, obs_frac_neg, null_frac_neg_median, null_frac_neg_q025, null_frac_neg_q975,
    obs_frac_neg_quantile, perm_p_frac_neg_2s)))
}

cat("\n(c) Hit count and min-p: Freedman-Lane vs column permutation\n")
for (i in seq_len(nrow(null_compare))) {
  with(null_compare[i, ], cat(sprintf(
    "  %-27s observed %d hits | null median %.1f, mean %.2f | perm p (hits) %.3f | perm p (min-p) %.3f\n",
    null, obs_hits_q10, null_hits_median, null_hits_mean, perm_p_hits, perm_p_min_p)))
}
cat("INTERPRETATION: the Freedman-Lane null keeps methylation's association",
    "with purity/subtype/age/stage, so it is the right null for 'associated",
    "with PFI given the covariates'. If its hit-count / min-p p-values are",
    "larger than the column-permutation ones, script 11's old null was anti-conservative.",
    "The two nulls also differ in exchangeability (within subtype vs across all",
    "patients), so such a difference reflects both changes together, not",
    "Freedman-Lane alone (S11).\n")

wy_table <- NULL
if (STORE_PMAT && length(candidates) > 0) {
  p_obs <- strat_res[, "p"]
  p_obs[is.na(p_obs)] <- 1
  p_mat_fl[is.na(p_mat_fl)] <- 1
  ord <- order(p_obs)
  q_mat <- apply(p_mat_fl[ord, , drop = FALSE], 2, function(x) rev(cummin(rev(x))))
  q_mat <- matrix(q_mat, nrow = length(ord))
  adj <- (rowSums(q_mat <= p_obs[ord]) + 1) / (n_perm_ok + 1)
  adj <- cummax(adj)
  wy_adj <- setNames(numeric(length(p_obs)), rownames(strat_res))
  wy_adj[ord] <- adj
  wy_table <- data.frame(promoter = candidates, strat_p = p_obs[candidates],
                         bh_q = p.adjust(p_obs, method = "BH")[candidates],
                         wy_minp_adj = wy_adj[candidates], stringsAsFactors = FALSE)
  cat("\n=== (d) Westfall-Young minP step-down adjusted p (Freedman-Lane null) ===\n")
  cat("Resolution is 1/(permutations+1) =", signif(1 / (n_perm_ok + 1), 3), "\n")
  print(wy_table, row.names = FALSE)
} else if (!STORE_PMAT) {
  cat("\n(d) skipped: N_PERM > 200, p-matrix not stored.\n")
}

with_pdf(file.path(fig_dir, "sign_fraction_null.pdf"), width = 10, height = 5, {
  par(mfrow = c(1, 2))
  xl <- range(c(null_col$frac_neg, null_fl$frac_neg, obs_fn))
  hist(null_col$frac_neg, breaks = 20, col = "grey80", border = "white", xlim = xl,
       xlab = "Negative-coefficient fraction", main = "Column-permutation null")
  abline(v = obs_fn, col = "tomato", lwd = 2, lty = 2)
  hist(null_fl$frac_neg, breaks = 20, col = "steelblue", border = "white", xlim = xl,
       xlab = "Negative-coefficient fraction", main = "Freedman-Lane null")
  abline(v = obs_fn, col = "tomato", lwd = 2, lty = 2)
})

write.csv(screen_table, file.path(paths$tables, "naive_strat_strat_pc1_screen.csv"),
          row.names = FALSE)
write.csv(null_compare, file.path(paths$tables, "null_comparison.csv"), row.names = FALSE)
write.csv(data.frame(permutation = perm_ids,
                     col = null_col, fl = null_fl),
          file.path(paths$tables, "permutation_null_results.csv"), row.names = FALSE)
if (!is.null(wy_table)) {
  write.csv(wy_table, file.path(paths$tables, "westfall_young_candidates.csv"),
            row.names = FALSE)
}

saveRDS(list(
  var_explained = var_expl, pc_scores = pca$x[complete, , drop = FALSE],
  pc_loadings = pca$rotation,
  screen_table = screen_table, obs_frac_neg = obs_frac,
  cor_coef_pc1_loading = cor_coef_load, frac_neg_by_tercile = by_tercile,
  null_col = null_col, null_fl = null_fl, null_compare = null_compare,
  p_mat_fl = p_mat_fl, westfall_young = wy_table
), file.path(paths$models, "sign_bias_diagnostics.rds"))

write_manifest(paths$run, list(
  script              = "14_sign_bias_diagnostics.R",
  n_permutations      = N_PERM,
  n_failed            = n_failed,
  n_failed_fits       = perm_fits$n_failed_fits,
  median_fit_fail_frac = perm_fits$median_fit_fail_frac,
  permutation_null    = "freedman_lane_within_strata vs column_permutation",
  n_promoters_screened = nrow(gm_s),
  candidate_source    = "derived: this run's stratified screen, BH q<0.10",
  candidates          = if (length(candidates) > 0) paste(candidates, collapse = ";") else "none",
  test_mode           = TEST_MODE,
  frac_neg_naive      = round(obs_frac[["naive"]], 4),
  frac_neg_stratified = round(obs_frac[["stratified"]], 4),
  frac_neg_strat_pc1  = round(obs_frac[["stratified_pc1"]], 4),
  seed                = 42
), inputs = METH_SCREEN_INPUTS,
   outputs = c(files_under(paths$run), files_under(fig_dir)))

cat("\n=== Done. Results in", paths$run, "and", fig_dir, "===\n")
