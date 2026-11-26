library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(doParallel)
library(foreach)

source("deployed_scripts_v1/00_utils.R")

paths <- make_run_dir(tag = "stratified_screen_validation", base = "02_analysis")
set.seed(42)
code_md5 <- script11_code_hash()

N_PERM <- full_or_test(100, 5)

m_mat         <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

gene_mat <- load_promoter_gene_mat(m_mat)

meth_clinical <- load_subtype_purity(meth_clinical)

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)

meth_c     <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c     <- surv_obj[complete]

cat("Analysis cohort:", sum(complete), "patients,",
    sum(meth_c$event), "events\n")

cat("\nRunning observed stratified screen on", nrow(gene_mat_c), "promoters...\n")
t0 <- Sys.time()
obs_res <- run_full_screen(gene_mat_c, surv_c, meth_c)
cat("Done in", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

obs_padj    <- p.adjust(obs_res[, "p"], method = "BH")
obs_n_hits  <- sum(obs_padj < 0.10, na.rm = TRUE)
obs_min_p   <- min(obs_res[, "p"], na.rm = TRUE)

cat("Observed: ", obs_n_hits, "promoters at BH q<0.10\n")

cat("\n=== TEST A: Freedman-Lane permutation null (", N_PERM, "permutations) ===\n")

n_cores <- n_workers()
cat("Running on", n_cores, "cores...\n")

fl <- freedman_lane_prepare(gene_mat_c, meth_c)

t0 <- Sys.time()
perm_out <- run_fl_permutations(fl, surv_c, meth_c, N_PERM, n_cores,
                                run_full_screen, stratified_screen_one,
                                freedman_lane_perm, paths$run)
cat("Permutation run:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

perm_results <- perm_out$results
n_failed     <- perm_out$n_failed
n_perm_ok    <- nrow(perm_results)

perm_p_hitcount <- (sum(perm_results$n_hits >= obs_n_hits) + 1) / (n_perm_ok + 1)
perm_p_minp     <- (sum(perm_results$min_p <= obs_min_p) + 1) / (n_perm_ok + 1)

cat("\nNull distribution of hit counts (BH q<0.10) across", n_perm_ok, "permutations:\n")
print(summary(perm_results$n_hits))
cat(sprintf("\nObserved: %d hits. Permutation p-value: %.3f\n",
            obs_n_hits, perm_p_hitcount))
cat(sprintf("Observed smallest p-value: %.3g. Permutation p-value: %.3f\n",
            obs_min_p, perm_p_minp))
cat("INTERPRETATION: a small perm_p_hitcount (e.g. <0.05) means the real",
    "screen finds more hits than chance alone typically produces under this",
    "exact procedure -- the candidate list clears the bar. If perm_p_hitcount",
    "is large, the observed hit count is unremarkable noise and the candidates should not be",
    "trusted without further validation.\n")

pdf(file.path(paths$diagnostics, "permutation_null_hit_counts.pdf"),
    width = 7, height = 5)
hist(perm_results$n_hits, breaks = max(10, n_perm_ok / 5), col = "grey80",
     border = "white", xlab = "Promoters at BH q<0.10 (null permutations)",
     main = sprintf("Permutation null (n=%d) vs observed (perm p=%.3f)",
                     n_perm_ok, perm_p_hitcount),
     xlim = range(c(perm_results$n_hits, obs_n_hits)))
abline(v = obs_n_hits, col = "tomato", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Observed = %d", obs_n_hits),
       col = "tomato", lty = 2, lwd = 2, bty = "n")
dev.off()

cat("\n=== TEST B: stability with Her2 excluded ===\n")

candidates <- screen_candidates(obs_res[, "p"])
cat("Candidates (this run's stratified screen, BH q<0.10):",
    if (length(candidates) > 0) paste(candidates, collapse = ", ") else "none", "\n")

no_her2 <- meth_c$subtype != "Her2"
meth_no_her2     <- meth_c[no_her2, ]
meth_no_her2$subtype <- droplevels(meth_no_her2$subtype)
gene_mat_no_her2 <- gene_mat_c[, no_her2, drop = FALSE]
surv_no_her2     <- surv_c[no_her2]

cat("Cohort with Her2 excluded:", sum(no_her2), "patients,",
    sum(meth_no_her2$event), "events\n")

if (length(candidates) == 0) {
  cat("No candidates at BH q<0.10 -- TEST B skipped.\n")
  stability <- NULL
} else {
  stability <- data.frame(
    promoter          = candidates,
    coef_with_her2    = obs_res[candidates, "coef"],
    p_with_her2       = obs_res[candidates, "p"],
    padj_with_her2    = obs_padj[candidates],
    coef_without_her2 = NA_real_,
    p_without_her2    = NA_real_,
    stringsAsFactors  = FALSE
  )

  for (g in candidates) {
    r <- stratified_screen_one(gene_mat_no_her2[g, ], surv_no_her2, meth_no_her2)
    stability[stability$promoter == g, "coef_without_her2"] <- r["coef"]
    stability[stability$promoter == g, "p_without_her2"]    <- r["p"]
  }

  stability$same_direction <- sign(stability$coef_with_her2) ==
    sign(stability$coef_without_her2)
  stability$still_p05 <- stability$p_without_her2 < 0.05

  cat("\n")
  print(stability)
  cat("\nINTERPRETATION: same_direction = FALSE or still_p05 = FALSE flags a",
      "candidate whose result is not stable to removing the Her2 stratum --",
      "treat it as substantially weaker evidence than its q-value alone",
      "suggests.\n")
}

write.csv(perm_results,
          file.path(paths$tables, "permutation_null_results.csv"),
          row.names = FALSE)
if (!is.null(stability)) {
  write.csv(stability, file.path(paths$tables, "her2_exclusion_stability.csv"),
            row.names = FALSE)
}

write_manifest(paths$run, list(
  script                = "11_stratified_screen_validation.R",
  n_permutations        = N_PERM,
  n_failed              = n_failed,
  n_failed_fits         = perm_out$n_failed_fits,
  median_fit_fail_frac  = perm_out$median_fit_fail_frac,
  permutation_null      = "freedman_lane_within_strata",
  n_cohort              = sum(complete),
  n_promoters           = nrow(gene_mat_c),
  screen_input_md5      = screen_input_hash(gene_mat_c, surv_c, meth_c),
  script_md5            = code_md5,
  observed_hits_fdr10   = obs_n_hits,
  permutation_p_hitcount = round(perm_p_hitcount, 4),
  permutation_p_minp    = round(perm_p_minp, 4),
  candidate_source      = "derived: this run's stratified screen, BH q<0.10",
  candidates            = if (length(candidates) > 0) paste(candidates, collapse = ";") else "none",
  n_candidates_tested   = length(candidates),
  n_stable_both_checks  = if (is.null(stability)) 0 else
                            sum(stability$same_direction & stability$still_p05),
  seed                  = 42
), inputs = METH_SCREEN_INPUTS)

cat("\n=== Done. Results in", paths$run, "===\n")
