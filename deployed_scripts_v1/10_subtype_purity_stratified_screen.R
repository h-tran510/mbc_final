library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

source("deployed_scripts_v1/00_utils.R")

paths <- make_run_dir(tag = "subtype_purity_stratified_screen", base = "02_analysis")
set.seed(42)

m_mat         <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

gene_mat <- load_promoter_gene_mat(m_mat)

cat("Promoters:", nrow(gene_mat), " Samples:", ncol(gene_mat), "\n")

meth_clinical <- load_subtype_purity(meth_clinical)

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)

cat("\nComplete cases (subtype + purity + age + stage):",
    sum(complete), "of", nrow(meth_clinical), "\n")

meth_c   <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c   <- surv_obj[complete]

meth_c$subtype <- droplevels(meth_c$subtype)

cat("\nPatients and EVENTS per subtype (check before trusting per-promoter results):\n")
print(table(subtype = meth_c$subtype, event = meth_c$event))
event_counts <- tapply(meth_c$event, meth_c$subtype, sum)
stopifnot("Empty subtype stratum after complete-case filtering" =
            all(!is.na(event_counts)))
if (any(event_counts < 10)) {
  cat("\nWARNING: at least one subtype has fewer than 10 events.",
      "Coefficients for promoters whose signal is concentrated in that",
      "stratum will be unstable -- treat individual hits from that group",
      "with caution regardless of p-value.\n")
}

n_prom <- nrow(gene_mat_c)
cat("\nRunning naive screen on", n_prom, "promoters...\n")
t0 <- Sys.time()
naive_res <- run_full_screen(gene_mat_c, surv_c, meth_c, naive_screen_one)
cat("Naive screen:", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

cat("\nRunning subtype+purity-stratified screen on", n_prom, "promoters...\n")
t0 <- Sys.time()
strat_res <- run_full_screen(gene_mat_c, surv_c, meth_c)
cat("Stratified screen:", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

results <- data.frame(
  promoter        = rownames(gene_mat_c),
  naive_coef      = naive_res[, "coef"],
  naive_p         = naive_res[, "p"],
  stratified_coef = strat_res[, "coef"],
  stratified_p    = strat_res[, "p"],
  stringsAsFactors = FALSE
)

results$naive_padj      <- p.adjust(results$naive_p, method = "BH")
results$stratified_padj <- p.adjust(results$stratified_p, method = "BH")

cat("\n=== Sign bias ===\n")
naive_signs <- sign(results$naive_coef[!is.na(results$naive_coef)])
strat_signs <- sign(results$stratified_coef[!is.na(results$stratified_coef)])

bt_naive <- binom.test(sum(naive_signs < 0), length(naive_signs), p = 0.5)
bt_strat <- binom.test(sum(strat_signs < 0), length(strat_signs), p = 0.5)

cat(sprintf("Naive:      %d/%d negative (%.1f%%), binom p = %.3g\n",
            sum(naive_signs < 0), length(naive_signs),
            100 * mean(naive_signs < 0), bt_naive$p.value))
cat(sprintf("Stratified: %d/%d negative (%.1f%%), binom p = %.3g\n",
            sum(strat_signs < 0), length(strat_signs),
            100 * mean(strat_signs < 0), bt_strat$p.value))
cat("INTERPRETATION: if stratified % moves toward 50% and binom p rises",
    "toward non-significant, subtype+purity explain the sign bias.",
    "If it stays lopsided, something else is still driving it.\n")

cat("\n=== Hit count / pi0 ===\n")
pi0_naive_p <- results$naive_p[!is.na(results$naive_p)]
pi0_strat_p <- results$stratified_p[!is.na(results$stratified_p)]
lambda_cut  <- 0.5
pi0_naive   <- min(mean(pi0_naive_p > lambda_cut) / (1 - lambda_cut), 1)
pi0_strat   <- min(mean(pi0_strat_p > lambda_cut) / (1 - lambda_cut), 1)

cat(sprintf("Naive:      raw p<0.05: %d | BH q<0.10: %d | pi0 est: %.3f\n",
            sum(results$naive_p < 0.05, na.rm = TRUE),
            sum(results$naive_padj < 0.10, na.rm = TRUE), pi0_naive))
cat(sprintf("Stratified: raw p<0.05: %d | BH q<0.10: %d | pi0 est: %.3f\n",
            sum(results$stratified_p < 0.05, na.rm = TRUE),
            sum(results$stratified_padj < 0.10, na.rm = TRUE), pi0_strat))
cat("INTERPRETATION: pi0 closer to 1 means fewer promoters look genuinely",
    "non-null. A stratified hit count that collapses toward a small,",
    "biologically plausible number (tens, not thousands) is the signature",
    "of the confound being successfully conditioned out.\n")

cat("\n=== Candidates: survive stratification AND pass FDR control ===\n")
candidates <- results[match(screen_candidates(setNames(results$stratified_p, results$promoter)),
                            results$promoter), ]
cat("Candidates at BH q < 0.10:", nrow(candidates), "\n")
if (nrow(candidates) > 0) print(head(candidates, 25))

cat("\nOverlap: promoters significant (raw p<0.05) in BOTH screens vs naive-only:\n")
naive_hits <- results$promoter[!is.na(results$naive_p) & results$naive_p < 0.05]
strat_hits <- results$promoter[!is.na(results$stratified_p) & results$stratified_p < 0.05]
cat("Naive-only hits (subtype/purity likely explains these):",
    length(setdiff(naive_hits, strat_hits)), "\n")
cat("Hits in both (candidates for real, subtype-independent signal):",
    length(intersect(naive_hits, strat_hits)), "\n")

pdf(file.path(paths$diagnostics, "naive_vs_stratified_p_histograms.pdf"),
    width = 10, height = 5)
par(mfrow = c(1, 2))
hist(pi0_naive_p, breaks = 30, col = "steelblue", border = "white",
     xlab = "p-value", main = sprintf("Naive screen (pi0 = %.2f)", pi0_naive))
hist(pi0_strat_p, breaks = 30, col = "tomato", border = "white",
     xlab = "p-value",
     main = sprintf("Subtype+purity stratified (pi0 = %.2f)", pi0_strat))
dev.off()

write.csv(results, file.path(paths$tables, "naive_vs_stratified_screen_full.csv"),
          row.names = FALSE)
write.csv(candidates, file.path(paths$tables, "stratified_screen_candidates_fdr10.csv"),
          row.names = FALSE)

write_manifest(paths$run, list(
  script            = "10_subtype_purity_stratified_screen.R",
  n_promoters       = n_prom,
  n_complete_cases  = sum(complete),
  min_events_per_subtype = min(event_counts),
  naive_hits_p05    = sum(results$naive_p < 0.05, na.rm = TRUE),
  stratified_hits_p05 = sum(results$stratified_p < 0.05, na.rm = TRUE),
  stratified_candidates_fdr10 = nrow(candidates),
  seed              = 42
), inputs = METH_SCREEN_INPUTS)

cat("\n=== Done. Results in", paths$run, "===\n")
