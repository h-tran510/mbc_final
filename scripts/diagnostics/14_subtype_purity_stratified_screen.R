# =============================================================================
# 14_subtype_purity_stratified_screen.R
#
# Purpose: extend B2's single-axis test (does PC1 survive strata(subtype)?)
# into a per-promoter screen. B2 showed the DOMINANT shared axis is largely
# explained by subtype + purity. That does not mean no individual promoter
# carries real, non-redundant prognostic information -- it means the naive,
# unstratified screen cannot tell the difference between "this promoter is
# real biology" and "this promoter is just riding on subtype/purity like
# thousands of others." This script reruns the univariate screen twice --
# once naively (matching uni_screen() in 00_utils.R), once conditioned on
# subtype and purity -- on the SAME patients, and compares:
#   1. does the D10 sign bias (originally ~2:1 negative) flatten toward 50/50?
#   2. does the D11 hit-count inflation (~1,415 "significant" promoters) drop
#      toward something small and plausible?
#   3. which specific promoters, if any, survive conditioning -- these are
#      the real candidates for "independent of the confound," not PC1 itself.
#
# IMPORTANT CAVEAT, read before trusting any single promoter's result:
# strata(subtype) needs enough EVENTS within each subtype, not just enough
# patients, for a stable estimate. Her2 has only 33 patients in this cohort --
# check the per-subtype event count printed below before interpreting.
#
# BEFORE RUNNING: write your prediction -- do you expect the sign bias to
# flatten, and roughly how many promoters (if any) you'd expect to survive --
# in your notes, with today's date and commit hash.
#
# Run from the project root.
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(readxl)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "subtype_purity_stratified_screen", base = "02_analysis")
set.seed(42)


# - - - 1. Load methylation data and rebuild gene_mat - - -
# Identical to 13_pc1_biological_characterization.R sections 1 and 1c.

m_mat         <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

ann_matched <- cache_or_compute("01_processed_data/ann_matched_cache.rds", function() {
  ann_df <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  ann_df[match(rownames(m_mat), rownames(ann_df)), ]
})

is_promoter     <- grepl("TSS200|TSS1500", ann_matched$UCSC_RefGene_Group)
m_mat_promoters <- m_mat[is_promoter, ]
gene_labels     <- ann_matched$UCSC_RefGene_Name[is_promoter]

gene_labels_clean <- sapply(strsplit(gene_labels, split = ";"), function(x) x[1])
gene_sums         <- rowsum(m_mat_promoters, group = gene_labels_clean)
gene_counts       <- table(gene_labels_clean)
gene_mat          <- gene_sums / as.vector(gene_counts)

cat("Promoters:", nrow(gene_mat), " Samples:", ncol(gene_mat), "\n")


# - - - 2. Load PAM50 subtype - - -
# Identical to 13_pc1_biological_characterization.R section B2.

pam50_path <- "01_processed_data/PAM50_classification.txt"
subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("PAM50 file should have 2 columns" = ncol(subtype_tab) == 2)
names(subtype_tab) <- c("barcode", "subtype")
subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)

dup_patients <- unique(subtype_tab$patient[duplicated(subtype_tab$patient)])
if (length(dup_patients) > 0) {
  disagree <- vapply(dup_patients, function(p) {
    length(unique(subtype_tab$subtype[subtype_tab$patient == p])) > 1
  }, logical(1))
  if (any(disagree)) {
    cat("CONFLICTING PAM50 calls -- inspect before trusting these patients:\n")
    print(subtype_tab[subtype_tab$patient %in% dup_patients[disagree], ])
  }
}
subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]

meth_patient <- substr(meth_clinical$sample_barcode, 1, 12)
meth_clinical$subtype <- factor(subtype_tab$subtype[match(meth_patient, subtype_tab$patient)])


# - - - 3. Load tumor purity - - -
# Identical to 13_pc1_biological_characterization.R section B3.

purity_path <- "01_processed_data/aran_purity_tcga.xlsx"
purity_tab  <- as.data.frame(read_excel(purity_path, sheet = "Supp Data 1", skip = 2))
names(purity_tab)[names(purity_tab) == "Sample ID"] <- "sample_barcode"
purity_tab$CPE <- as.numeric(purity_tab$CPE)
purity_tab <- purity_tab[purity_tab$`Cancer type` == "BRCA" &
                           !is.na(purity_tab$`Cancer type`), ]
purity_tab$sample_id_15 <- substr(purity_tab$sample_barcode, 1, 15)

meth_clinical$purity <- purity_tab$CPE[
  match(substr(meth_clinical$sample_barcode, 1, 15), purity_tab$sample_id_15)
]


# - - - 4. Define the analysis cohort - - -
# Both screens run on IDENTICAL patients, so the before/after comparison is
# apples to apples. The stratified screen additionally requires subtype.

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)

cat("\nComplete cases (subtype + purity + age + stage):",
    sum(complete), "of", nrow(meth_clinical), "\n")

meth_c   <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c   <- surv_obj[complete]

# CRITICAL CHECK: events per subtype, not just patients per subtype. A
# stratum with few events makes that subtype's contribution to the partial
# likelihood unstable, regardless of how many patients sit in it.
cat("\nPatients and EVENTS per subtype (check before trusting per-promoter results):\n")
print(table(subtype = meth_c$subtype, event = meth_c$event))
event_counts <- tapply(meth_c$event, meth_c$subtype, sum)
if (any(event_counts < 10)) {
  cat("\nWARNING: at least one subtype has fewer than 10 events.",
      "Coefficients for promoters whose signal is concentrated in that",
      "stratum will be unstable -- treat individual hits from that group",
      "with caution regardless of p-value.\n")
}


# - - - 5. Define the two screens - - -

# Naive: matches uni_screen() in 00_utils.R exactly (no subtype term).
naive_screen_one <- function(g) {
  fit <- try(coxph(surv_c ~ scale(g) + meth_c$age_years + meth_c$stage +
                      scale(meth_c$purity)), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

# Stratified: adds strata(subtype), restricting every risk-set comparison to
# within-subtype. This is B2's logic, applied to ONE promoter at a time
# instead of to PC1.
stratified_screen_one <- function(g) {
  fit <- try(coxph(surv_c ~ scale(g) + meth_c$age_years + meth_c$stage +
                      scale(meth_c$purity) + strata(meth_c$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}


# - - - 6. Run both screens across all promoters - - -
# Sequential on purpose: ~4,100 simple Cox fits complete in a few minutes and
# this is a one-off diagnostic, not production infrastructure. If it drags,
# the screen_train()-style foreach/doParallel pattern from script 12 drops in
# cleanly here -- each promoter's fit is independent of every other's.

n_prom <- nrow(gene_mat_c)
cat("\nRunning naive screen on", n_prom, "promoters...\n")
t0 <- Sys.time()
naive_res <- t(apply(gene_mat_c, 1, naive_screen_one))
cat("Naive screen:", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

cat("\nRunning subtype+purity-stratified screen on", n_prom, "promoters...\n")
t0 <- Sys.time()
strat_res <- t(apply(gene_mat_c, 1, stratified_screen_one))
cat("Stratified screen:", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

results <- data.frame(
  promoter        = rownames(gene_mat_c),
  naive_coef      = naive_res[, "coef"],
  naive_p         = naive_res[, "p"],
  stratified_coef = strat_res[, "coef"],
  stratified_p    = strat_res[, "p"],
  stringsAsFactors = FALSE
)

# Multiple testing: raw p < 0.05 across ~4,100 tests is not a defensible
# "significant" claim on its own -- this is the same issue D11 flagged.
results$naive_padj      <- p.adjust(results$naive_p, method = "BH")
results$stratified_padj <- p.adjust(results$stratified_p, method = "BH")


# - - - 7. Does the D10/D11 anomaly resolve? - - -

cat("\n=== Sign bias (D10 analog) ===\n")
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

cat("\n=== Hit count / pi0 (D11 analog) ===\n")
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
candidates <- results[!is.na(results$stratified_padj) &
                         results$stratified_padj < 0.10, ]
candidates <- candidates[order(candidates$stratified_padj), ]
cat("Candidates at BH q < 0.10:", nrow(candidates), "\n")
if (nrow(candidates) > 0) print(head(candidates, 25))

cat("\nOverlap: promoters significant (raw p<0.05) in BOTH screens vs naive-only:\n")
naive_hits <- results$promoter[!is.na(results$naive_p) & results$naive_p < 0.05]
strat_hits <- results$promoter[!is.na(results$stratified_p) & results$stratified_p < 0.05]
cat("Naive-only hits (subtype/purity likely explains these):",
    length(setdiff(naive_hits, strat_hits)), "\n")
cat("Hits in both (candidates for real, subtype-independent signal):",
    length(intersect(naive_hits, strat_hits)), "\n")


# - - - 8. Diagnostic plot - - -

pdf(file.path(paths$diagnostics, "naive_vs_stratified_p_histograms.pdf"),
    width = 10, height = 5)
par(mfrow = c(1, 2))
hist(pi0_naive_p, breaks = 30, col = "steelblue", border = "white",
     xlab = "p-value", main = sprintf("Naive screen (pi0 = %.2f)", pi0_naive))
hist(pi0_strat_p, breaks = 30, col = "tomato", border = "white",
     xlab = "p-value",
     main = sprintf("Subtype+purity stratified (pi0 = %.2f)", pi0_strat))
dev.off()


# - - - 9. Save outputs - - -

write.csv(results, file.path(paths$tables, "naive_vs_stratified_screen_full.csv"),
          row.names = FALSE)
write.csv(candidates, file.path(paths$tables, "stratified_screen_candidates_fdr10.csv"),
          row.names = FALSE)

write_manifest(paths$run, list(
  script            = "14_subtype_purity_stratified_screen.R",
  n_promoters       = n_prom,
  n_complete_cases  = sum(complete),
  min_events_per_subtype = min(event_counts),
  naive_hits_p05    = sum(results$naive_p < 0.05, na.rm = TRUE),
  stratified_hits_p05 = sum(results$stratified_p < 0.05, na.rm = TRUE),
  stratified_candidates_fdr10 = nrow(candidates),
  seed              = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
