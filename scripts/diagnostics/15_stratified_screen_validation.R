# =============================================================================
# 15_stratified_screen_validation.R
#
# Purpose: script 14 produced 5 candidate promoters that survive subtype+purity
# stratification at BH q<0.10 (C1orf86, PLA2G3, CBX5, KLK9, RAB6B), out of
# 4,141 tests, in a 566-patient / 85-event cohort split across four uneven
# strata (Her2 has only 6 events). That is an exploratory hit list, not a
# validated one. This script runs two checks before trusting it:
#
#   TEST A (permutation / negative control, D5-style): if promoter identity
#   is scrambled relative to patient outcome -- breaking any real
#   promoter-survival relationship while leaving the real, established
#   age/stage/purity/subtype-survival relationships untouched -- how many
#   promoters cross BH q<0.10 by chance alone, using the EXACT same
#   screen-and-correct procedure as the real run? If 5 hits is unremarkable
#   under this null, the candidate list does not clear the bar.
#
#   TEST B (Her2-exclusion stability): refit each of the 5 candidates with
#   Her2 patients dropped. A real, broadly-supported effect should survive
#   losing 33 patients / 6 events. A candidate whose significance depends on
#   the Her2 stratum specifically is not a stable finding.
#
# Run from the project root, after 14_subtype_purity_stratified_screen.R.
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(readxl)
library(doParallel)
library(foreach)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "stratified_screen_validation", base = "02_analysis")
set.seed(42)

N_PERM <- 100   # ~12 sec/permutation single-core; parallelized below.
                # 100 gives permutation p-values to a resolution of 0.01.
                # Drop to 20-30 for a quick sanity check before committing
                # to a full run.


# - - - 1. Load methylation data and rebuild gene_mat - - -
# Identical to scripts 13 and 14.

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


# - - - 2. Load PAM50 subtype and purity - - -
# Identical to script 14.

pam50_path <- "01_processed_data/PAM50_classification.txt"
subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)
names(subtype_tab) <- c("barcode", "subtype")
subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)
subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]

meth_patient <- substr(meth_clinical$sample_barcode, 1, 12)
meth_clinical$subtype <- factor(subtype_tab$subtype[match(meth_patient, subtype_tab$patient)])

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


# - - - 3. Define the analysis cohort - - -
# Identical to script 14 -- same 566 patients, so results are comparable.

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)

meth_c     <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c     <- surv_obj[complete]

cat("Analysis cohort:", sum(complete), "patients,",
    sum(meth_c$event), "events\n")


# - - - 4. Screening function (explicit args, not closures) - - -
# Takes the promoter value vector AND the survival/covariate data as
# arguments, rather than reading meth_c/surv_c from the enclosing
# environment. This matters for TEST A: the whole point of the permutation
# is to scramble which patient's row of gene_mat gets compared against
# which patient's survival/covariates, so the function needs to accept
# whatever (g, surv_x, meth_x) triple it's given, real or permuted.

stratified_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity) + strata(meth_x$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

run_full_screen <- function(gm, surv_x, meth_x) {
  res <- t(apply(gm, 1, stratified_screen_one, surv_x = surv_x, meth_x = meth_x))
  colnames(res) <- c("coef", "z", "p")
  res
}


# - - - 5. Recompute the observed (real) screen - - -
# Recomputed fresh here (not loaded from script 14's output) so this
# validation is guaranteed to be testing the exact same code path.

cat("\nRunning observed stratified screen on", nrow(gene_mat_c), "promoters...\n")
t0 <- Sys.time()
obs_res <- run_full_screen(gene_mat_c, surv_c, meth_c)
cat("Done in", round(as.numeric(Sys.time() - t0, units = "secs")), "sec\n")

obs_padj    <- p.adjust(obs_res[, "p"], method = "BH")
obs_n_hits  <- sum(obs_padj < 0.10, na.rm = TRUE)
obs_min_p   <- min(obs_res[, "p"], na.rm = TRUE)

cat("Observed: ", obs_n_hits, "promoters at BH q<0.10\n")


# =============================================================================
# TEST A: Permutation null
# =============================================================================
# Shuffle which patient's methylation profile (columns of gene_mat_c) gets
# compared against which patient's survival/covariates. This breaks any real
# promoter-survival relationship while leaving age, stage, purity, and
# subtype correctly paired with their true outcomes -- so the established
# real confound effects stay intact, and only the promoter term is forced
# to be null. Re-run the IDENTICAL screen-and-BH-correct procedure used on
# the real data, and count how many promoters cross q<0.10 by chance.

cat("\n=== TEST A: permutation null (", N_PERM, "permutations) ===\n")

n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Running on", n_cores, "cores...\n")

t0 <- Sys.time()
perm_results <- foreach(
  r = seq_len(N_PERM),
  .packages = c("survival"),
  .combine  = rbind
) %dopar% {
  set.seed(1000 + r)
  perm_idx <- sample(ncol(gene_mat_c))
  gm_perm  <- gene_mat_c[, perm_idx, drop = FALSE]
  # meth_c and surv_c are NOT touched -- only the gene_mat <-> patient
  # pairing is scrambled, so age/stage/purity/subtype keep their real,
  # true relationship to survival.

  res_p  <- run_full_screen(gm_perm, surv_c, meth_c)
  padj_p <- p.adjust(res_p[, "p"], method = "BH")

  c(n_hits = sum(padj_p < 0.10, na.rm = TRUE),
    min_p  = min(res_p[, "p"], na.rm = TRUE))
}
stopCluster(cl)
cat("Permutation run:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

perm_results <- as.data.frame(perm_results)

# Permutation p-value: how often does chance alone produce a hit count at
# least as extreme as what was actually observed? The "+1" in numerator and
# denominator is standard practice -- it accounts for the observed data
# itself being one valid draw from the null, and guarantees the p-value is
# never reported as exactly zero regardless of N_PERM.
perm_p_hitcount <- (sum(perm_results$n_hits >= obs_n_hits) + 1) / (N_PERM + 1)
perm_p_minp     <- (sum(perm_results$min_p <= obs_min_p) + 1) / (N_PERM + 1)

cat("\nNull distribution of hit counts (BH q<0.10) across", N_PERM, "permutations:\n")
print(summary(perm_results$n_hits))
cat(sprintf("\nObserved: %d hits. Permutation p-value: %.3f\n",
            obs_n_hits, perm_p_hitcount))
cat(sprintf("Observed smallest p-value: %.3g. Permutation p-value: %.3f\n",
            obs_min_p, perm_p_minp))
cat("INTERPRETATION: a small perm_p_hitcount (e.g. <0.05) means the real",
    "screen finds more hits than chance alone typically produces under this",
    "exact procedure -- the candidate list clears the bar. If perm_p_hitcount",
    "is large, 5 hits is unremarkable noise and the candidates should not be",
    "trusted without further validation.\n")

pdf(file.path(paths$diagnostics, "permutation_null_hit_counts.pdf"),
    width = 7, height = 5)
hist(perm_results$n_hits, breaks = max(10, N_PERM / 5), col = "grey80",
     border = "white", xlab = "Promoters at BH q<0.10 (null permutations)",
     main = sprintf("Permutation null (n=%d) vs observed (perm p=%.3f)",
                     N_PERM, perm_p_hitcount),
     xlim = range(c(perm_results$n_hits, obs_n_hits)))
abline(v = obs_n_hits, col = "tomato", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Observed = %d", obs_n_hits),
       col = "tomato", lty = 2, lwd = 2, bty = "n")
dev.off()


# =============================================================================
# TEST B: Her2-exclusion stability
# =============================================================================
# Her2 carries only 6 events in this cohort (see script 14's warning). Any
# candidate whose apparent significance depends on that stratum specifically
# is not a stable finding. Refit each candidate with Her2 dropped entirely.

cat("\n=== TEST B: stability with Her2 excluded ===\n")

candidates <- c("C1orf86", "PLA2G3", "CBX5", "KLK9", "RAB6B")
missing_candidates <- setdiff(candidates, rownames(gene_mat_c))
if (length(missing_candidates) > 0) {
  cat("WARNING: these candidates are not in gene_mat_c and will be skipped",
      "(check for a promoter-naming mismatch):",
      paste(missing_candidates, collapse = ", "), "\n")
}
candidates <- intersect(candidates, rownames(gene_mat_c))

no_her2 <- meth_c$subtype != "Her2"
meth_no_her2     <- meth_c[no_her2, ]
meth_no_her2$subtype <- droplevels(meth_no_her2$subtype)
gene_mat_no_her2 <- gene_mat_c[, no_her2, drop = FALSE]
surv_no_her2     <- surv_c[no_her2]

cat("Cohort with Her2 excluded:", sum(no_her2), "patients,",
    sum(meth_no_her2$event), "events\n")

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

# Same-direction coefficient AND still nominally significant without Her2 is
# the bar for "not obviously Her2-dependent" -- NOT full FDR control again,
# since we're now testing only 5 pre-selected promoters, not screening 4,141.
stability$same_direction <- sign(stability$coef_with_her2) ==
  sign(stability$coef_without_her2)
stability$still_p05 <- stability$p_without_her2 < 0.05

cat("\n")
print(stability)
cat("\nINTERPRETATION: same_direction = FALSE or still_p05 = FALSE flags a",
    "candidate whose result is not stable to removing the Her2 stratum --",
    "treat it as substantially weaker evidence than its q-value alone",
    "suggests.\n")


# - - - 6. Save outputs - - -

write.csv(data.frame(permutation = seq_len(N_PERM), perm_results),
          file.path(paths$tables, "permutation_null_results.csv"),
          row.names = FALSE)
write.csv(stability, file.path(paths$tables, "her2_exclusion_stability.csv"),
          row.names = FALSE)

write_manifest(paths$run, list(
  script                = "15_stratified_screen_validation.R",
  n_permutations        = N_PERM,
  observed_hits_fdr10   = obs_n_hits,
  permutation_p_hitcount = round(perm_p_hitcount, 4),
  permutation_p_minp    = round(perm_p_minp, 4),
  n_candidates_tested   = length(candidates),
  n_stable_both_checks  = sum(stability$same_direction & stability$still_p05),
  seed                  = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
