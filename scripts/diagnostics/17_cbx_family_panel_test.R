# =============================================================================
# 17_cbx_family_panel_test.R
#
# Purpose: test whether the CBX (chromobox) gene family adds prognostic
# information beyond subtype + purity, motivated by CBX5 surviving scripts
# 14-15's confound testing.
#
# CIRCULARITY WARNING, read before interpreting any result below:
# CBX5 was DISCOVERED by screening all 4,141 promoters against this exact
# 566-patient cohort in script 14. Any result reported here for CBX5 is
# CONFIRMATORY on the same data it was found in, not independent evidence --
# scripts 14/15's permutation and Her2 checks were legitimate robustness
# checks, but they never left this dataset either. CBX1 and CBX3 are
# genuinely blind: they were never screened, they're included purely because
# they're HP1's other two paralogs (same chromo-shadow domain, same general
# heterochromatin-reading function as CBX5/HP1alpha). Any signal from CBX1/
# CBX3 is real new evidence. Any signal attributed to CBX5 is not -- report
# these two categories separately, not as one undifferentiated "signature."
#
# Two panels, kept separate because they are mechanistically distinct:
#   hp1_paralogs <- CBX1, CBX3, CBX5   (HP1 family: heterochromatin/H3K9me)
#   polycomb_cbx <- CBX2, CBX4, CBX6, CBX7, CBX8   (Polycomb/H3K27me3, PRC1)
#
# With only 85 events total, this script deliberately avoids putting
# multiple genes into one Cox model as separate terms (severe overfitting
# risk). Genes are tested one at a time, and the "signature" is a single
# combined score, not a multi-term model -- the same low-dimensional
# discipline as the mean_cv_risk_score approach in 08_methylation_survival_v3.
#
# Run from the project root, after scripts 14 and 15.
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(readxl)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "cbx_family_panel", base = "02_analysis")
set.seed(42)

N_PERM <- 100


# - - - 1. Load data, gene_mat, subtype, purity - - -
# Identical to scripts 13-15.

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
gene_sums   <- rowsum(m_mat_promoters, group = gene_labels_clean)
gene_counts <- table(gene_labels_clean)
gene_mat    <- gene_sums / as.vector(gene_counts)

pam50_path <- "01_processed_data/PAM50_classification.txt"
subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)
names(subtype_tab) <- c("barcode", "subtype")
subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)
subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]
meth_patient <- substr(meth_clinical$sample_barcode, 1, 12)
meth_clinical$subtype <- factor(subtype_tab$subtype[match(meth_patient, subtype_tab$patient)])

purity_tab <- as.data.frame(read_excel("01_processed_data/aran_purity_tcga.xlsx",
                                       sheet = "Supp Data 1", skip = 2))
names(purity_tab)[names(purity_tab) == "Sample ID"] <- "sample_barcode"
purity_tab$CPE <- as.numeric(purity_tab$CPE)
purity_tab <- purity_tab[purity_tab$`Cancer type` == "BRCA" &
                           !is.na(purity_tab$`Cancer type`), ]
purity_tab$sample_id_15 <- substr(purity_tab$sample_barcode, 1, 15)
meth_clinical$purity <- purity_tab$CPE[
  match(substr(meth_clinical$sample_barcode, 1, 15), purity_tab$sample_id_15)
]

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)
meth_c     <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c     <- surv_obj[complete]

cat("Analysis cohort:", sum(complete), "patients,", sum(meth_c$event), "events\n")


# - - - 2. Define the panels and check promoter coverage - - -

hp1_paralogs <- c("CBX1", "CBX3", "CBX5")
polycomb_cbx <- c("CBX2", "CBX4", "CBX6", "CBX7", "CBX8")

hp1_present      <- intersect(hp1_paralogs, rownames(gene_mat_c))
polycomb_present <- intersect(polycomb_cbx, rownames(gene_mat_c))

cat("\nHP1 paralogs with a promoter in gene_mat:", paste(hp1_present, collapse = ", "), "\n")
cat("Missing:", paste(setdiff(hp1_paralogs, hp1_present), collapse = ", "), "\n")
cat("\nPolycomb CBX members with a promoter in gene_mat:",
    paste(polycomb_present, collapse = ", "), "\n")
cat("Missing:", paste(setdiff(polycomb_cbx, polycomb_present), collapse = ", "), "\n")


# - - - 3. Screening function (same as scripts 14-15) - - -

stratified_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity) + strata(meth_x$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}


# - - - 4. Individual-gene tests (the actual blind test, for CBX1/CBX3) - - -

cat("\n=== Individual genes, subtype+purity-stratified ===\n")

all_panel_genes <- unique(c(hp1_present, polycomb_present))
indiv_res <- t(sapply(all_panel_genes, function(g) {
  stratified_screen_one(gene_mat_c[g, ], surv_c, meth_c)
}))
indiv_res <- data.frame(
  gene   = all_panel_genes,
  coef   = indiv_res[, "coef"],
  p      = indiv_res[, "p"],
  status = ifelse(all_panel_genes == "CBX5",
                   "CONFIRMATORY (discovered in this data via script 14)",
                   "BLIND (never screened before this run)"),
  stringsAsFactors = FALSE
)
indiv_res$padj_within_panel <- p.adjust(indiv_res$p, method = "BH")
print(indiv_res[order(indiv_res$p), ])


# - - - 5. Combined score: HP1 paralogs only - - -
# Simple, transparent, low-dimensional: mean of z-scored promoter values
# across the panel, entered as ONE covariate. Deliberately not a LASSO-
# weighted combination -- with only 3 genes and 85 events, a weighted fit
# would just be overfitting noise into the weights.

cat("\n=== Combined HP1-paralog score ===\n")

if (length(hp1_present) >= 2) {
  # rowMeans, not colMeans: see the note in script 20 -- t(...) is patients x
  # genes, scale() z-scores each gene (column), and averaging ACROSS genes
  # per patient needs rowMeans. colMeans would collapse to one number per
  # GENE instead of one score per patient (the same bug already caught once
  # in script 13's B1 proliferation section).
  hp1_score <- rowMeans(scale(t(gene_mat_c[hp1_present, , drop = FALSE])))
  stopifnot("hp1_score must have one value per patient" =
              length(hp1_score) == ncol(gene_mat_c))
  meth_c$hp1_score <- hp1_score

  r <- stratified_screen_one(hp1_score, surv_c, meth_c)
  cat(sprintf("HP1 combined score: coef=%.3f, p=%.3g\n", r["coef"], r["p"]))
  cat("NOTE: this score includes CBX5, so its apparent strength is not fully",
      "independent of the data it was mined from -- compare against the",
      "CBX1/CBX3-only score below to see what the blind genes contribute",
      "on their own.\n")

  blind_only <- setdiff(hp1_present, "CBX5")
  if (length(blind_only) >= 1) {
    blind_score <- if (length(blind_only) == 1) {
      scale(gene_mat_c[blind_only, ])[, 1]
    } else {
      rowMeans(scale(t(gene_mat_c[blind_only, , drop = FALSE])))
    }
    stopifnot("blind_score must have one value per patient" =
                length(blind_score) == ncol(gene_mat_c))
    meth_c$blind_score <- blind_score
    r_blind <- stratified_screen_one(blind_score, surv_c, meth_c)
    cat(sprintf("\nCBX1+CBX3 ONLY (fully blind) score: coef=%.3f, p=%.3g\n",
                r_blind["coef"], r_blind["p"]))
    cat("THIS is the genuinely independent result -- no gene in this score",
        "was selected using this outcome data.\n")
  }
} else {
  cat("Fewer than 2 HP1 paralogs present -- skipping combined score.\n")
}


# - - - 6. Permutation validation, rebuilt for this panel's scale - - -
# Same logic as script 15 Test A, resized: shuffle methylation<->patient
# assignment, rerun the panel-only screen + BH correction, count hits.

cat("\n=== Permutation null, rescaled to this panel (", N_PERM, "perms) ===\n")

obs_padj   <- indiv_res$padj_within_panel
obs_n_hits <- sum(obs_padj < 0.10, na.rm = TRUE)
obs_min_p  <- min(indiv_res$p, na.rm = TRUE)

perm_n_hits <- integer(N_PERM)
perm_min_p  <- numeric(N_PERM)

for (r in seq_len(N_PERM)) {
  set.seed(2000 + r)
  perm_idx <- sample(ncol(gene_mat_c))
  gm_perm  <- gene_mat_c[all_panel_genes, perm_idx, drop = FALSE]
  res_p <- t(sapply(all_panel_genes, function(g) {
    stratified_screen_one(gm_perm[g, ], surv_c, meth_c)
  }))
  padj_p <- p.adjust(res_p[, "p"], method = "BH")
  perm_n_hits[r] <- sum(padj_p < 0.10, na.rm = TRUE)
  perm_min_p[r]  <- min(res_p[, "p"], na.rm = TRUE)
}

perm_p_hitcount <- (sum(perm_n_hits >= obs_n_hits) + 1) / (N_PERM + 1)
perm_p_minp     <- (sum(perm_min_p <= obs_min_p) + 1) / (N_PERM + 1)

cat(sprintf("Observed hits (q<0.10, within %d-gene panel): %d\n",
            length(all_panel_genes), obs_n_hits))
cat(sprintf("Permutation p (hit count): %.3f | Permutation p (min p): %.3f\n",
            perm_p_hitcount, perm_p_minp))
cat("NOTE: with only", length(all_panel_genes), "genes tested, BH correction",
    "is far less severe than in the 4,141-promoter screen -- do not expect,",
    "or need, the same q-value thresholds to feel as strict.\n")


# - - - 7. Her2-exclusion stability, per gene - - -

cat("\n=== Her2-exclusion stability ===\n")

no_her2 <- meth_c$subtype != "Her2"
meth_no_her2 <- meth_c[no_her2, ]
meth_no_her2$subtype <- droplevels(meth_no_her2$subtype)
gene_mat_no_her2 <- gene_mat_c[, no_her2, drop = FALSE]
surv_no_her2 <- surv_c[no_her2]

stability <- data.frame(
  gene = all_panel_genes, p_with_her2 = indiv_res$p,
  p_without_her2 = NA_real_, same_direction = NA
)
for (g in all_panel_genes) {
  r <- stratified_screen_one(gene_mat_no_her2[g, ], surv_no_her2, meth_no_her2)
  i <- which(stability$gene == g)
  stability$p_without_her2[i] <- r["p"]
  stability$same_direction[i] <- sign(indiv_res$coef[indiv_res$gene == g]) == sign(r["coef"])
}
print(stability)


# - - - 8. Save outputs - - -

write.csv(indiv_res, file.path(paths$tables, "cbx_family_individual_results.csv"),
          row.names = FALSE)
write.csv(stability, file.path(paths$tables, "cbx_family_her2_stability.csv"),
          row.names = FALSE)
write.csv(data.frame(permutation = seq_len(N_PERM), n_hits = perm_n_hits,
                     min_p = perm_min_p),
          file.path(paths$tables, "cbx_family_permutation_null.csv"),
          row.names = FALSE)

write_manifest(paths$run, list(
  script = "17_cbx_family_panel_test.R",
  hp1_genes_present = paste(hp1_present, collapse = ";"),
  polycomb_genes_present = paste(polycomb_present, collapse = ";"),
  n_permutations = N_PERM,
  observed_hits = obs_n_hits,
  permutation_p_hitcount = round(perm_p_hitcount, 4),
  seed = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
cat("Remember: CBX5 results are confirmatory, not independent. The blind",
    "test is CBX1/CBX3 (and the Polycomb genes, if present) considered on",
    "their own.\n")
