# =============================================================================
# 13_pc1_biological_characterization.R
#
# Purpose: determine WHAT PC1 (the dominant shared factor in promoter
# methylation) actually is, biologically. PC1 is already known to be:
#   - pervasive across promoters (drives the inflated univariate hit count)
#   - independently associated with PFI
#   - correlated with global mean methylation (r = -0.711 in the v3 run)
#
# This script does not build a predictive model. It tests competing
# explanations for PC1:
#   B1  proliferation           (RNA-seq)      <- highest priority
#   B2  PAM50 subtype
#   B3  tumor purity            (Aran et al. 2015 CPE)
#   B4  genomic context of PC1  (island vs open sea)  <- highest priority
#   B5  curated CIMP panel      (Wang et al. 2017)    <- needs manual input
#   B6  mutation burden / TP53
#   B7  integration: what is left of PC1 once B1-B6 are accounted for
#
# BEFORE RUNNING: write your predictions for each test in your project notes,
# with today's date and the current commit hash. Do not revise them afterward;
# append instead, so the gap between prediction and result stays visible.
#
# Run from the project root (same working directory as 08_methylation_survival_v3.R).
# =============================================================================


# - - - 0. Setup - - -

library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(readxl)

source("scripts/00_utils.R")

paths <- make_run_dir(tag = "pc1_B1_fixed", base = "02_analysis")

# Toggles for the expensive / externally-dependent sections
run_b6_mutations <- FALSE   # requires a TCGA mutation download (see B6)
set.seed(42)

# Collects one row per test so everything lands in a single results table
results_log <- list()
log_result <- function(test, metric, value, note = "") {
  # unname(): coefficient extractions carry names like "scale(pc1)", which
  # otherwise leak into the results table as row names
  results_log[[length(results_log) + 1]] <<- data.frame(
    test = test, metric = metric, value = unname(value), note = note,
    stringsAsFactors = FALSE, row.names = NULL
  )
}


# - - - 1. Load methylation data and rebuild gene_mat - - -
# Identical to sections 1 and 1c of 08_methylation_survival_v3.R, so PC1 here
# is the same object as in the pipeline. If you change aggregation there,
# change it here too (or factor this block into 00_utils.R).

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


# - - - 1b. Compute PC1 and fix its sign - - -
# PCA sign is arbitrary: prcomp can return PC1 or -PC1 depending on numerical
# details, which silently flips every hazard ratio and correlation below.
# Anchor it to global mean methylation so "high PC1 = more methylated" always.

pca <- prcomp(t(gene_mat), scale. = TRUE)
pc1 <- pca$x[, 1]

global_meth <- colMeans(gene_mat)

if (cor(pc1, global_meth) < 0) {
  pc1 <- -pc1
  pca$rotation[, 1] <- -pca$rotation[, 1]
  cat("PC1 sign flipped so that higher PC1 = higher global methylation.\n")
}

var_explained <- summary(pca)$importance[2, 1]
cat(sprintf("PC1 variance explained: %.1f%%\n", 100 * var_explained))
cat(sprintf("cor(PC1, global mean methylation) = %.3f\n", cor(pc1, global_meth)))

log_result("B0", "pc1_var_explained", var_explained)
log_result("B0", "cor_pc1_globalmeth", cor(pc1, global_meth))

meth_clinical$pc1 <- pc1
meth_clinical$global_meth <- global_meth


# =============================================================================
# B1. Is PC1 a proliferation readout?
# =============================================================================
# Venet & Boutros (2011) showed random gene signatures look prognostic in breast
# cancer because proliferation is pervasive across the transcriptome. If PC1 is
# tightly correlated with proliferation, your methylation confound is the same
# phenomenon surfacing in a second molecular layer. If it is NOT, methylation
# has its own distinct pervasive confound -- the more novel result.

cat("\n=== B1. PC1 vs proliferation ===\n")

expr_path <- "01_processed_data/expr_matrix_vst.rds"

if (file.exists(expr_path)) {
  expr_matrix <- readRDS(expr_path)

  # Rownames may be Ensembl IDs (ENSG..., possibly versioned) or gene symbols.
  # Build a symbol -> row index lookup that works for either.
  expr_rownames <- rownames(expr_matrix)

  if (any(grepl("^ENSG", expr_rownames))) {
    cat("Expression rownames look like Ensembl IDs; mapping to symbols.\n")
    brca_se  <- readRDS("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds")
    rd       <- as.data.frame(SummarizedExperiment::rowData(brca_se))
    # gene_name is the standard column from GDCprepare; fall back if renamed
    symbol_col <- if ("gene_name" %in% names(rd)) "gene_name" else names(rd)[1]
    ens_strip  <- sub("\\..*$", "", expr_rownames)
    rd_strip   <- sub("\\..*$", "", rownames(rd))
    sym_lookup <- rd[[symbol_col]][match(ens_strip, rd_strip)]
    rm(brca_se); gc()
  } else {
    sym_lookup <- expr_rownames
  }

  # Proliferation panel. MKI67 alone is the crude version; the panel is more
  # stable. Alternative: pull HALLMARK_E2F_TARGETS via msigdbr::msigdbr().
  prolif_genes <- c("MKI67", "PCNA", "TOP2A", "BIRC5", "CCNB1", "AURKA", "RRM2")
  prolif_idx   <- which(sym_lookup %in% prolif_genes)

  cat("Proliferation genes found:", length(prolif_idx), "of", length(prolif_genes), "\n")

  if (length(prolif_idx) >= 3) {
    # Align expression samples to methylation samples via the 15-char barcode.
    # RNA aliquots carry analyte code R, methylation carries D, so the full
    # 28-char barcodes will never match -- truncation is required, not optional.
    expr_id15 <- substr(colnames(expr_matrix), 1, 15)
    meth_id15 <- substr(meth_clinical$sample_barcode, 1, 15)
    expr_col  <- match(meth_id15, expr_id15)

    cat("Patients with both methylation and expression:",
        sum(!is.na(expr_col)), "of", nrow(meth_clinical), "\n")

    prolif_sub <- expr_matrix[prolif_idx, , drop = FALSE]

    # t() gives samples x genes. scale() then z-scores each GENE (column).
    # We want one score per SAMPLE, so average across genes = rowMeans.
    # BUG FIXED (2026-08-18): this was colMeans, which averages across SAMPLES
    # and returns 7 numbers (one per gene). Indexing that length-7 vector by
    # sample position yielded ~all NA, so the first run silently computed the
    # correlation on 6 patients instead of 755.
    prolif_score <- rowMeans(scale(t(prolif_sub)))
    stopifnot("prolif_score must have one value per sample" =
                length(prolif_score) == ncol(expr_matrix))

    meth_clinical$prolif <- prolif_score[expr_col]        # NA where unmatched

    ok <- !is.na(meth_clinical$prolif)
    cat("Patients contributing to the proliferation correlation:", sum(ok), "\n")

    ct <- cor.test(meth_clinical$pc1[ok], meth_clinical$prolif[ok],
                   method = "spearman")
    print(ct)
    log_result("B1", "spearman_rho_pc1_prolif", unname(ct$estimate),
               sprintf("p = %.3g, n = %d", ct$p.value, sum(ok)))

    # MKI67 alone, as a sanity check that the panel is not doing something odd
    mki_idx <- which(sym_lookup == "MKI67")
    if (length(mki_idx) == 1) {
      mki <- expr_matrix[mki_idx, expr_col]
      ct_mki <- cor.test(meth_clinical$pc1[ok], mki[ok], method = "spearman")
      cat(sprintf("MKI67 alone: rho = %.3f, p = %.3g\n",
                  ct_mki$estimate, ct_mki$p.value))
      log_result("B1", "spearman_rho_pc1_MKI67", unname(ct_mki$estimate))
    }

    pdf(file.path(paths$diagnostics, "b1_pc1_vs_proliferation.pdf"),
        width = 6, height = 5)
    plot(meth_clinical$prolif[ok], meth_clinical$pc1[ok],
         xlab = "Proliferation score (mean z-scored expression)",
         ylab = "PC1 (higher = more methylated)",
         main = sprintf("PC1 vs proliferation (rho = %.2f)", ct$estimate),
         pch = 16, col = rgb(0, 0, 0, 0.35))
    lines(lowess(meth_clinical$prolif[ok], meth_clinical$pc1[ok]),
          col = "red", lwd = 2)
    dev.off()

    # INTERPRETATION
    # |rho| > 0.5  -> PC1 is largely the proliferation confound in another layer
    # |rho| < 0.2  -> methylation carries its own distinct pervasive confound
    # in between   -> partial; quantify with the variance partition in B7

    rm(expr_matrix); gc()
  } else {
    cat("Too few proliferation genes matched -- check the symbol mapping.\n")
  }
} else {
  cat("No expression matrix at", expr_path, "-- skipping B1.\n")
}


# =============================================================================
# B2. Is PC1 just PAM50 subtype?
# =============================================================================

cat("\n=== B2. PC1 vs PAM50 subtype ===\n")

pam50_path <- "01_processed_data/PAM50_classification.txt"

subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)

# The header in this file reads "Particpant ID" (sic -- typo in the source), so
# do NOT match column names by pattern. Take the two columns positionally.
stopifnot("PAM50 file should have 2 columns" = ncol(subtype_tab) == 2)
names(subtype_tab) <- c("barcode", "subtype")

# Barcodes here are 16-char sample+vial (TCGA-A1-A0SK-01A), not the 15-char
# sample barcode used elsewhere in this project. PAM50 is assigned per patient,
# so collapse to the 12-char participant barcode and join on that -- this also
# sidesteps the 01A/01B vial-letter mismatch.
subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)

cat("PAM50 rows:", nrow(subtype_tab), "| unique patients:",
    length(unique(subtype_tab$patient)), "\n")

# A handful of participants appear twice (01A and 01B vials of the same tumor).
# Check the duplicate calls AGREE before collapsing -- a disagreement means the
# subtype call is unstable for that tumor, which is worth knowing rather than
# silently resolving by taking whichever row came first.
dup_patients <- unique(subtype_tab$patient[duplicated(subtype_tab$patient)])
if (length(dup_patients) > 0) {
  disagree <- vapply(dup_patients, function(p) {
    length(unique(subtype_tab$subtype[subtype_tab$patient == p])) > 1
  }, logical(1))
  cat("Duplicate participants:", length(dup_patients),
      "| with conflicting subtype calls:", sum(disagree), "\n")
  if (any(disagree)) {
    cat("CONFLICTING CALLS -- inspect before trusting these patients:\n")
    print(subtype_tab[subtype_tab$patient %in% dup_patients[disagree], ])
  }
}
subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]

meth_patient <- substr(meth_clinical$sample_barcode, 1, 12)
meth_clinical$subtype <- factor(
  subtype_tab$subtype[match(meth_patient, subtype_tab$patient)]
)

cat("Methylation patients with a PAM50 call:",
    sum(!is.na(meth_clinical$subtype)), "of", nrow(meth_clinical), "\n")

cat("\nSubtype counts (check these BEFORE trusting any per-subtype test):\n")
print(table(meth_clinical$subtype, useNA = "ifany"))

sub_ok <- !is.na(meth_clinical$subtype)

if (sum(sub_ok) > 50 && nlevels(droplevels(meth_clinical$subtype[sub_ok])) > 1) {

  # Does PC1 differ across subtypes? Kruskal-Wallis is the safer default here:
  # subtype group sizes are very unequal, which strains ANOVA's assumptions.
  kw <- kruskal.test(pc1 ~ subtype, data = meth_clinical[sub_ok, ])
  print(kw)
  log_result("B2", "kruskal_p_pc1_by_subtype", kw$p.value)

  print(aggregate(pc1 ~ subtype, data = meth_clinical[sub_ok, ],
                  FUN = function(x) c(n = length(x), mean = mean(x), sd = sd(x))))

  pdf(file.path(paths$diagnostics, "b2_pc1_by_subtype.pdf"), width = 7, height = 5)
  boxplot(pc1 ~ subtype, data = meth_clinical[sub_ok, ],
          xlab = "PAM50 subtype", ylab = "PC1",
          main = "PC1 distribution by subtype", col = "grey85")
  dev.off()

  # Does PC1 still predict survival WITHIN subtype?
  # strata() lets each subtype have its own baseline hazard, so the PC1
  # coefficient is estimated from within-subtype variation only.
  cat("\nPC1 unadjusted for subtype:\n")
  m_plain <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage,
                   data = meth_clinical[sub_ok, ])
  print(summary(m_plain)$coefficients["scale(pc1)", ])

  cat("\nPC1 stratified by subtype:\n")
  m_strat <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage +
                     strata(subtype), data = meth_clinical[sub_ok, ])
  print(summary(m_strat)$coefficients["scale(pc1)", ])

  log_result("B2", "pc1_HR_unstratified",
             exp(coef(m_plain)["scale(pc1)"]))
  log_result("B2", "pc1_HR_subtype_stratified",
             exp(coef(m_strat)["scale(pc1)"]))
  log_result("B2", "pc1_p_subtype_stratified",
             summary(m_strat)$coefficients["scale(pc1)", 5])

  # INTERPRETATION
  # PC1 varies strongly by subtype AND loses significance under strata()
  #   -> PC1 is largely a subtype proxy
  # PC1 survives stratification -> it carries information subtype does not
} else {
  cat("Insufficient subtype coverage -- skipping B2 tests.\n")
}


# =============================================================================
# B3. Tumor purity, handled properly
# =============================================================================
# NOTE: v3 gates purity behind purity_available <- !any(is.na(purity)), which is
# all-or-nothing: 6 missing patients disable purity everywhere. Here we use
# complete-case analysis per model instead, which is what coxph() would do
# anyway, and we test purity against PC1 directly -- a check v3 skips because
# gene_mat does not exist yet when its purity section runs.

cat("\n=== B3. PC1 vs tumor purity ===\n")

purity_path <- "01_processed_data/aran_purity_tcga.xlsx"

if (file.exists(purity_path)) {
  purity_tab <- as.data.frame(read_excel(purity_path, sheet = "Supp Data 1", skip = 2))
  names(purity_tab)[names(purity_tab) == "Sample ID"] <- "sample_barcode"
  purity_tab$CPE <- as.numeric(purity_tab$CPE)
  purity_tab <- purity_tab[purity_tab$`Cancer type` == "BRCA" &
                             !is.na(purity_tab$`Cancer type`), ]
  purity_tab$sample_id_15 <- substr(purity_tab$sample_barcode, 1, 15)

  meth_clinical$purity <- purity_tab$CPE[
    match(substr(meth_clinical$sample_barcode, 1, 15), purity_tab$sample_id_15)
  ]

  pok <- !is.na(meth_clinical$purity)
  cat("Patients with purity:", sum(pok), "of", nrow(meth_clinical), "\n")

  ct_p <- cor.test(meth_clinical$pc1[pok], meth_clinical$purity[pok],
                   method = "spearman")
  print(ct_p)
  log_result("B3", "spearman_rho_pc1_purity", unname(ct_p$estimate),
             sprintf("p = %.3g", ct_p$p.value))

  cat("\nPC1 adjusted for purity:\n")
  m_pur <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage + scale(purity),
                 data = meth_clinical)   # coxph drops incomplete rows itself
  print(summary(m_pur)$coefficients)
  log_result("B3", "pc1_HR_purity_adjusted", exp(coef(m_pur)["scale(pc1)"]))

  pdf(file.path(paths$diagnostics, "b3_pc1_vs_purity.pdf"), width = 6, height = 5)
  plot(meth_clinical$purity[pok], meth_clinical$pc1[pok],
       xlab = "Tumor purity (CPE)", ylab = "PC1",
       main = sprintf("PC1 vs purity (rho = %.2f)", ct_p$estimate),
       pch = 16, col = rgb(0, 0, 0, 0.35))
  lines(lowess(meth_clinical$purity[pok], meth_clinical$pc1[pok]),
        col = "red", lwd = 2)
  dev.off()

  # INTERPRETATION: strong correlation means PC1 partly reflects how much
  # normal tissue is in the sample, not tumor biology. Recall purity alone was
  # NOT prognostic (p = 0.28), so this is about what PC1 IS, not about whether
  # purity predicts outcome.
} else {
  cat("No purity file -- skipping B3.\n")
}


# =============================================================================
# B4. What KIND of methylation change is PC1?
# =============================================================================
# This separates two near-opposite phenomena that get conflated:
#   global hypomethylation -> loss of methylation in open sea / gene deserts,
#     mechanically linked to replication, so it tracks proliferation
#   CIMP -> GAIN of methylation at CpG islands, a different biology entirely
#
# Rather than a probe-level PCA (expensive on 400k x 755), correlate each
# probe's M-value with the PC1 score. Same information, far cheaper, and it
# avoids the promoter-aggregation bias -- promoters are island-enriched by
# construction, so promoter-level loadings alone cannot answer this question.

cat("\n=== B4. Genomic context of PC1 ===\n")

# Fast per-probe correlation with pc1: standardize both, then one crossprod.
pc1_z  <- as.vector(scale(pc1))
m_cent <- m_mat - rowMeans(m_mat)
m_sd   <- sqrt(rowSums(m_cent^2))
probe_cor <- as.vector((m_cent %*% pc1_z) / (m_sd * sqrt(length(pc1_z) - 1)))
names(probe_cor) <- rownames(m_mat)

probe_cor[!is.finite(probe_cor)] <- NA   # constant probes

context <- factor(ann_matched$Relation_to_Island)
cat("\nProbes per genomic context:\n")
print(table(context))

cat("\nCorrelation with PC1 by context (positive = gains methylation as PC1 rises):\n")
ctx_summary <- aggregate(probe_cor ~ context,
                         data = data.frame(probe_cor, context),
                         FUN = function(x) c(n = length(x),
                                             mean = mean(x),
                                             median = median(x)))
print(ctx_summary)

kw_ctx <- kruskal.test(probe_cor ~ context,
                       data = data.frame(probe_cor, context))
print(kw_ctx)
log_result("B4", "kruskal_p_cor_by_context", kw_ctx$p.value)

# Where do the STRONGEST loadings sit?
strong <- !is.na(probe_cor) & abs(probe_cor) > 0.5
cat("\nProbes with |cor| > 0.5:", sum(strong), "\n")
if (sum(strong) > 0) {
  cat("Context breakdown of strong probes:\n")
  print(table(context[strong]))
  cat("Direction of strong probes (negative = loses methylation):\n")
  print(table(context[strong], sign(probe_cor[strong])))

  # Enrichment vs the array background, so "lots of open sea probes" is not
  # mistaken for signal when open sea is simply the biggest category.
  obs <- table(context[strong])
  exp_prop <- table(context) / sum(table(context))
  enrich <- obs / (sum(obs) * as.vector(exp_prop))
  cat("\nEnrichment (observed / expected) among strong probes:\n")
  print(round(enrich, 2))
}

pdf(file.path(paths$diagnostics, "b4_pc1_correlation_by_context.pdf"),
    width = 7, height = 5)
boxplot(probe_cor ~ context, data = data.frame(probe_cor, context),
        xlab = "Relation to CpG island", ylab = "Probe correlation with PC1",
        main = "Where does PC1 act?", col = "grey85", outline = FALSE)
abline(h = 0, lty = 2, col = "red")
dev.off()

# INTERPRETATION
#  Strong NEGATIVE correlations concentrated in OpenSea
#    -> classic global hypomethylation / partially-methylated-domain phenotype,
#       replication-linked, consistent with the B1 proliferation hypothesis
#  Strong POSITIVE correlations concentrated in Island
#    -> CIMP-like hypermethylation, a different biology
#  Both -> PC1 is a composite axis; report it as such rather than naming it

write.csv(data.frame(probe = names(probe_cor), cor_with_pc1 = probe_cor,
                     context = as.character(context)),
          file.path(paths$tables, "b4_probe_correlation_with_pc1.csv"),
          row.names = FALSE)


# =============================================================================
# B5. Data-driven PC1 vs curated CIMP panel (Wang et al. 2017)
# =============================================================================
# Wang H, et al. "Survival differences of CIMP subtypes integrated with CNA
# information in human breast cancer." Oncotarget 2017;8(30):48807-48819.
# doi:10.18632/oncotarget.16178 -- 765 TCGA breast cancers, 450K array,
# 15-marker CpG panel. Closest prior work to your PC1 finding, and the natural
# point of departure: curated panel vs data-driven latent factor, same cohort.
#
# ACTION REQUIRED: the 15 probe IDs are in that paper's supplement. Fill them
# in below. Left empty deliberately -- guessing probe IDs would silently
# produce a meaningless comparison.

cat("\n=== B5. PC1 vs curated CIMP panel ===\n")

cimp_probes <- c()   # <- paste the 15 cg######## IDs from the paper here

if (length(cimp_probes) > 0) {
  present <- intersect(cimp_probes, rownames(m_mat))
  cat("Panel probes found on your array:", length(present), "of", length(cimp_probes), "\n")

  if (length(present) >= 5) {
    cimp_score <- colMeans(m_mat[present, , drop = FALSE])
    meth_clinical$cimp <- cimp_score

    ct_c <- cor.test(meth_clinical$pc1, cimp_score, method = "spearman")
    print(ct_c)
    log_result("B5", "spearman_rho_pc1_cimp", unname(ct_c$estimate))

    # Head-to-head prognostic comparison on identical patients
    m_pc1  <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage,
                    data = meth_clinical)
    m_cimp <- coxph(Surv(time, event) ~ scale(cimp) + age_years + stage,
                    data = meth_clinical)
    m_both <- coxph(Surv(time, event) ~ scale(pc1) + scale(cimp) + age_years + stage,
                    data = meth_clinical)

    cat(sprintf("\nC-index  PC1: %.3f | CIMP panel: %.3f | both: %.3f\n",
                summary(m_pc1)$concordance["C"],
                summary(m_cimp)$concordance["C"],
                summary(m_both)$concordance["C"]))
    cat("\nWhen both are in the model, which survives?\n")
    print(summary(m_both)$coefficients)

    log_result("B5", "c_index_pc1",  summary(m_pc1)$concordance["C"])
    log_result("B5", "c_index_cimp", summary(m_cimp)$concordance["C"])
  }
} else {
  cat("cimp_probes is empty -- fill in the panel from the Wang 2017 supplement.\n")
}


# =============================================================================
# B6. Genomic instability: mutation burden and TP53
# =============================================================================
# Off by default: requires a fresh GDC download (slow, and the first run will
# prompt for disk space). Turn on once B1-B4 have narrowed the hypothesis.

cat("\n=== B6. PC1 vs mutation burden ===\n")

if (run_b6_mutations) {
  maf_cache <- "01_processed_data/brca_maf.rds"

  maf <- cache_or_compute(maf_cache, function() {
    query <- TCGAbiolinks::GDCquery(
      project = "TCGA-BRCA",
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
    )
    TCGAbiolinks::GDCdownload(query)
    TCGAbiolinks::GDCprepare(query)
  })

  maf$patient <- substr(maf$Tumor_Sample_Barcode, 1, 12)

  # Mutation burden: count of non-silent mutations per patient
  nonsilent <- maf[maf$Variant_Classification != "Silent", ]
  tmb_tab <- as.data.frame(table(nonsilent$patient))
  names(tmb_tab) <- c("patient", "n_mut")

  meth_patient <- substr(meth_clinical$sample_barcode, 1, 12)
  meth_clinical$tmb <- tmb_tab$n_mut[match(meth_patient, tmb_tab$patient)]

  tok <- !is.na(meth_clinical$tmb)
  ct_t <- cor.test(meth_clinical$pc1[tok], log1p(meth_clinical$tmb[tok]),
                   method = "spearman")
  print(ct_t)
  log_result("B6", "spearman_rho_pc1_log_tmb", unname(ct_t$estimate))

  # TP53 status: the most common driver, tied to the aggressive phenotype
  tp53_patients <- unique(nonsilent$patient[nonsilent$Hugo_Symbol == "TP53"])
  meth_clinical$tp53 <- factor(ifelse(meth_patient %in% tp53_patients,
                                      "mutant", "wildtype"))
  print(table(meth_clinical$tp53))

  wt <- wilcox.test(pc1 ~ tp53, data = meth_clinical)
  print(wt)
  log_result("B6", "wilcox_p_pc1_by_tp53", wt$p.value)

  pdf(file.path(paths$diagnostics, "b6_pc1_by_tp53.pdf"), width = 6, height = 5)
  boxplot(pc1 ~ tp53, data = meth_clinical, col = "grey85",
          xlab = "TP53 status", ylab = "PC1", main = "PC1 by TP53 status")
  dev.off()
} else {
  cat("run_b6_mutations = FALSE -- skipping (set TRUE once B1-B4 are done).\n")
}


# =============================================================================
# B7. Integration: what is left of PC1?
# =============================================================================
# Throw every candidate explanation in at once and ask what absorbs PC1.
# Two separate questions here, kept separate on purpose:
#   (a) how much of PC1's VARIANCE is explained by known factors? (linear model)
#   (b) does PC1 still predict survival once they are all adjusted for? (Cox)

cat("\n=== B7. Integration ===\n")

# Assemble whichever predictors actually got built above
candidates <- c("prolif", "purity", "subtype", "tmb")
available  <- candidates[candidates %in% names(meth_clinical)]
available  <- available[sapply(available, function(v) sum(!is.na(meth_clinical[[v]])) > 50)]
cat("Predictors available for integration:", paste(available, collapse = ", "), "\n")

if (length(available) > 0) {

  # (a) Variance partition -- how much of PC1 is "known"?
  f_lm <- as.formula(paste("pc1 ~", paste(available, collapse = " + ")))
  lm_pc1 <- lm(f_lm, data = meth_clinical)
  print(summary(lm_pc1))

  r2 <- summary(lm_pc1)$r.squared
  cat(sprintf("\nPC1 is %.1f%% explained by known factors; %.1f%% unexplained.\n",
              100 * r2, 100 * (1 - r2)))
  log_result("B7", "r2_pc1_by_known_factors", r2)

  # Per-predictor contribution. Order matters in sequential ANOVA, so this is
  # descriptive, not a clean decomposition -- report it as such.
  cat("\nSequential ANOVA (order-dependent -- descriptive only):\n")
  print(anova(lm_pc1))

  # Collinearity check: these predictors are likely entangled with each other,
  # which inflates standard errors and makes single coefficients unreliable.
  if (requireNamespace("car", quietly = TRUE) && length(available) > 1) {
    cat("\nVariance inflation factors:\n")
    print(car::vif(lm_pc1))
  }

  # (b) Does PC1 still predict survival after full adjustment?
  f_cox_full <- as.formula(paste("Surv(time, event) ~ scale(pc1) + age_years + stage +",
                                 paste(available, collapse = " + ")))
  f_cox_null <- as.formula(paste("Surv(time, event) ~ age_years + stage +",
                                 paste(available, collapse = " + ")))

  cox_full <- coxph(f_cox_full, data = meth_clinical)
  cox_null <- coxph(f_cox_null, data = meth_clinical)

  cat("\nFull model:\n")
  print(summary(cox_full))

  cat("\nLikelihood ratio test -- does PC1 add anything beyond known factors?\n")
  print(anova(cox_null, cox_full))

  log_result("B7", "pc1_HR_fully_adjusted", exp(coef(cox_full)["scale(pc1)"]))
  log_result("B7", "pc1_p_fully_adjusted",
             summary(cox_full)$coefficients["scale(pc1)", 5])
  log_result("B7", "c_index_full_model", summary(cox_full)$concordance["C"])

  # INTERPRETATION
  # PC1 loses significance -> it was a composite of known confounds
  # PC1 survives           -> there is a real, unexplained epigenetic axis here,
  #                           and THAT is the finding worth writing up
} else {
  cat("No integration predictors available -- run B1-B3 first.\n")
}


# =============================================================================
# 9. Save outputs
# =============================================================================

results_df <- do.call(rbind, results_log)
write.csv(results_df, file.path(paths$tables, "pc1_characterization_results.csv"),
          row.names = FALSE)

saveRDS(meth_clinical, file.path(paths$run, "meth_clinical_with_pc1.rds"))
saveRDS(pca,           file.path(paths$models, "promoter_pca.rds"))

write_manifest(paths$run, list(
  script          = "13_pc1_biological_characterization.R",
  n_samples       = ncol(gene_mat),
  n_promoters     = nrow(gene_mat),
  pc1_var_expl    = round(var_explained, 4),
  b6_mutations    = run_b6_mutations,
  cimp_panel_set  = length(cimp_probes) > 0,
  seed            = 42
))

cat("\n=== Done. Results in", paths$run, "===\n")
print(results_df)
