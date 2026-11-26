library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)

source("deployed_scripts_v1/00_utils.R")

paths <- make_run_dir(tag = "pc1_B1_fixed", base = "02_analysis")

run_b6_mutations <- FALSE
set.seed(42)

results_log <- list()

log_result <- function(test, metric, value, note = "", n = NA_integer_) {
  results_log[[length(results_log) + 1]] <<- data.frame(
    test = test, metric = metric, value = unname(value),
    n = as.integer(n), note = note,
    stringsAsFactors = FALSE, row.names = NULL
  )
}

model_n <- function(fit) tryCatch(as.integer(fit$n), error = function(e) NA_integer_)

m_mat         <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

ann_matched <- load_matched_annotation(m_mat)
gene_mat    <- build_promoter_gene_mat(m_mat, ann_matched)

cat("Promoters:", nrow(gene_mat), " Samples:", ncol(gene_mat), "\n")

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

cat("\n=== B1. PC1 vs proliferation ===\n")

expr_path <- "01_processed_data/expr_matrix_vst.rds"
read_rna_se <- FALSE

if (file.exists(expr_path)) {
  expr_matrix <- readRDS(expr_path)

  expr_rownames <- rownames(expr_matrix)

  if (any(grepl("^ENSG", expr_rownames))) {
    cat("Expression rownames look like Ensembl IDs; mapping to symbols.\n")
    brca_se  <- readRDS("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds")
    read_rna_se <- TRUE
    rd       <- as.data.frame(SummarizedExperiment::rowData(brca_se))
    symbol_col <- if ("gene_name" %in% names(rd)) "gene_name" else names(rd)[1]
    ens_strip  <- sub("\\..*$", "", expr_rownames)
    rd_strip   <- sub("\\..*$", "", rownames(rd))
    sym_lookup <- rd[[symbol_col]][match(ens_strip, rd_strip)]
    rm(brca_se); gc()
  } else {
    sym_lookup <- expr_rownames
  }

  prolif_genes <- c("MKI67", "PCNA", "TOP2A", "BIRC5", "CCNB1", "AURKA", "RRM2")
  prolif_idx   <- which(sym_lookup %in% prolif_genes)

  cat("Proliferation genes found:", length(prolif_idx), "of", length(prolif_genes), "\n")

  if (length(prolif_idx) >= 3) {
    assert_one_per_patient(colnames(expr_matrix), "VST expression matrix columns")
    assert_one_per_patient(meth_clinical$sample_barcode, "meth_clinical RNA barcodes")
    expr_col  <- match(meth_clinical$sample_barcode, colnames(expr_matrix))
    stale <- is.na(expr_col) & substr(meth_clinical$sample_barcode, 1, 12) %in%
      substr(colnames(expr_matrix), 1, 12)
    if (any(stale)) {
      stop(sum(stale), " patient(s) have a different RNA aliquot in the VST matrix ",
           "than in meth_clinical (e.g. ", meth_clinical$sample_barcode[which(stale)[1]],
           ") -- rerun from 02.", call. = FALSE)
    }

    cat("Patients with both methylation and expression:",
        sum(!is.na(expr_col)), "of", nrow(meth_clinical), "\n")

    prolif_sub <- expr_matrix[prolif_idx, , drop = FALSE]

    prolif_score <- rowMeans(scale(t(prolif_sub)))
    stopifnot("prolif_score must have one value per sample" =
                length(prolif_score) == ncol(expr_matrix))

    meth_clinical$prolif <- prolif_score[expr_col]

    ok <- !is.na(meth_clinical$prolif)
    cat("Patients contributing to the proliferation correlation:", sum(ok), "\n")

    ct <- cor.test(meth_clinical$pc1[ok], meth_clinical$prolif[ok],
                   method = "spearman")
    print(ct)
    log_result("B1", "spearman_rho_pc1_prolif", unname(ct$estimate),
               sprintf("p = %.3g, n = %d", ct$p.value, sum(ok)))

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

    rm(expr_matrix); gc()
  } else {
    cat("Too few proliferation genes matched -- check the symbol mapping.\n")
  }
} else {
  cat("No expression matrix at", expr_path, "-- skipping B1.\n")
}

cat("\n=== B2. PC1 vs PAM50 subtype ===\n")

sp <- load_subtype_purity(meth_clinical, need_purity = FALSE)
meth_clinical$subtype <- sp$subtype
purity_loaded <- attr(sp, "purity_loaded")
purity_by_row <- sp$purity
rm(sp)

cat("Methylation patients with a PAM50 call:",
    sum(!is.na(meth_clinical$subtype)), "of", nrow(meth_clinical), "\n")

cat("\nSubtype counts (check these BEFORE trusting any per-subtype test):\n")
print(table(meth_clinical$subtype, useNA = "ifany"))

sub_ok <- !is.na(meth_clinical$subtype)

if (sum(sub_ok) > 50 && nlevels(droplevels(meth_clinical$subtype[sub_ok])) > 1) {

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

  cat("\nPC1 unadjusted for subtype:\n")
  m_plain <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage,
                   data = meth_clinical[sub_ok, ])
  print(summary(m_plain)$coefficients["scale(pc1)", ])

  cat("\nPC1 stratified by subtype:\n")
  m_strat <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage +
                     strata(subtype), data = meth_clinical[sub_ok, ])
  print(summary(m_strat)$coefficients["scale(pc1)", ])

  log_result("B2", "pc1_HR_unstratified",
             exp(coef(m_plain)["scale(pc1)"]), n = model_n(m_plain))
  log_result("B2", "pc1_HR_subtype_stratified",
             exp(coef(m_strat)["scale(pc1)"]), n = model_n(m_strat))
  log_result("B2", "pc1_p_subtype_stratified",
             summary(m_strat)$coefficients["scale(pc1)", 5], n = model_n(m_strat))

} else {
  cat("Insufficient subtype coverage -- skipping B2 tests.\n")
}

cat("\n=== B3. PC1 vs tumor purity ===\n")

if (purity_loaded) {
  meth_clinical$purity <- purity_by_row

  pok <- !is.na(meth_clinical$purity)
  cat("Patients with purity:", sum(pok), "of", nrow(meth_clinical), "\n")

  ct_p <- cor.test(meth_clinical$pc1[pok], meth_clinical$purity[pok],
                   method = "spearman")
  print(ct_p)
  log_result("B3", "spearman_rho_pc1_purity", unname(ct_p$estimate),
             sprintf("p = %.3g", ct_p$p.value), n = sum(pok))

  cat("\nPC1 adjusted for purity:\n")
  m_pur <- coxph(Surv(time, event) ~ scale(pc1) + age_years + stage + scale(purity),
                 data = meth_clinical)
  print(summary(m_pur)$coefficients)
  log_result("B3", "pc1_HR_purity_adjusted", exp(coef(m_pur)["scale(pc1)"]),
             note = "complete-case on purity; cohort differs from B2",
             n = model_n(m_pur))

  pdf(file.path(paths$diagnostics, "b3_pc1_vs_purity.pdf"), width = 6, height = 5)
  plot(meth_clinical$purity[pok], meth_clinical$pc1[pok],
       xlab = "Tumor purity (CPE)", ylab = "PC1",
       main = sprintf("PC1 vs purity (rho = %.2f)", ct_p$estimate),
       pch = 16, col = rgb(0, 0, 0, 0.35))
  lines(lowess(meth_clinical$purity[pok], meth_clinical$pc1[pok]),
        col = "red", lwd = 2)
  dev.off()

} else {
  cat("No purity file -- skipping B3.\n")
}

cat("\n=== B4. Genomic context of PC1 ===\n")

pc1_z  <- as.vector(scale(pc1))
m_cent <- m_mat - rowMeans(m_mat)
m_sd   <- sqrt(rowSums(m_cent^2))
probe_cor <- as.vector((m_cent %*% pc1_z) / (m_sd * sqrt(length(pc1_z) - 1)))
names(probe_cor) <- rownames(m_mat)

probe_cor[!is.finite(probe_cor)] <- NA

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

strong <- !is.na(probe_cor) & abs(probe_cor) > 0.5
cat("\nProbes with |cor| > 0.5:", sum(strong), "\n")
if (sum(strong) > 0) {
  cat("Context breakdown of strong probes:\n")
  print(table(context[strong]))
  cat("Direction of strong probes (negative = loses methylation):\n")
  print(table(context[strong], sign(probe_cor[strong])))

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

write.csv(data.frame(probe = names(probe_cor), cor_with_pc1 = probe_cor,
                     context = as.character(context)),
          file.path(paths$tables, "b4_probe_correlation_with_pc1.csv"),
          row.names = FALSE)

cat("\n=== B5. PC1 vs curated CIMP panel ===\n")

cimp_probes <- c()

if (length(cimp_probes) > 0) {
  present <- intersect(cimp_probes, rownames(m_mat))
  cat("Panel probes found on your array:", length(present), "of", length(cimp_probes), "\n")

  if (length(present) >= 5) {
    cimp_score <- colMeans(m_mat[present, , drop = FALSE])
    meth_clinical$cimp <- cimp_score

    ct_c <- cor.test(meth_clinical$pc1, cimp_score, method = "spearman")
    print(ct_c)
    log_result("B5", "spearman_rho_pc1_cimp", unname(ct_c$estimate))

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

    log_result("B5", "c_index_pc1",  summary(m_pc1)$concordance["C"],
               n = model_n(m_pc1))
    log_result("B5", "c_index_cimp", summary(m_cimp)$concordance["C"],
               n = model_n(m_cimp))
    log_result("B5", "c_index_both", summary(m_both)$concordance["C"],
               n = model_n(m_both))
    if (!identical(model_n(m_pc1), model_n(m_cimp))) {
      warning("B5: PC1 and CIMP models were fit on different numbers of ",
              "patients -- their C-indices are not directly comparable.",
              call. = FALSE)
    }
  }
} else {
  cat("cimp_probes is empty -- fill in the panel from the Wang 2017 supplement.\n")
}

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

cat("\n=== B7. Integration ===\n")

candidates <- c("prolif", "purity", "subtype", "tmb")
available <- candidates[candidates %in% names(meth_clinical)]
enough_data <- sapply(available, function(predictor) {
  sum(!is.na(meth_clinical[[predictor]])) > 50
})
available <- available[enough_data]
cat("Predictors available for integration:", paste(available, collapse = ", "), "\n")

if (length(available) > 0) {
  f_lm <- as.formula(paste("pc1 ~", paste(available, collapse = " + ")))
  lm_pc1 <- lm(f_lm, data = meth_clinical)
  print(summary(lm_pc1))

  r2 <- summary(lm_pc1)$r.squared
  cat(sprintf("\nPC1 is %.1f%% explained by known factors; %.1f%% unexplained.\n",
              100 * r2, 100 * (1 - r2)))
  log_result("B7", "r2_pc1_by_known_factors", r2,
             n = tryCatch(as.integer(nobs(lm_pc1)), error = function(e) NA_integer_))

  cat("\nSequential ANOVA (order-dependent -- descriptive only):\n")
  print(anova(lm_pc1))

  if (requireNamespace("car", quietly = TRUE) && length(available) > 1) {
    cat("\nVariance inflation factors:\n")
    print(car::vif(lm_pc1))
  }

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

  log_result("B7", "pc1_HR_fully_adjusted", exp(coef(cox_full)["scale(pc1)"]),
             n = model_n(cox_full))
  log_result("B7", "pc1_p_fully_adjusted",
             summary(cox_full)$coefficients["scale(pc1)", 5], n = model_n(cox_full))
  log_result("B7", "c_index_full_model", summary(cox_full)$concordance["C"],
             n = model_n(cox_full))

} else {
  cat("No integration predictors available -- run B1-B3 first.\n")
}

results_df <- do.call(rbind, results_log)
write.csv(results_df, file.path(paths$tables, "pc1_characterization_results.csv"),
          row.names = FALSE)

saveRDS(meth_clinical, file.path(paths$run, "meth_clinical_with_pc1.rds"))
saveRDS(pca,           file.path(paths$models, "promoter_pca.rds"))

write_manifest(paths$run, list(
  script          = "09_pc1_biological_characterization.R",
  n_samples       = ncol(gene_mat),
  n_promoters     = nrow(gene_mat),
  pc1_var_expl    = round(var_explained, 4),
  b6_mutations    = run_b6_mutations,
  cimp_panel_set  = length(cimp_probes) > 0,
  seed            = 42
), inputs = c(METH_SCREEN_INPUTS, expr_path,
              if (read_rna_se) "00_raw_data/TCGA_BRCA_SummarizedExperiment.rds",
              if (run_b6_mutations) "01_processed_data/brca_maf.rds"))

cat("\n=== Done. Results in", paths$run, "===\n")
print(results_df)
