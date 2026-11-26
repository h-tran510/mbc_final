library(survival)

source("deployed_scripts_v1/00_utils.R")

GENES      <- c("C1orf86", "PLA2G3")
ALPHA      <- 0.05
POWER      <- 0.80
N_TN       <- full_or_test(Inf, 40)
GSE_PATH   <- "00_raw_data/external/GSE141441/GSE141441_series_matrix.txt.gz"
ANN_PATH   <- ANN_CACHE_PATH
M_PATH     <- "01_processed_data/meth_m_matrix.rds"
FIG_DIR    <- "03_figures/external_validation"

paths <- make_run_dir(tag = "external_validation", base = "02_analysis")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

tcga_m      <- readRDS(M_PATH)
ann         <- load_matched_annotation(tcga_m, cache_path = ANN_PATH)
tcga_probes <- rownames(tcga_m)
rm(tcga_m)
invisible(gc())

gse   <- load_geo_series_matrix(GSE_PATH)
pheno <- gse$pheno
cat("GSE141441:", nrow(gse$beta), "probes x", ncol(gse$beta), "samples\n")

is_tn  <- grepl(":TN[0-9]+$", pheno$title)
is_ctx <- grepl(":(Case|Ctrl)_[0-9]+$", pheno$title)
stopifnot("A sample title is neither TN nor Case/Ctrl" = all(is_tn | is_ctx))
cat("Excluded from the primary cohort (Case/Ctrl, outcome-selected):",
    sum(!is_tn), "\n")

tn_idx  <- head(which(is_tn), N_TN)
ctx_idx <- which(is_ctx)

stage_raw <- gsub(" ", "", pheno$ajcc)
stage3    <- sub("^(I|II|III)[ABC]?$", "\\1", stage_raw)
stopifnot("An AJCC value does not collapse to I/II/III" =
            all(is.na(stage3) | stage3 %in% c("I", "II", "III")))

clin <- data.frame(
  id    = pheno$geo_accession,
  time  = as.numeric(pheno$dfsmo),
  event = as.integer(pheno$relapse),
  age   = as.numeric(pheno$age),
  stage = factor(stage3, levels = c("I", "II", "III")),
  er    = as.numeric(pheno$er),
  her2  = as.numeric(pheno$her2),
  stringsAsFactors = FALSE
)
stopifnot("Missing outcome in GSE141441" =
            !anyNA(clin$time) && !anyNA(clin$event) && !anyNA(clin$age))

clin$event5 <- ifelse(clin$event == 1 & clin$time <= 60, 1L,
                      ifelse(clin$time >= 60, 0L, NA_integer_))
clin$not_tn <- (!is.na(clin$er) & clin$er != 0) | (!is.na(clin$her2) & clin$her2 != 0)

tn <- clin[tn_idx, ]
cat(sprintf("Primary cohort: n = %d, events = %d, missing stage = %d, ER/HER2-flagged = %d\n",
            nrow(tn), sum(tn$event), sum(is.na(tn$stage)), sum(tn$not_tn)))

to_m <- function(b) {
  b <- pmin(pmax(b, 0.001), 0.999)
  log2(b / (1 - b))
}
beta_all <- gse$beta[, c(tn_idx, ctx_idx), drop = FALSE]
rm(gse)
sample_median <- apply(beta_all, 2, function(b) median(to_m(b), na.rm = TRUE))

present     <- tcga_probes %in% rownames(beta_all)
cat("TCGA screen probes present in GSE141441:", sum(present), "of",
    length(tcga_probes), "\n")

m_tcga   <- to_m(beta_all[tcga_probes[present], , drop = FALSE])
rm(beta_all)
ann_tcga <- ann[present, ]
gene_mat <- build_promoter_gene_mat(m_tcga, ann_tcga)

gene_mat_mc   <- build_promoter_gene_mat(sweep(m_tcga, 2, sample_median), ann_tcga,
                                         verbose = FALSE)
rm(m_tcga)

n_probes <- attr(gene_mat, "n_probes_per_gene")
genes    <- GENES[GENES %in% rownames(gene_mat)]
if (length(genes) < length(GENES)) {
  msg <- paste("No promoter probe in the TCGA probe set for:",
               paste(setdiff(GENES, genes), collapse = ", "))
  if (!TEST_MODE) stop(msg)
  warning(msg, " (test mode: its probe set is the sample's, not the full cohort's)",
          call. = FALSE)
}
stopifnot("None of GENES has a promoter probe" = length(genes) > 0)

rownames(clin) <- clin$id
for (g in genes) {
  clin[[g]]                  <- NA_real_
  clin[[paste0(g, "_mc")]]   <- NA_real_
  clin[colnames(gene_mat), g] <- gene_mat[g, ]
  clin[colnames(gene_mat), paste0(g, "_mc")] <- gene_mat_mc[g, ]
}
stopifnot("Gene matrix columns are not the cohort samples" =
            identical(colnames(gene_mat), clin$id[c(tn_idx, ctx_idx)]))
tn  <- clin[tn_idx, ]
ctx <- clin[ctx_idx, ]

fit_row <- function(gene, model, d, mvar, rhs, logistic = FALSE) {
  if (logistic) {
    lhs <- "event5"
    event_column <- "event5"
    effect_label <- "OR"
    se_column <- "Std. Error"
  } else {
    lhs <- "Surv(time, event)"
    event_column <- "event"
    effect_label <- "HR"
    se_column <- "se(coef)"
  }

  d$M <- d[[mvar]]
  f <- as.formula(paste(lhs, "~", rhs))
  d <- droplevels(d[complete.cases(d[, all.vars(f)]), ])
  ev <- sum(d[[event_column]])
  out <- data.frame(gene = gene, model = model, effect = effect_label,
                     n = nrow(d), events = ev, coef = NA_real_, HR = NA_real_,
                     CI_low = NA_real_, CI_high = NA_real_, p = NA_real_, note = "")
  warns <- character(0)
  fit <- tryCatch(withCallingHandlers(
    if (logistic) glm(f, family = binomial, data = d) else coxph(f, data = d),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) {
    out$note <- conditionMessage(fit)
    return(out)
  }
  out$note <- paste(unique(warns), collapse = "; ")
  if (nzchar(out$note)) warning(gene, " / ", model, ": ", out$note, call. = FALSE)
  cs <- coef(summary(fit))
  if (!"M" %in% rownames(cs) || is.na(cs["M", 1])) {
    out$note <- paste(c("M coefficient not estimable (aliased)", out$note[nzchar(out$note)]),
                      collapse = "; ")
    return(out)
  }
  s  <- cs["M", ]
  se <- s[[se_column]]
  out$coef    <- s[[1]]
  out$HR      <- exp(s[[1]])
  out$CI_low  <- exp(s[[1]] - qnorm(0.975) * se)
  out$CI_high <- exp(s[[1]] + qnorm(0.975) * se)
  out$p       <- s[[length(s)]]
  out
}

adjusted_terms <- "M + age + stage"
results_by_gene <- lapply(genes, function(gene) {
  rbind(
    fit_row(gene, "primary", tn, gene, adjusted_terms),
    fit_row(gene, "median_centred", tn, paste0(gene, "_mc"), adjusted_terms),
    fit_row(gene, "A_logistic_relapse_5y", tn, gene, adjusted_terms, logistic = TRUE),
    fit_row(gene, "B_excl_ER_HER2_flagged", tn[!tn$not_tn, ], gene, adjusted_terms),
    fit_row(gene, "C_unadjusted", tn, gene, "M"),
    fit_row(gene, "context_case_ctrl", ctx, gene, adjusted_terms)
  )
})
results <- do.call(rbind, results_by_gene)
results$p_lt_nominal    <- results$p < ALPHA
results$p_lt_bonferroni <- results$p < ALPHA / length(GENES)

primary <- results[results$model == "primary", ]
verdict_error <- NULL
if (anyNA(primary$p) || any(nzchar(primary$note))) {
  bad <- primary[is.na(primary$p) | nzchar(primary$note), ]
  verdict_error <- paste0("Primary model failed or is degenerate for ",
                          paste(bad$gene, collapse = ", "), ": ",
                          paste(bad$note, collapse = "; "))
}
replicated   <- primary$coef < 0 & primary$p < ALPHA
contradicted <- primary$coef > 0 & primary$p < ALPHA
verdict <- if (!is.null(verdict_error)) {
  NA_character_
} else if (any(contradicted)) {
  "CONTRADICTED"
} else if (all(replicated) && length(genes) == length(GENES)) {
  "REPLICATED"
} else if (any(replicated)) {
  "PARTIAL"
} else {
  "NOT REPLICATED"
}

z_sum <- qnorm(1 - ALPHA / 2) + qnorm(POWER)
detectable <- function(sd, D) exp(-z_sum / (sd * sqrt(D)))
power_tab <- do.call(rbind, lapply(genes, function(g) {
  sd_m <- sd(tn[[g]], na.rm = TRUE)
  D_pr <- primary$events[primary$gene == g]
  data.frame(gene = g, n_missing_M = sum(is.na(tn[[g]])), sd_M = sd_m,
             events_cohort = sum(tn$event),
             events_primary_fit = D_pr,
             HR_per_SD_cohort = detectable(1, sum(tn$event)),
             HR_per_unit_M_cohort = detectable(sd_m, sum(tn$event)),
             HR_per_SD_primary_fit = detectable(1, D_pr),
             HR_per_unit_M_primary_fit = detectable(sd_m, D_pr))
}))

for (g in genes) {
  d    <- tn[!is.na(tn[[g]]), ]
  brks <- quantile(d[[g]], c(0, 1 / 3, 2 / 3, 1))
  if (anyDuplicated(brks)) {
    warning("Tied tertile cut-points for ", g, "; KM plot skipped.", call. = FALSE)
    next
  }
  d$tertile <- cut(d[[g]], brks, include.lowest = TRUE, labels = c("low", "mid", "high"))
  km <- survfit(Surv(time, event) ~ tertile, data = d)
  lr <- survdiff(Surv(time, event) ~ tertile, data = d)
  lr_p <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  with_pdf(file.path(FIG_DIR, paste0("km_tertile_", g, ".pdf")), {
    plot(km, col = c("#1b9e77", "#7570b3", "#d95f02"), lwd = 2, mark.time = TRUE,
         xlab = "Months (dfsmo)", ylab = "Relapse-free",
         main = sprintf("GSE141441 TN, %s promoter M tertiles (n = %d, %d events)\nlog-rank p = %.3g",
                        g, nrow(d), sum(d$event), lr_p))
    legend("bottomleft", legend = paste(levels(d$tertile), "M"), bty = "n",
           col = c("#1b9e77", "#7570b3", "#d95f02"), lwd = 2)
  })
}

git <- MANIFEST_GIT
write.csv(results, file.path(paths$run, "results.csv"), row.names = FALSE)
jsonlite::write_json(list(
  dataset              = "GSE141441 (GPL13534 450K)",
  test_mode            = TEST_MODE,
  genes_requested      = GENES,
  genes_tested         = genes,
  cohort_n             = nrow(tn),
  cohort_events        = sum(tn$event),
  n_excluded_case_ctrl = sum(!is_tn),
  context_n            = nrow(ctx),
  context_events       = sum(ctx$event),
  n_missing_stage      = sum(is.na(tn$stage)),
  n_er_her2_flagged    = sum(tn$not_tn),
  n_censored_before_60m = sum(is.na(tn$event5)),
  tcga_probes          = length(tcga_probes),
  tcga_probes_present  = sum(present),
  probes_per_gene      = as.list(setNames(as.integer(n_probes[genes]), genes)),
  alpha                = ALPHA,
  bonferroni_alpha     = ALPHA / length(GENES),
  power                = POWER,
  detectable_hr_formula = "HR = exp(-(qnorm(1 - alpha/2) + qnorm(power)) / (sd * sqrt(events)))  [Schoenfeld 1983; Hsieh & Lavori 2000]",
  detectable_hr        = power_tab,
  verdict              = verdict,
  verdict_error        = if (is.null(verdict_error)) "" else verdict_error,
  verdict_rule         = "per gene on the primary model: coef<0 & p<alpha replicates, coef>0 & p<alpha contradicts; any contradiction -> CONTRADICTED; all replicated -> REPLICATED; some -> PARTIAL; none -> NOT REPLICATED",
  git_sha              = if (identical(git$sha, "NOT_A_GIT_REPO")) NA_character_ else git$sha,
  git_dirty            = git$dirty,
  run_time             = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
), file.path(paths$run, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
write_manifest(paths$run, list(genes = paste(genes, collapse = ","), verdict = verdict,
                               test_mode = TEST_MODE),
               inputs = c(GSE_PATH, ANN_PATH, M_PATH),
               outputs = c(files_under(paths$run), files_under(FIG_DIR)))

cat(sprintf("\n=== GSE141441 external validation%s ===\n", if (TEST_MODE) " (TEST MODE)" else ""))
cat(sprintf("n = %d, events = %d (missing stage %d; Case/Ctrl excluded %d)\n",
            nrow(tn), sum(tn$event), sum(is.na(tn$stage)), sum(!is_tn)))
cat("Probes per gene:", paste(genes, n_probes[genes], sep = "=", collapse = ", "), "\n\n")
print(results[, c("gene", "model", "effect", "n", "events", "coef", "HR", "CI_low",
                  "CI_high", "p", "note")], digits = 3, row.names = FALSE)
cat("\nDetectable HR at", POWER * 100, "% power, alpha", ALPHA, ":\n")
print(power_tab, digits = 3, row.names = FALSE)
cat("\nBonferroni threshold:", ALPHA / length(GENES), "\n")
cat("VERDICT:", verdict, "\n")
cat("Written:", paths$run, "and", FIG_DIR, "\n")
if (!is.null(verdict_error)) stop(verdict_error, " -- outputs written, no verdict")
