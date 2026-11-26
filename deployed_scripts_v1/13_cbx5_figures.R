library(survival)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(doParallel)
library(foreach)

source("deployed_scripts_v1/00_utils.R")

fig_dir <- "03_figures/cbx5"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

m_mat         <- readRDS("01_processed_data/meth_m_matrix.rds")
meth_clinical <- readRDS("01_processed_data/meth_clinical.rds")

match_idx <- run_alignment_checks(m_mat, meth_clinical)
m_mat <- m_mat[, match_idx]

surv_obj <- Surv(time = meth_clinical$time, event = meth_clinical$event)

gene_mat <- load_promoter_gene_mat(m_mat)

stopifnot("CBX5 must have a promoter in gene_mat" = "CBX5" %in% rownames(gene_mat))

meth_clinical <- load_subtype_purity(meth_clinical)

complete <- !is.na(meth_clinical$subtype) & !is.na(meth_clinical$purity) &
  !is.na(meth_clinical$age_years) & !is.na(meth_clinical$stage)
meth_c     <- meth_clinical[complete, ]
gene_mat_c <- gene_mat[, complete, drop = FALSE]
surv_c     <- surv_obj[complete]
meth_c$cbx5 <- gene_mat_c["CBX5", ]

cat("Analysis cohort:", sum(complete), "patients,", sum(meth_c$event), "events\n")

cat("\n=== Figure 1: forest plot ===\n")

m_naive <- coxph(surv_c ~ scale(cbx5) + age_years + stage + scale(purity),
                  data = meth_c)

m_strat <- coxph(surv_c ~ scale(cbx5) + age_years + stage + scale(purity) +
                    strata(subtype), data = meth_c)

no_her2      <- meth_c$subtype != "Her2"
meth_no_her2 <- meth_c[no_her2, ]
meth_no_her2$subtype <- droplevels(meth_no_her2$subtype)
surv_no_her2 <- surv_c[no_her2]

m_noher2 <- coxph(surv_no_her2 ~ scale(cbx5) + age_years + stage + scale(purity) +
                     strata(subtype), data = meth_no_her2)

extract_hr <- function(fit, label, n, events) {
  s <- summary(fit)$coefficients["scale(cbx5)", ]
  ci <- summary(fit)$conf.int["scale(cbx5)", ]
  data.frame(
    model  = label,
    hr     = unname(ci["exp(coef)"]),
    lower  = unname(ci["lower .95"]),
    upper  = unname(ci["upper .95"]),
    p      = unname(s["Pr(>|z|)"]),
    n      = n,
    events = events,
    stringsAsFactors = FALSE
  )
}

forest_df <- rbind(
  extract_hr(m_naive,  "Naive (no subtype term)",        nrow(meth_c),        sum(meth_c$event)),
  extract_hr(m_strat,  "Subtype-stratified",              nrow(meth_c),        sum(meth_c$event)),
  extract_hr(m_noher2, "Subtype-stratified, Her2 excluded", nrow(meth_no_her2), sum(meth_no_her2$event))
)
forest_df <- forest_df[nrow(forest_df):1, ]

print(forest_df)
write.csv(forest_df, file.path(fig_dir, "cbx5_forest_data.csv"), row.names = FALSE)

with_pdf(file.path(fig_dir, "cbx5_forest_plot.pdf"), {
  par(mar = c(4.5, 15, 2, 2))
  xlim <- c(min(forest_df$lower) * 0.85, max(forest_df$upper) * 1.15)
  plot(NA, xlim = log(xlim), ylim = c(0.5, nrow(forest_df) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "Hazard ratio (log scale)", ylab = "",
       main = "CBX5 promoter methylation: hazard ratio across model specifications")
  axis_ticks <- pretty(xlim, n = 6)
  axis_ticks <- axis_ticks[axis_ticks > 0]
  axis(1, at = log(axis_ticks), labels = axis_ticks)
  axis(2, at = seq_len(nrow(forest_df)),
       labels = sprintf("%s\n(n=%d, events=%d)", forest_df$model, forest_df$n, forest_df$events),
       las = 1, cex.axis = 0.8)
  abline(v = log(1), lty = 2, col = "grey50")
  segments(log(forest_df$lower), seq_len(nrow(forest_df)),
           log(forest_df$upper), seq_len(nrow(forest_df)), lwd = 2, col = "steelblue")
  points(log(forest_df$hr), seq_len(nrow(forest_df)), pch = 15, cex = 1.4, col = "steelblue")
  text(log(forest_df$upper), seq_len(nrow(forest_df)),
       labels = sprintf("  HR=%.2f, p=%.2g", forest_df$hr, forest_df$p),
       pos = 4, cex = 0.75, xpd = NA)
}, width = 7.5, height = 4)

cat("Saved:", file.path(fig_dir, "cbx5_forest_plot.pdf"), "\n")

cat("\n=== Figure 2: Kaplan-Meier by subtype ===\n")

meth_c$cbx5_group <- factor(
  ifelse(meth_c$cbx5 >= median(meth_c$cbx5), "High methylation", "Low methylation"),
  levels = c("Low methylation", "High methylation")
)

subtypes_present <- levels(droplevels(meth_c$subtype))
n_panels <- length(subtypes_present)
n_col <- 2
n_row <- ceiling(n_panels / n_col)

grp_cols <- c("Low methylation" = "steelblue", "High methylation" = "tomato")
km_summary <- list()

with_pdf(file.path(fig_dir, "cbx5_km_by_subtype.pdf"), {
  par(mfrow = c(n_row, n_col))

  for (st in subtypes_present) {
    sub_df <- meth_c[meth_c$subtype == st, ]

    sub_df$cbx5_group <- droplevels(sub_df$cbx5_group)

    if (nlevels(sub_df$cbx5_group) < 2) {
      plot.new()
      title(main = sprintf("%s (n=%d): all patients on one side of the\nglobal median -- no comparison shown",
                           st, nrow(sub_df)), cex.main = 0.9)
      km_summary[[st]] <- data.frame(subtype = st, n = nrow(sub_df),
                                      events = sum(sub_df$event), logrank_p = NA_real_)
      next
    }

    sub_surv  <- Surv(time = sub_df$time, event = sub_df$event)
    cols_here <- unname(grp_cols[levels(sub_df$cbx5_group)])

    km_fit <- survfit(sub_surv ~ cbx5_group, data = sub_df)
    lr <- tryCatch(survdiff(sub_surv ~ cbx5_group, data = sub_df), error = function(e) NULL)
    lr_p <- if (!is.null(lr)) {
      pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
    } else {
      NA_real_
    }

    plot(km_fit, col = cols_here, lty = 1, lwd = 2,
         xlab = "Time (days)", ylab = "Progression-free probability",
         main = sprintf("%s (n=%d, events=%d)", st, nrow(sub_df), sum(sub_df$event)),
         mark.time = TRUE)
    legend("bottomleft",
           legend = sprintf("%s (n=%d)", levels(sub_df$cbx5_group),
                            as.integer(table(sub_df$cbx5_group))),
           col = cols_here, lty = 1, lwd = 2, bty = "n", cex = 0.8)
    if (!is.na(lr_p)) {
      mtext(paste0("Log-rank p = ", format.pval(lr_p, digits = 2)), side = 3, line = 0.2, cex = 0.8)
    }

    km_summary[[st]] <- data.frame(subtype = st, n = nrow(sub_df),
                                    events = sum(sub_df$event), logrank_p = lr_p)
  }
}, width = 5 * n_col, height = 4.5 * n_row)

km_summary_df <- do.call(rbind, km_summary)
print(km_summary_df)
write.csv(km_summary_df, file.path(fig_dir, "cbx5_km_by_subtype_summary.csv"), row.names = FALSE)

cat("Saved:", file.path(fig_dir, "cbx5_km_by_subtype.pdf"), "\n")

cat("\n=== Figure 3: permutation null histogram ===\n")

list_validation_runs <- function() {
  if (!dir.exists("02_analysis")) return(character(0))
  run_dirs <- list.dirs("02_analysis", recursive = FALSE)
  run_dirs <- run_dirs[grepl("stratified_screen_validation$", basename(run_dirs))]
  run_dirs[order(basename(run_dirs), decreasing = TRUE)]
}

N_PERM <- full_or_test(100, 5)

this_input_md5  <- screen_input_hash(gene_mat_c, surv_c, meth_c)
this_script_md5 <- script11_code_hash()
validation_run_matches <- function(run_dir) {
  perm_csv      <- file.path(run_dir, "tables", "permutation_null_results.csv")
  manifest_path <- file.path(run_dir, "run_manifest.txt")
  if (!file.exists(perm_csv) || !file.exists(manifest_path)) {
    cat("Not reusing", run_dir, "-- permutation table or manifest missing.\n")
    return(FALSE)
  }
  lines <- readLines(manifest_path)
  get <- function(key) tryCatch(read_manifest_value(lines, key),
                                error = function(e) NA_character_)
  same_input  <- identical(get("screen_input_md5"), this_input_md5)
  same_script <- identical(get("script_md5"), this_script_md5)
  need   <- N_PERM - max_perm_failures(N_PERM)
  n_perm <- if (same_input && same_script) {
    tryCatch(nrow(read.csv(perm_csv)), error = function(e) NA_integer_)
  } else NA_integer_
  enough <- !is.na(n_perm) && n_perm >= need
  null_ok <- identical(get("permutation_null"), "freedman_lane_within_strata")
  ok <- same_input && same_script && enough && null_ok
  if (!ok) {
    cat(sprintf(paste0("Not reusing %s -- provenance mismatch: screen input %s, script 11 %s,",
                       " completed permutations %s (need >= %d), permutation_null %s (need freedman_lane_within_strata).",
                       " A missing key means an old manifest.\n"),
                run_dir, if (same_input) "match" else "differ/missing",
                if (same_script) "match" else "differ/missing",
                if (is.na(n_perm)) "not read" else n_perm, need,
                get("permutation_null")))
  }
  ok
}

latest_run <- NULL
for (run_dir in list_validation_runs()) {
  if (validation_run_matches(run_dir)) {
    latest_run <- run_dir
    break
  }
}
reuse_existing <- !is.null(latest_run)
if (reuse_existing) {
  cat("Reusing existing permutation results from:", latest_run, "\n")
  perm_results   <- read.csv(file.path(latest_run, "tables", "permutation_null_results.csv"))
  manifest_lines <- readLines(file.path(latest_run, "run_manifest.txt"))

  obs_n_hits      <- as.integer(read_manifest_value(manifest_lines, "observed_hits_fdr10"))
  perm_p_hitcount <- as.numeric(read_manifest_value(manifest_lines, "permutation_p_hitcount"))
  stopifnot(
    "observed_hits_fdr10 did not parse to a finite number" = is.finite(obs_n_hits),
    "permutation_p_hitcount did not parse to a finite number" = is.finite(perm_p_hitcount)
  )
  get11 <- function(key) tryCatch(read_manifest_value(manifest_lines, key),
                                  error = function(e) NA_character_)
  perm_fail <- list(n_failed = get11("n_failed"), n_failed_fits = get11("n_failed_fits"),
                    median_fit_fail_frac = get11("median_fit_fail_frac"))
} else {
  cat("No reusable validation run found -- recomputing (this takes a few minutes)...\n")

  obs_res    <- run_full_screen(gene_mat_c, surv_c, meth_c)
  obs_padj   <- p.adjust(obs_res[, "p"], method = "BH")
  obs_n_hits <- sum(obs_padj < 0.10, na.rm = TRUE)

  n_cores <- n_workers()

  fl <- freedman_lane_prepare(gene_mat_c, meth_c)
  perm_out <- run_fl_permutations(
    fl, surv_c, meth_c, N_PERM, n_cores, run_full_screen, stratified_screen_one,
    freedman_lane_perm,
    file.path("02_analysis", "logs", "13_cbx5_figures"))
  perm_results <- perm_out$results
  perm_fail    <- perm_out[c("n_failed", "n_failed_fits", "median_fit_fail_frac")]
  perm_p_hitcount <- (sum(perm_results$n_hits >= obs_n_hits) + 1) / (nrow(perm_results) + 1)
}

with_pdf(file.path(fig_dir, "cbx5_permutation_null.pdf"), {
  hist(perm_results$n_hits, breaks = max(10, nrow(perm_results) / 5), col = "grey80",
       border = "white", xlab = "Promoters at BH q<0.10 (null permutations)",
       main = sprintf("Permutation null (n=%d) vs observed (perm p=%.3f)",
                       nrow(perm_results), perm_p_hitcount),
       xlim = range(c(perm_results$n_hits, obs_n_hits)))
  abline(v = obs_n_hits, col = "tomato", lwd = 2, lty = 2)
  legend("topright", legend = sprintf("Observed = %d", obs_n_hits),
         col = "tomato", lty = 2, lwd = 2, bty = "n")
}, width = 7, height = 5)

cat(sprintf("Observed hits: %d | permutation p (hit count): %.3f\n", obs_n_hits, perm_p_hitcount))
cat("Saved:", file.path(fig_dir, "cbx5_permutation_null.pdf"), "\n")

write_manifest(fig_dir, list(
  script               = "13_cbx5_figures.R",
  figure3_source       = if (reuse_existing) latest_run else "recomputed in this run",
  n_permutations       = N_PERM,
  n_permutations_used  = nrow(perm_results),
  n_failed             = perm_fail$n_failed,
  n_failed_fits        = perm_fail$n_failed_fits,
  median_fit_fail_frac = perm_fail$median_fit_fail_frac,
  observed_hits_fdr10  = obs_n_hits,
  permutation_p_hitcount = round(perm_p_hitcount, 4)
), inputs = c(METH_SCREEN_INPUTS,
              if (reuse_existing) file.path(latest_run, c("tables/permutation_null_results.csv",
                                                          "run_manifest.txt"))),
   outputs = c(files_under(fig_dir),
               files_under(file.path("02_analysis", "logs", "13_cbx5_figures"))))
cat("NOTE: this is the whole-screen permutation null (validates the screen's",
    "BH q<0.10 candidate list collectively), not a CBX5-specific permutation test --",
    "consistent with how this test was originally run in 11_stratified_screen_validation.R.\n")

cat("\n=== Done. All CBX5 figures saved to", fig_dir, "===\n")
