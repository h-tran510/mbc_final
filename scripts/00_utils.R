# Reusable helper functions for this project
library(survival)

# Saves results into a new timestamped directory

make_run_dir <- function(tag, base) {
  full_path <- file.path(base, paste0(Sys.Date(), "_", tag))
  subdirs <- c("models", "diagnostics", "tables")
  for (subdir in subdirs) {
    dir.create(file.path(full_path, subdir), recursive = TRUE, showWarnings = FALSE)
  }
  list(
    run = full_path,
    models = file.path(full_path, "models"),
    diagnostics = file.path(full_path, "diagnostics"),
    tables = file.path(full_path, "tables")
  )
}

# Writes run_manifest.txt, which includes info about run
# Params includes K, n_repeats, K_inner, prescreen, endpoint, seed

write_manifest <- function(run, params) {
  lines <- paste0(names(params), " = ", unlist(params))
  git_hash <- system("git rev-parse --short HEAD", intern = TRUE)
  session_text <- capture.output(utils::sessionInfo())
  combined_output <- c(
    "=== GIT COMMIT HASH ===", # nolint: indentation_linter.
    git_hash,
    "",
    "=== R SESSION INFO ===",
    session_text,
    "",
    "=== PARAMETERS ===",
    lines)
  writeLines(combined_output, file.path(run, "run_manifest.txt"))
}

# save orr reload objects 
cache_or_compute <- function(path, compute_fn, force = FALSE) {
  if (!force && file.exists(path)) {
    cat("Loading cached:", path, "\n")
    return(readRDS(path))
  }
  cat("Computing (not cached):", path, "\n")
  result <- compute_fn()
  saveRDS(result, path)
  result
}

# check for sample misalignment
run_alignment_checks <- function(m_mat, meth_clinical, id_len = 15) {
  meth_id <- substr(colnames(m_mat), 1, id_len)
  clin_id <- substr(meth_clinical$sample_barcode, 1, id_len)
 
  stopifnot(
    "Duplicate methylation sample IDs after truncation" =
      sum(duplicated(meth_id)) == 0,
    "Duplicate clinical sample IDs after truncation" =
      sum(duplicated(clin_id)) == 0,
    "Methylation and clinical do not share the same sample set" =
      setequal(meth_id, clin_id)
  )
 
  match_idx <- match(clin_id, meth_id)
  stopifnot("Unmatched clinical rows after join" = all(!is.na(match_idx)))
 
  cat("Alignment check passed. Reorder was",
      if (identical(match_idx, seq_along(match_idx))) "a no-op"
      else "necessary and has been applied",
      "\n")
 
  match_idx
}


# univariate cox screen
uni_screen <- function(g, surv_obj, age, stage, purity = NULL) {
  fit <- try(
    if (is.null(purity)) {
      coxph(surv_obj ~ scale(g) + age + stage)
    } else {
      coxph(surv_obj ~ scale(g) + age + stage + scale(purity))
    },
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}
 
# extract probes from the pipeline
extract_probes <- function(fit, s, clin_colnames) {
  coefs        <- as.matrix(coef(fit, s = s))
  selected_all <- rownames(coefs)[coefs[, 1] != 0]
  selected_all[!selected_all %in% clin_colnames]
}