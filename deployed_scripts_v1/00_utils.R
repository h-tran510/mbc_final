library(survival)

TEST_MODE <- as.integer(Sys.getenv("MBC_SUBSET_N", "0")) > 0
full_or_test <- function(full, test) if (TEST_MODE) test else full

make_run_dir <- function(tag, base, allow_existing = FALSE) {
  stamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  full_path <- file.path(base, paste0(stamp, "_", tag))

  if (dir.exists(full_path) && !allow_existing) {
    stop("Run directory already exists, refusing to overwrite: ", full_path,
         "\nPass allow_existing = TRUE if writing into it is intended.")
  }

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

git_out <- function(args) {
  tryCatch(system2("git", shQuote(c("--no-optional-locks", args)), stdout = TRUE, stderr = FALSE,
                   timeout = 60),
           error = function(e) NULL, warning = function(w) NULL)
}

git_provenance <- function() {
  sha <- git_out(c("rev-parse", "HEAD"))
  if (length(sha) == 0 || !nzchar(sha[1])) {
    return(list(sha = "NOT_A_GIT_REPO", dirty = NA))
  }
  status <- git_out(c("status", "--porcelain", "--untracked-files=all", "--",
                      ":/deployed_scripts_v1/*.R", ":/test_mode/*.R", ":/test_mode/*.sh"))
  list(sha = sha[1], dirty = if (is.null(status)) NA else length(status) > 0)
}

MANIFEST_START_TIME <- Sys.time()
MANIFEST_GIT <- git_provenance()

files_under <- function(dir, since = MANIFEST_START_TIME) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  written_since(f, since)
}

written_since <- function(paths, since = MANIFEST_START_TIME) {
  m <- as.numeric(file.mtime(paths))
  s <- as.numeric(since)
  paths[!is.na(m) & m >= ifelse(m == floor(m), floor(s), s)]
}

file_hash_lines <- function(paths, mark_written = FALSE) {
  paths <- unique(as.character(paths))
  if (length(paths) == 0) return("(none)")
  vapply(paths, function(p) {
    if (!file.exists(p)) return(paste(p, "absent", "absent", sep = "\t"))
    if (dir.exists(p)) return(paste(p, "directory", "directory", sep = "\t"))
    md5 <- tryCatch(unname(tools::md5sum(p)), error = function(e) NA_character_)
    line <- if (is.na(md5)) paste(p, "unreadable", "unreadable", sep = "\t") else
      paste(p, sprintf("%.0f", file.size(p)), md5, sep = "\t")
    if (mark_written && length(written_since(p)) > 0) {
      line <- paste(line, "written_during_run", sep = "\t")
    }
    line
  }, character(1), USE.NAMES = FALSE)
}

manifest_script_name <- function(script = NULL) {
  if (!is.null(script)) return(script)
  for (fr in rev(sys.frames())) {
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      f <- get("ofile", envir = fr, inherits = FALSE)
      if (is.character(f) && length(f) == 1) return(basename(f))
    }
  }
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) basename(f[1]) else "unknown (not run by Rscript <file> or source())"
}

ANN_CACHE_PATH <- "01_processed_data/ann_matched_cache.rds"
SUBTYPE_PURITY_INPUTS <- c(pam50  = "01_processed_data/PAM50_classification.txt",
                           purity = "01_processed_data/aran_purity_tcga.xlsx")
METH_SCREEN_INPUTS <- unname(c("01_processed_data/meth_m_matrix.rds",
                               "01_processed_data/meth_clinical.rds",
                               ANN_CACHE_PATH, SUBTYPE_PURITY_INPUTS))

write_manifest <- function(run, params, inputs = character(0), outputs = NULL,
                           script = NULL, file = "run_manifest.txt") {
  lines <- paste0(names(params), " = ", unlist(params))
  manifest_path <- file.path(run, file)
  if (is.null(outputs)) outputs <- files_under(run)
  outputs <- setdiff(outputs, manifest_path)

  git <- MANIFEST_GIT
  if (identical(git$sha, "NOT_A_GIT_REPO")) {
    message("WARNING: could not read a git commit hash; manifest provenance ",
            "will read NOT_A_GIT_REPO. Results are unaffected.")
  }

  session_text <- capture.output(utils::sessionInfo())
  stamp <- function(t) format(t, "%Y-%m-%d %H:%M:%S %Z")
  combined_output <- c(
    "=== GIT COMMIT HASH ===",
    git$sha,
    "",
    "=== PROVENANCE ===",
    paste0("git_dirty: ", git$dirty),
    paste0("script: ", manifest_script_name(script)),
    paste0("start_time: ", stamp(MANIFEST_START_TIME)),
    paste0("end_time: ", stamp(Sys.time())),
    "",
    "=== INPUT FILES (path, bytes, md5 at end of run) ===",
    file_hash_lines(inputs, mark_written = TRUE),
    "",
    "=== OUTPUT FILES (path, bytes, md5 at end of run) ===",
    file_hash_lines(outputs),
    "",
    "=== R SESSION INFO ===",
    session_text,
    "",
    "=== PARAMETERS ===",
    lines)
  writeLines(combined_output, manifest_path)
}

read_manifest_value <- function(lines, key) {
  hit <- grep(paste0("^", key, " = "), lines, value = TRUE)
  if (length(hit) != 1L) {
    stop(sprintf("Expected exactly one '%s' line in the manifest, found %d.",
                 key, length(hit)))
  }
  sub(paste0("^", key, " = "), "", hit)
}

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

env_int <- function(name, default) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) return(default)
  n <- suppressWarnings(as.integer(v))
  if (is.na(n) || n < 1 || !grepl("^[0-9]+$", trimws(v))) {
    stop(name, " must be a positive integer, got '", v, "'", call. = FALSE)
  }
  n
}

n_workers <- function(test_cap = 5L) {
  cores <- parallel::detectCores()
  default <- if (is.na(cores)) {
    warning("parallel::detectCores() returned NA -- using 2 workers (set MBC_N_CORES to override).",
            call. = FALSE)
    2L
  } else max(1L, as.integer(cores) - 1L)
  n <- env_int("MBC_N_CORES", default)
  if (TEST_MODE) min(n, test_cap) else n
}

make_worker_cluster <- function(n, log_dir) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(log_dir, "workers.log")
  file.create(log_path)
  cat("Starting", n, "worker(s); worker output ->", log_path, "\n")
  parallel::makeCluster(n, outfile = log_path)
}

max_perm_failures <- function(n_total, max_frac = 0.05) {
  max(1L, as.integer(floor(max_frac * n_total)))
}

check_perm_failures <- function(failed_ids, n_total, max_frac = 0.05) {
  n_failed <- length(failed_ids)
  if (n_failed > 0) {
    cat(sprintf("%d of %d permutations failed and were dropped (permutations %s; see workers.log).\n",
                n_failed, n_total, paste(failed_ids, collapse = ", ")))
  }
  limit <- max_perm_failures(n_total, max_frac)
  if (n_failed > limit) {
    stop(sprintf("%d of %d permutations failed (more than %d allowed) -- the permutation null is not trusted. See workers.log.",
                 n_failed, n_total, limit), call. = FALSE)
  }
  invisible(n_failed)
}

stratified_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity) + strata(meth_x$subtype)),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

run_full_screen <- function(gm, surv_x, meth_x, one = stratified_screen_one) {
  res <- t(apply(gm, 1, one, surv_x = surv_x, meth_x = meth_x))
  colnames(res) <- c("coef", "z", "p")
  res
}

screen_candidates <- function(p, q_max = 0.10) {
  q    <- p.adjust(p, method = "BH")
  hits <- names(q)[!is.na(q) & q < q_max]
  hits[order(q[hits])]
}

freedman_lane_prepare <- function(gm, meth_x) {
  meth_x$subtype <- droplevels(factor(meth_x$subtype))
  meth_x$stage   <- droplevels(factor(meth_x$stage))
  Z <- model.matrix(~ age_years + stage + scale(purity) + subtype, data = meth_x)
  stopifnot("Freedman-Lane design lost rows (NA covariates?)" = nrow(Z) == ncol(gm))
  qrZ <- qr(Z)
  Y   <- t(gm)
  list(fitted = qr.fitted(qrZ, Y), resid = qr.resid(qrZ, Y), dimnames = dimnames(gm),
       strata = split(seq_len(ncol(gm)), meth_x$subtype))
}

strata_perm <- function(perm_idx, strata) {
  out <- seq_along(perm_idx)
  for (w in strata) out[w] <- w[rank(perm_idx[w])]
  out
}

freedman_lane_perm <- function(fl, perm_idx) {
  perm_w <- strata_perm(perm_idx, fl$strata)
  gm_fl <- t(fl$fitted + fl$resid[perm_w, , drop = FALSE])
  dimnames(gm_fl) <- fl$dimnames
  gm_fl
}

run_fl_permutations <- function(fl, surv_c, meth_c, N_PERM, n_cores,
                                run_full_screen, stratified_screen_one,
                                freedman_lane_perm, worker_log_dir,
                                fit_fail_max = PERM_FIT_FAIL_MAX) {
  force(fit_fail_max)
  cl <- make_worker_cluster(n_cores, worker_log_dir)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  doParallel::registerDoParallel(cl)

  n_pat <- nrow(fl$resid)
  `%dopar%` <- foreach::`%dopar%`
  res <- foreach::foreach(
    r = seq_len(N_PERM),
    .packages = c("survival"),
    .export   = c("run_full_screen", "stratified_screen_one", "freedman_lane_perm",
                  "strata_perm"),
    .combine  = rbind
  ) %dopar% {
    tryCatch({
      set.seed(1000 + r)
      perm_idx <- sample(n_pat)
      res_p    <- run_full_screen(freedman_lane_perm(fl, perm_idx), surv_c, meth_c)
      n_fit_failed <- sum(is.na(res_p[, "p"]))
      if (n_fit_failed > fit_fail_max * nrow(res_p)) {
        message(sprintf("perm %d seed %d failed: %d of %d Cox fits returned NA",
                        r, 1000 + r, n_fit_failed, nrow(res_p)))
        c(permutation = r, n_hits = NA, min_p = NA,
          n_fit_failed = n_fit_failed, n_fits = nrow(res_p))
      } else {
        padj_p <- p.adjust(res_p[, "p"], method = "BH")
        c(permutation = r,
          n_hits = sum(padj_p < 0.10, na.rm = TRUE),
          min_p  = min(res_p[, "p"], na.rm = TRUE),
          n_fit_failed = n_fit_failed, n_fits = nrow(res_p))
      }
    }, error = function(e) {
      message(sprintf("perm %d seed %d failed: %s", r, 1000 + r, conditionMessage(e)))
      c(permutation = r, n_hits = NA, min_p = NA, n_fit_failed = NA, n_fits = NA)
    })
  }

  res <- as.data.frame(res)
  failed <- is.na(res$n_hits)
  fits   <- perm_fit_failures(res$n_fit_failed, res$n_fits, failed)
  n_failed <- check_perm_failures(res$permutation[failed], N_PERM)
  c(list(results = res[!failed, ], n_failed = n_failed), fits)
}

PERM_FIT_FAIL_MAX <- 0.5

perm_fit_failures <- function(n_fit_failed, n_fits, failed) {
  frac <- n_fit_failed / n_fits
  out <- list(n_failed_fits        = sum(failed & !is.na(frac)),
              median_fit_fail_frac = if (all(is.na(frac))) NA_real_
                                     else round(stats::median(frac, na.rm = TRUE), 4))
  cat(sprintf("Cox fit failures per permutation: median %.2f%% of promoters; %d permutation(s) failed on > %.0f%% NA fits (R30).\n",
              100 * out$median_fit_fail_frac, out$n_failed_fits, 100 * PERM_FIT_FAIL_MAX))
  out
}

gdc_cached <- function(path, fetch, max_age_days = 7,
                       tries = env_int("MBC_GDC_TRIES", 10L),
                       wait_s = env_int("MBC_GDC_WAIT_S", 60L)) {
  age_days <- if (file.exists(path)) {
    as.numeric(difftime(Sys.time(), file.mtime(path), units = "days"))
  } else Inf
  if (age_days <= max_age_days) {
    cat(sprintf("Reading local copy (%.1f days old): %s\n", age_days, path))
    return(readRDS(path))
  }
  cat(if (is.finite(age_days)) sprintf("Local copy is %.1f days old -- refreshing from GDC: %s\n", age_days, path)
      else sprintf("No local copy -- downloading from GDC: %s\n", path))

  for (i in seq_len(tries)) {
    cat(sprintf("[%s] GDC attempt %d of %d\n", format(Sys.time(), "%H:%M:%S"), i, tries))
    value <- tryCatch(fetch(), error = function(e) e)
    if (!inherits(value, "error")) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(value, path)
      return(value)
    }
    if (!grepl("resolve hostname|server down|timed out|Timeout|Could not|curl|SSL|connection",
               conditionMessage(value), ignore.case = TRUE)) {
      stop("GDC fetch failed with a non-network error (not retried): ",
           conditionMessage(value), call. = FALSE)
    }
    cat(sprintf("GDC attempt %d of %d failed: %s\n", i, tries, conditionMessage(value)))
    if (i < tries) Sys.sleep(wait_s)
  }
  if (file.exists(path)) {
    warning(sprintf("GDC refresh failed after %d attempts; using the %.1f-day-old local copy: %s",
                    tries, age_days, path), call. = FALSE)
    return(readRDS(path))
  }
  stop("GDC download failed after ", tries, " attempts and there is no local copy: ", path)
}

cache_or_compute_keyed <- function(path, key, compute_fn, force = FALSE) {
  if (!force && file.exists(path)) {
    cached <- try(readRDS(path), silent = TRUE)
    if (!inherits(cached, "try-error") && is.list(cached) &&
          identical(names(cached), c("key", "value")) &&
          identical(cached$key, key)) {
      cat("Loading cached (key matched):", path, "\n")
      return(cached$value)
    }
    cat("Cache at", path, "was built from different input -- recomputing.\n")
  }
  cat("Computing (not cached):", path, "\n")
  value <- compute_fn()
  saveRDS(list(key = key, value = value), path)
  value
}

promoter_gene_pairs <- function(ref_name, ref_group,
                                tss_groups = c("TSS200", "TSS1500")) {
  nm <- strsplit(ifelse(is.na(ref_name), "", ref_name), ";", fixed = TRUE)
  gr <- strsplit(ifelse(is.na(ref_group), "", ref_group), ";", fixed = TRUE)
  stopifnot("UCSC_RefGene_Name and UCSC_RefGene_Group are not positionally paired" =
              identical(lengths(nm), lengths(gr)))
  pairs <- data.frame(probe = rep(seq_along(nm), lengths(nm)),
                      gene  = unlist(nm, use.names = FALSE),
                      group = unlist(gr, use.names = FALSE),
                      stringsAsFactors = FALSE)
  pairs <- pairs[pairs$group %in% tss_groups & nzchar(pairs$gene), c("probe", "gene")]
  pairs <- pairs[!duplicated(pairs), ]
  rownames(pairs) <- NULL
  pairs
}

build_promoter_gene_mat <- function(m_mat, ann_matched, verbose = TRUE) {
  stopifnot(
    "Annotation is not row-aligned to m_mat (row counts differ)" =
      nrow(ann_matched) == nrow(m_mat),
    "Annotation is not row-aligned to m_mat (probe IDs differ)" =
      isTRUE(all.equal(rownames(ann_matched), rownames(m_mat)))
  )

  pairs <- promoter_gene_pairs(ann_matched$UCSC_RefGene_Name,
                               ann_matched$UCSC_RefGene_Group)
  if (verbose) {
    cat("Promoter probes:", length(unique(pairs$probe)),
        "| probe-gene pairs:", nrow(pairs),
        "| probes that are a TSS probe of more than one gene:",
        sum(tabulate(pairs$probe) > 1), "\n")
  }

  sums   <- rowsum(m_mat[pairs$probe, , drop = FALSE], group = pairs$gene)
  counts <- table(pairs$gene)

  stopifnot(
    "rowsum/table group ordering diverged -- promoter means would be misaligned" =
      identical(rownames(sums), names(counts))
  )

  gene_mat <- sums / as.numeric(counts)
  attr(gene_mat, "n_probes_per_gene") <- counts
  gene_mat
}

load_matched_annotation <- function(m_mat,
                                    cache_path = ANN_CACHE_PATH) {
  ann <- cache_or_compute_keyed(
    path = cache_path,
    key  = rownames(m_mat),
    compute_fn = function() {
      ann_df <- minfi::getAnnotation(
        IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19
      )
      ann_df[match(rownames(m_mat), rownames(ann_df)), ]
    }
  )
  stopifnot("annotation rows are not the m_mat probes in order" =
              identical(as.character(ann$Name), rownames(m_mat)))
  ann
}

load_promoter_gene_mat <- function(m_mat,
                                   cache_path = ANN_CACHE_PATH,
                                   verbose = TRUE) {
  ann_matched <- load_matched_annotation(m_mat, cache_path = cache_path)
  build_promoter_gene_mat(m_mat, ann_matched, verbose = verbose)
}

assert_annotation_cache <- function(m_mat,
                                    cache_path = ANN_CACHE_PATH) {
  cached <- readRDS(cache_path)
  stopifnot(
    "annotation cache is not in the keyed list(key, value) format" =
      is.list(cached) && identical(names(cached), c("key", "value")),
    "annotation cache key differs from rownames(m_mat) -- stale cache" =
      identical(cached$key, rownames(m_mat))
  )
  cat("Annotation cache key matches m_mat:", nrow(m_mat), "probes\n")
  invisible(TRUE)
}

with_pdf <- function(path, expr, width = 7, height = 6) {
  pdf(path, width = width, height = height)
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  force(expr)
  invisible(path)
}

row_vars <- function(m) {
  v <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowVars(as.matrix(m))
  } else {
    apply(m, 1, var)
  }
  names(v) <- rownames(m)
  v
}

screen_input_hash <- function(gene_mat_c, surv_c, meth_c) {
  x <- list(gene_mat_c, unclass(surv_c),
            as.data.frame(meth_c)[, c("sample_barcode", "age_years", "stage",
                                      "purity", "subtype")])
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb")
  serialize(x, con, version = 2)
  close(con)
  unname(tools::md5sum(f))
}

script11_code_hash <- function() {
  paste(unname(tools::md5sum(c("deployed_scripts_v1/11_stratified_screen_validation.R",
                               "deployed_scripts_v1/00_utils.R"))), collapse = "+")
}

RNA_ALIQUOT_RULE <- paste(
  "one RNA aliquot per patient: frozen before FFPE (colData preservation_method),",
  "then the lowest vial letter (01A before 01B), then the highest plate id",
  "(barcode field 6, the most recent sequencing plate), then the highest full barcode")

pick_one_aliquot <- function(barcodes, ffpe) {
  stopifnot("duplicate barcodes" = !anyDuplicated(barcodes),
            "ffpe must be TRUE/FALSE per barcode" =
              length(ffpe) == length(barcodes) && !anyNA(ffpe))
  f <- strsplit(barcodes, "-", fixed = TRUE)
  stopifnot("every barcode needs the 7 TCGA fields" = all(lengths(f) == 7))
  patient <- substr(barcodes, 1, 12)
  type <- substr(barcodes, 14, 15)
  mixed <- unique(patient[duplicated(patient) &
                            !duplicated(paste(patient, type))])
  if (length(mixed) > 0) {
    stop(length(mixed), " patient(s) have RNA aliquots of more than one sample ",
         "type (", paste(head(mixed, 10), collapse = ", "), ")", call. = FALSE)
  }
  vial    <- substr(vapply(f, `[`, "", 4), 3, 3)
  plate   <- vapply(f, `[`, "", 6)
  o <- order(patient, ffpe, vial, plate, barcodes, method = "radix",
             decreasing = c(FALSE, FALSE, FALSE, TRUE, TRUE))
  keep <- logical(length(barcodes))
  keep[o[!duplicated(patient[o])]] <- TRUE
  multi  <- patient %in% patient[duplicated(patient)]
  choice <- data.frame(patient_id = patient, barcode = barcodes, ffpe = ffpe,
                       chosen = keep)[multi, ]
  choice <- choice[order(choice$patient_id, !choice$chosen, choice$barcode,
                         method = "radix"), ]
  rownames(choice) <- NULL
  list(keep = keep, choice = choice)
}

assert_one_per_patient <- function(barcodes, what, fix = "rerun from 02") {
  p   <- substr(barcodes, 1, 12)
  dup <- unique(p[duplicated(p)])
  if (length(dup) > 0) {
    stop(what, ": ", length(dup), " patient(s) appear more than once (",
         paste(head(dup, 10), collapse = ", "), if (length(dup) > 10) ", ...",
         "). One aliquot per patient is required -- ", fix, ".", call. = FALSE)
  }
  invisible(TRUE)
}

run_alignment_checks <- function(m_mat, meth_clinical) {
  stopifnot("meth_clinical has no meth_sample_barcode column -- rerun 07" =
              "meth_sample_barcode" %in% names(meth_clinical))
  assert_one_per_patient(colnames(m_mat), "methylation matrix columns", "rerun 07")
  assert_one_per_patient(meth_clinical$meth_sample_barcode, "meth_clinical rows", "rerun 07")

  match_idx <- match(meth_clinical$meth_sample_barcode, colnames(m_mat))
  stopifnot(
    "Methylation and clinical do not share the same sample set" =
      !anyNA(match_idx) && length(match_idx) == ncol(m_mat))

  cat("Alignment check passed. Reorder was",
      if (identical(match_idx, seq_along(match_idx))) "a no-op"
      else "necessary and has been applied",
      "\n")

  match_idx
}

uni_screen <- function(g, surv_obj, age, stage, purity = NULL, subtype = NULL) {
  terms <- c("scale(g)", "age", "stage",
             if (!is.null(purity))  "scale(purity)",
             if (!is.null(subtype)) "strata(subtype)")
  fit <- try(coxph(reformulate(terms, response = "surv_obj")), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

load_subtype_purity <- function(meth_clinical, barcode_col = "sample_barcode",
                                pam50_path  = SUBTYPE_PURITY_INPUTS[["pam50"]],
                                purity_path = SUBTYPE_PURITY_INPUTS[["purity"]],
                                need_subtype = TRUE, need_purity = TRUE) {
  stopifnot("barcode_col not found in meth_clinical" = barcode_col %in% names(meth_clinical))
  barcodes <- meth_clinical[[barcode_col]]

  subtype_loaded <- file.exists(pam50_path)
  purity_loaded  <- file.exists(purity_path)
  if (need_subtype && !subtype_loaded) {
    stop("PAM50 file not found: ", pam50_path, call. = FALSE)
  }
  if (need_purity && !purity_loaded) {
    stop("Purity file not found: ", purity_path, " (Aran et al. 2015, ",
         "Nature Communications 6:8971, Supplementary Data 1)", call. = FALSE)
  }

  if (subtype_loaded) {
    subtype_tab <- read.delim(pam50_path, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot("PAM50 file should have 2 columns" = ncol(subtype_tab) == 2)
    names(subtype_tab) <- c("barcode", "subtype")
    subtype_tab$patient <- substr(subtype_tab$barcode, 1, 12)
    cat("PAM50 rows:", nrow(subtype_tab), "| unique patients:",
        length(unique(subtype_tab$patient)), "\n")

    dup_patients <- unique(subtype_tab$patient[duplicated(subtype_tab$patient)])
    if (length(dup_patients) > 0) {
      disagree <- vapply(dup_patients, function(p) {
        length(unique(subtype_tab$subtype[subtype_tab$patient == p])) > 1
      }, logical(1))
      cat("Duplicate participants:", length(dup_patients),
          "| with conflicting subtype calls:", sum(disagree), "\n")
      if (any(disagree)) {
        cat("CONFLICTING PAM50 calls -- inspect before trusting these patients:\n")
        print(subtype_tab[subtype_tab$patient %in% dup_patients[disagree], ])
      }
    }
    subtype_tab <- subtype_tab[!duplicated(subtype_tab$patient), ]
    meth_clinical$subtype <- factor(
      subtype_tab$subtype[match(substr(barcodes, 1, 12), subtype_tab$patient)])
  } else {
    cat("No PAM50 file at", pam50_path, "-- subtype set to NA.\n")
    meth_clinical$subtype <- factor(rep(NA_character_, nrow(meth_clinical)))
  }

  if (purity_loaded) {
    purity_tab <- as.data.frame(readxl::read_excel(purity_path, sheet = "Supp Data 1", skip = 2))
    names(purity_tab)[names(purity_tab) == "Sample ID"] <- "sample_barcode"
    purity_tab$CPE <- as.numeric(purity_tab$CPE)
    purity_tab <- purity_tab[purity_tab$`Cancer type` == "BRCA" &
                               !is.na(purity_tab$`Cancer type`), ]
    purity_tab$sample_id_15 <- substr(purity_tab$sample_barcode, 1, 15)
    meth_clinical$purity <- purity_tab$CPE[
      match(substr(barcodes, 1, 15), purity_tab$sample_id_15)]
  } else {
    cat("No purity file at", purity_path, "-- purity set to NA.\n")
    meth_clinical$purity <- NA_real_
  }

  attr(meth_clinical, "subtype_loaded") <- subtype_loaded
  attr(meth_clinical, "purity_loaded")  <- purity_loaded
  meth_clinical
}

naive_screen_one <- function(g, surv_x, meth_x) {
  fit <- try(coxph(surv_x ~ scale(g) + meth_x$age_years + meth_x$stage +
                      scale(meth_x$purity)), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(coef = NA, z = NA, p = NA))
  s <- summary(fit)$coefficients[1, ]
  c(coef = unname(s[1]), z = unname(s[4]), p = unname(s[5]))
}

extract_probes <- function(fit, s, clin_colnames) {
  coefs        <- as.matrix(coef(fit, s = s))
  selected_all <- rownames(coefs)[coefs[, 1] != 0]
  selected_all[!selected_all %in% clin_colnames]
}

load_geo_series_matrix <- function(path) {
  con <- gzfile(path, "r")
  on.exit(close(con))
  head_lines <- character(0)
  repeat {
    l <- readLines(con, n = 1)
    stopifnot("No !series_matrix_table_begin line found" = length(l) == 1)
    if (l == "!series_matrix_table_begin") break
    head_lines <- c(head_lines, l)
  }

  fields <- function(l) gsub('^"|"$', "", strsplit(l, "\t", fixed = TRUE)[[1]][-1])
  s_lines <- grep("^!Sample_", head_lines, value = TRUE)
  tags    <- sub("\t.*", "", sub("^!Sample_", "", s_lines))
  single  <- tags != "characteristics_ch1" & !duplicated(tags)

  pheno <- as.data.frame(lapply(s_lines[single], fields), col.names = tags[single],
                         stringsAsFactors = FALSE)
  for (l in s_lines[tags == "characteristics_ch1"]) {
    v   <- fields(l)
    key <- unique(sub(":.*", "", v[nzchar(v)]))
    stopifnot("A characteristics row mixes keys" = length(key) == 1)
    val <- ifelse(nzchar(v), sub("^[^:]*: ?", "", v), NA_character_)
    pheno[[key]] <- ifelse(val %in% "#N/A!", NA_character_, val)
  }

  tab <- withCallingHandlers(
    data.table::fread(path, skip = "\"ID_REF\"", header = TRUE,
                      na.strings = c("", "NA", "null", "NULL"),
                      colClasses = list(character = "ID_REF")),
    warning = function(w) if (grepl("Discarded single-line footer: <<!series_matrix_table_end>>",
                                    conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning"))
  tab  <- tab[!startsWith(tab$ID_REF, "!"), ]
  beta <- as.matrix(tab[, -1])
  rownames(beta) <- tab$ID_REF

  stopifnot("Sample order differs between the table and the header" =
              identical(colnames(beta), pheno$geo_accession),
            "Duplicate probe IDs" = !anyDuplicated(rownames(beta)),
            "No !Sample_data_row_count to check the table against" =
              !is.null(pheno$data_row_count),
            "Row count differs from !Sample_data_row_count" =
              all(as.integer(pheno$data_row_count) == nrow(beta)))
  list(pheno = pheno, beta = beta)
}
