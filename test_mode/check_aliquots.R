# =============================================================================
# test_mode/check_aliquots.R
#
# Test-mode fixture and check for the one-aliquot-per-patient rule
# (RNA_ALIQUOT_RULE / pick_one_aliquot() in deployed_scripts_v1/00_utils.R).
# The sample dataset has one RNA aliquot per patient, so test_mode/run_parallel.sh
# plants extra aliquots before 02 and checks 02's choice after it. Runs inside
# _test_run/ only (it rewrites that copy of the RNA SE):
#
#   Rscript test_mode/check_aliquots.R plant    # before 02
#   Rscript test_mode/check_aliquots.R verify   # after 02
#
# plant: for three patients, adds copies of their RNA column under new barcodes,
# placed BEFORE all original columns, so a "first in row order wins" rule would
# pick them:
#   patient 1  + 01A, FFPE, plate ZZZZ          -> original kept (frozen before FFPE,
#                                                  although the copy has the higher plate)
#   patient 2  + 01A, frozen, plate ZZZZ        -> the planted one kept (highest plate)
#   patient 3  + 01A, frozen, plate 0000, and
#              + 01C, frozen, plate ZZZZ        -> original kept (vial A, then highest plate)
# The three originals are the first frozen 01A aliquots in barcode order. The copies carry
# the same counts, so no downstream number changes.
# verify: merged_data has one row per patient, 02's choice table equals the
# expected picks, and the guards stop on duplicates.
# =============================================================================

suppressMessages(library(SummarizedExperiment))
source("deployed_scripts_v1/00_utils.R")

mode <- commandArgs(trailingOnly = TRUE)[1]
se_path   <- "00_raw_data/TCGA_BRCA_SummarizedExperiment.rds"
plan_path <- "run_logs/planted_aliquots.csv"
stopifnot("run inside _test_run/ in test mode (MBC_SUBSET_N > 0)" =
            TEST_MODE && basename(getwd()) == "_test_run")

with_field <- function(b, field, value) {
  f <- strsplit(b, "-", fixed = TRUE)[[1]]
  f[field] <- value
  paste(f, collapse = "-")
}

if (identical(mode, "plant")) {
  se <- readRDS(se_path)
  assert_one_per_patient(colnames(se), "sample RNA SE before planting",
                         "rebuild it with test_mode/make_sample_dataset.R")
  if (!"preservation_method" %in% names(colData(se))) se$preservation_method <- "Unknown"
  # originals: frozen 01A aliquots of patients 02's clinical merge keeps
  clin <- readRDS("00_raw_data/gdc_clinical.rds")
  frozen_01a <- colnames(se)[substr(colnames(se), 14, 16) == "01A" &
                               !se$preservation_method %in% "FFPE" &
                               substr(colnames(se), 1, 12) %in% clin$submitter_id]
  stopifnot("the fixture needs three frozen 01A originals" = length(frozen_01a) >= 3)
  orig <- sort(frozen_01a, method = "radix")[1:3]
  new <- c(with_field(orig[1], 4, "01A") |> with_field(6, "ZZZZ"),
           with_field(orig[2], 4, "01A") |> with_field(6, "ZZZZ"),
           with_field(orig[3], 4, "01A") |> with_field(6, "0000"),
           with_field(orig[3], 4, "01C") |> with_field(6, "ZZZZ"))
  from <- c(orig[1], orig[2], orig[3], orig[3])
  ffpe <- c(TRUE, FALSE, FALSE, FALSE)
  planted <- se[, match(from, colnames(se))]
  colnames(planted) <- new
  planted$preservation_method <- ifelse(ffpe, "FFPE", "Unknown")
  if ("barcode" %in% names(colData(planted))) planted$barcode <- new
  se <- cbind(planted, se)
  # test copy only: the first assay (the counts 02/03 read) and rowData (09) are
  # all the pipeline uses; uncompressed, so it is fast to write and read
  assays(se) <- assays(se)[1]
  saveRDS(se, se_path, compress = FALSE)
  expected <- data.frame(patient_id = substr(c(orig, from), 1, 12),
                         barcode = c(orig, new), ffpe = c(FALSE, FALSE, FALSE, ffpe))
  expected$chosen <- expected$barcode %in% c(orig[1], new[2], orig[3])
  expected <- expected[order(expected$patient_id, !expected$chosen, expected$barcode,
                             method = "radix"), ]
  rownames(expected) <- NULL
  write.csv(expected, plan_path, row.names = FALSE)
  cat("Planted", length(new), "extra aliquots for 3 patients:", ncol(se), "RNA columns\n")
  print(expected)
} else if (identical(mode, "verify")) {
  expected <- read.csv(plan_path, stringsAsFactors = FALSE)
  got      <- read.csv("02_analysis/rna_aliquot_choice.csv", stringsAsFactors = FALSE)
  merged   <- readRDS("01_processed_data/merged_data.rds")
  print(got)
  stopifnot(
    "02's aliquot choice differs from the expected picks" = identical(got, expected),
    "merged_data has a repeated patient" = !anyDuplicated(merged$patient_id),
    "merged_data does not carry the chosen aliquots" =
      all(expected$barcode[expected$chosen] %in% merged$sample_barcode) &&
      !any(expected$barcode[!expected$chosen] %in% merged$sample_barcode))
  # the guards fail loudly on a repeated patient
  # TRUE only when expr stops with the duplicate-patient message
  fails <- function(expr) tryCatch({ expr; FALSE }, error = function(e)
    grepl("appear more than once", conditionMessage(e), fixed = TRUE))
  dup <- c(expected$barcode[1:2], merged$sample_barcode[1])
  m_dup <- matrix(0, 1, 2, dimnames = list("cg1", expected$barcode[expected$patient_id ==
                                                               expected$patient_id[1]]))
  stopifnot(
    "assert_one_per_patient did not stop on a repeated patient" =
      fails(assert_one_per_patient(dup, "fixture")),
    "run_alignment_checks did not stop on a repeated patient" =
      fails(run_alignment_checks(m_dup, data.frame(meth_sample_barcode = colnames(m_dup)))))
  cat("Aliquot check passed:", nrow(merged), "patients,",
      length(unique(expected$patient_id)), "planted multi-aliquot patients resolved by the rule\n")
} else {
  stop("usage: Rscript test_mode/check_aliquots.R plant|verify")
}
