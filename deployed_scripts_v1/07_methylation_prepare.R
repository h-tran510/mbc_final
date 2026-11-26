library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readxl)

source("deployed_scripts_v1/00_utils.R")

meth_se     <- readRDS("00_raw_data/TCGA_BRCA_Methylation450_SE.rds")
merged_data <- readRDS("01_processed_data/merged_data.rds")
cdr_raw <- read_excel("00_raw_data/TCGAsupp/TCGA-CDR-SupplementalTableS1.xlsx")

cat("Raw methylation SE:", nrow(meth_se), "probes x", ncol(meth_se), "samples\n")

cdr_breast <- cdr_raw |>
  filter(type == "BRCA")
n_brca <- nrow(cdr_breast)
cdr_breast <- cdr_breast |>
  filter(is.na(Redaction)) |>
  select(bcr_patient_barcode, PFI, PFI.time)
cat("CDR BRCA rows dropped as Redacted:", n_brca - nrow(cdr_breast), "\n")
meth_barcodes    <- colnames(meth_se)
meth_patient_ids <- substr(meth_barcodes, 1, 12)

assert_one_per_patient(merged_data$sample_barcode, "merged_data")

meth_meta <- data.frame(
  meth_sample_barcode = meth_barcodes,
  patient_id          = meth_patient_ids,
  stringsAsFactors = FALSE
)

meth_clinical <- as.data.frame(
  merge(meth_meta, merged_data, by = "patient_id")
)

meth_clinical <- meth_clinical[order(meth_clinical$meth_sample_barcode), ]
meth_clinical <- meth_clinical[!duplicated(meth_clinical$patient_id), ]

cat("Patients matched to clinical data:", nrow(meth_clinical), "\n")

meth_clinical$age_years <- meth_clinical$age_at_diagnosis / 365.25

stage_map <- c(
  "Stage I"    = "I",   "Stage IA"   = "I",   "Stage IB"   = "I",
  "Stage II"   = "II",  "Stage IIA"  = "II",  "Stage IIB"  = "II",
  "Stage III"  = "III", "Stage IIIA" = "III", "Stage IIIB" = "III",
  "Stage IIIC" = "III",
  "Stage IV"   = "IV"
)
meth_clinical$stage <- factor(
  stage_map[meth_clinical$ajcc_pathologic_stage],
  levels = c("I", "II", "III", "IV")
)

meth_clinical <- left_join(meth_clinical, cdr_breast, by = "bcr_patient_barcode")

meth_clinical$time <- meth_clinical$PFI.time
meth_clinical$event <- meth_clinical$PFI

valid <- !is.na(meth_clinical$time) & meth_clinical$time > 0 &
         !is.na(meth_clinical$event) &
         !is.na(meth_clinical$age_years) &
         !is.na(meth_clinical$stage)

no_cdr_pfi <- unique(meth_clinical$bcr_patient_barcode[is.na(meth_clinical$PFI)])
meth_clinical <- meth_clinical[valid, ]
cat("Matched patients removed for having no CDR PFI (incl. Redacted):",
    length(no_cdr_pfi), "\n")

cat("Samples retained after clinical filtering:", nrow(meth_clinical), "\n")
cat("PFI events:", sum(meth_clinical$event), "\n")
cat("Median follow-up (days):",
    median(meth_clinical$time[meth_clinical$event == 0]), "\n")

meth_se_clean <- meth_se[, meth_clinical$meth_sample_barcode]

assay_names <- names(assays(meth_se_clean))

cat("Available assays:\n")
print(assay_names)

beta_assays <- intersect(c("Beta_value", "Beta_value_mean", "beta", "Beta"), assay_names)
if (length(beta_assays) > 0) {
  beta_mat <- assay(meth_se_clean, beta_assays[1])
} else {
  cat("No named beta-value assay found. Using the first assay by index.\n")
  beta_mat <- assay(meth_se_clean, 1)
}

cat("Beta matrix dimensions:", dim(beta_mat), "\n")

miss_frac <- rowMeans(is.na(beta_mat))
beta_mat  <- beta_mat[miss_frac <= 0.2, ]
cat("Probes after missing-value filter (>20% NA removed):", nrow(beta_mat), "\n")

na_rows <- which(rowSums(is.na(beta_mat)) > 0)
if (length(na_rows) > 0) {
  row_medians <- apply(beta_mat[na_rows, , drop = FALSE], 1, median, na.rm = TRUE)
  for (i in seq_along(na_rows)) {
    probe_row <- na_rows[i]
    missing_samples <- is.na(beta_mat[probe_row, ])
    beta_mat[probe_row, missing_samples] <- row_medians[i]
  }
  cat("Probes with imputed missing values:", length(na_rows), "\n")
}

probe_chr <- as.character(seqnames(rowRanges(meth_se_clean)[rownames(beta_mat)]))
sex_probes <- probe_chr %in% c("chrX", "chrY")
beta_mat   <- beta_mat[!sex_probes, ]
cat("Probes after sex chromosome removal:", nrow(beta_mat), "\n")

probe_var <- row_vars(beta_mat)
var_thresh <- quantile(probe_var, 0.9)
beta_filt  <- beta_mat[probe_var >= var_thresh, ]
cat("Probes after top-10% variance filter:", nrow(beta_filt), "\n")

beta_clipped <- pmin(pmax(beta_filt, 0.001), 0.999)
m_mat        <- log2(beta_clipped / (1 - beta_clipped))

cat("M-value matrix dimensions:", dim(m_mat), "\n")
cat("M-value range:", round(range(m_mat), 2), "\n")

stopifnot(identical(colnames(m_mat), meth_clinical$meth_sample_barcode))

saveRDS(m_mat,         "01_processed_data/meth_m_matrix.rds")
saveRDS(meth_clinical, "01_processed_data/meth_clinical.rds")

write_manifest("01_processed_data", list(n_probes = nrow(m_mat), n_samples = ncol(m_mat)),
               inputs  = c("00_raw_data/TCGA_BRCA_Methylation450_SE.rds",
                           "01_processed_data/merged_data.rds",
                           "00_raw_data/TCGAsupp/TCGA-CDR-SupplementalTableS1.xlsx"),
               outputs = c("01_processed_data/meth_m_matrix.rds",
                           "01_processed_data/meth_clinical.rds"),
               script  = "07_methylation_prepare.R",
               file    = "run_manifest_07_methylation_prepare.txt")

cat("\nSaved:\n")
cat("  01_processed_data/meth_m_matrix.rds\n")
cat("  01_processed_data/meth_clinical.rds\n")
