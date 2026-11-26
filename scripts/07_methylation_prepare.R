library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(readxl)

# Loads the HM450 methylation SummarizedExperiment, loads in TCGA-CDR data,
# Aligns samples to clinical, merges TCGA-CDR with TCGA BRCA,
# Runs probe-level QC, applies M-value transformation, and saves the
# Processed matrix and matched clinical dataframe for use in survival model.

# ── 1. Load raw data ──────────────────────────────────────────────────────────
meth_se     <- readRDS("00_raw_data/TCGA_BRCA_Methylation450_SE.rds")
merged_data <- readRDS("01_processed_data/merged_data.rds")
cdr_raw <- read_excel("00_raw_data/TCGAsupp/TCGA-CDR-SupplementalTableS1.xlsx")

cat("Raw methylation SE:", nrow(meth_se), "probes x", ncol(meth_se), "samples\n")

# ── 2. Align methylation samples with clinical data ───────────────────────────

# Filter by BRCA into object cdr_breast
cdr_breast <- cdr_raw |>
  filter(type == "BRCA") |>
  select(bcr_patient_barcode, PFI, PFI.time)
meth_barcodes    <- colnames(meth_se)
meth_patient_ids <- substr(meth_barcodes, 1, 12)

# One clinical row per patient (deduplicated from the RNA-seq merged_data)
clinical_dedup <- merged_data[!duplicated(merged_data$patient_id), ]

# Rename methylation sample barcode to avoid merge() creating sample_barcode.x/y
meth_meta <- data.frame(
  meth_sample_barcode = meth_barcodes,
  patient_id          = meth_patient_ids,
  stringsAsFactors = FALSE
)

meth_clinical <- as.data.frame(
  merge(meth_meta, clinical_dedup, by = "patient_id")
)

# One methylation sample per patient — keep the lexicographically first barcode
meth_clinical <- meth_clinical[order(meth_clinical$meth_sample_barcode), ]
meth_clinical <- meth_clinical[!duplicated(meth_clinical$patient_id), ]

cat("Patients matched to clinical data:", nrow(meth_clinical), "\n")

# ── 3. Clinical covariate preparation ────────────────────────────────────────
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

# Left-join cdr_breast with clinical covariates
meth_clinical <- left_join(meth_clinical, cdr_breast, by = "bcr_patient_barcode")

# Override time and event with PFI.time and PFI
meth_clinical$time <- meth_clinical$PFI.time
meth_clinical$event <- meth_clinical$PFI

# ── 4. Exclude samples with invalid survival data or missing covariates ───────
# time <= 0: same-day or data entry errors (mirrors RNA-seq pipeline exclusion)
valid <- !is.na(meth_clinical$time) & meth_clinical$time > 0 &
         !is.na(meth_clinical$event) &
         !is.na(meth_clinical$age_years) &
         !is.na(meth_clinical$stage)

meth_clinical <- meth_clinical[valid, ]

cat("Samples retained after clinical filtering:", nrow(meth_clinical), "\n")
cat("PFI events:", sum(meth_clinical$event), "\n")
cat("Median follow-up (days):",
    median(meth_clinical$time[meth_clinical$event == 0]), "\n")

# Subset SE to retained, matched samples
meth_se_clean <- meth_se[, meth_clinical$meth_sample_barcode]

# ── 5. Extract beta value matrix ──────────────────────────────────────────────
assay_names <- names(assays(meth_se_clean))

cat("Available assays:\n")
print(assay_names)

if (!is.null(assay_names) && "Beta_value" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "Beta_value")

} else if (!is.null(assay_names) && "Beta_value_mean" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "Beta_value_mean")

} else if (!is.null(assay_names) && "beta" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "beta")

} else if (!is.null(assay_names) && "Beta" %in% assay_names) {
  beta_mat <- assay(meth_se_clean, "Beta")

} else {
  cat("No named beta-value assay found. Using the first assay by index.\n")
  beta_mat <- assay(meth_se_clean, 1)
}

cat("Beta matrix dimensions:", dim(beta_mat), "\n")

# ── 6. Probe-level QC filtering ───────────────────────────────────────────────

# (a) Remove probes missing in more than 20% of samples
miss_frac <- rowMeans(is.na(beta_mat))
beta_mat  <- beta_mat[miss_frac <= 0.2, ]
cat("Probes after missing-value filter (>20% NA removed):", nrow(beta_mat), "\n")

# (b) Impute remaining missing values with each probe's row median
na_rows <- which(rowSums(is.na(beta_mat)) > 0)
if (length(na_rows) > 0) {
  row_medians <- apply(beta_mat[na_rows, , drop = FALSE], 1, median, na.rm = TRUE)
  for (i in seq_along(na_rows)) {
    idx <- na_rows[i]
    beta_mat[idx, is.na(beta_mat[idx, ])] <- row_medians[i]
  }
  cat("Probes with imputed missing values:", length(na_rows), "\n")
}

# (c) Remove sex chromosome probes to prevent sex from driving the risk score
probe_chr <- as.character(seqnames(rowRanges(meth_se_clean)[rownames(beta_mat)]))
sex_probes <- probe_chr %in% c("chrX", "chrY")
beta_mat   <- beta_mat[!sex_probes, ]
cat("Probes after sex chromosome removal:", nrow(beta_mat), "\n")

# (d) Variance filter — retain top 10% most variable probes
probe_var  <- apply(beta_mat, 1, var)
var_thresh <- quantile(probe_var, 0.9)
beta_filt  <- beta_mat[probe_var >= var_thresh, ]
cat("Probes after top-10% variance filter:", nrow(beta_filt), "\n")

# ── 7. M-value transformation ─────────────────────────────────────────────────
# M = log2(beta / (1 - beta)). Beta clipped to (0.001, 0.999) to avoid ±Inf.
# M-values are more statistically well-behaved than beta values for regression:
# they are approximately normal and homoscedastic.
beta_clipped <- pmin(pmax(beta_filt, 0.001), 0.999)
m_mat        <- log2(beta_clipped / (1 - beta_clipped))

cat("M-value matrix dimensions:", dim(m_mat), "\n")
cat("M-value range:", round(range(m_mat), 2), "\n")

# Confirm column order of m_mat matches row order of meth_clinical
stopifnot(identical(colnames(m_mat), meth_clinical$meth_sample_barcode))

# ── 8. Save processed outputs ─────────────────────────────────────────────────
saveRDS(m_mat,         "01_processed_data/meth_m_matrix.rds")
saveRDS(meth_clinical, "01_processed_data/meth_clinical.rds")

cat("\nSaved:\n")
cat("  01_processed_data/meth_m_matrix.rds\n")
cat("  01_processed_data/meth_clinical.rds\n")