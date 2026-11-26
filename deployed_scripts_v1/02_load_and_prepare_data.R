library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

source("deployed_scripts_v1/00_utils.R")

brca_data <- readRDS("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds")

colData(brca_data)

ffpe <- if ("preservation_method" %in% names(colData(brca_data))) {
  print(table(preservation_method = brca_data$preservation_method, useNA = "ifany"))
  brca_data$preservation_method %in% "FFPE"
} else {
  message("WARNING: colData has no preservation_method: the FFPE step of the ",
          "aliquot rule is skipped (every aliquot treated as frozen)")
  rep(FALSE, ncol(brca_data))
}
cat("RNA aliquots flagged FFPE:", sum(ffpe), "\n")
aliquots <- pick_one_aliquot(colnames(brca_data), ffpe)
cat("Rule:", RNA_ALIQUOT_RULE, "\n")
cat("RNA samples:", ncol(brca_data), "| patients:", sum(aliquots$keep),
    "| multi-aliquot patients:", length(unique(aliquots$choice$patient_id)),
    "| aliquots dropped:", sum(!aliquots$keep), "\n")
print(aliquots$choice)
dir.create("02_analysis", showWarnings = FALSE)
write.csv(aliquots$choice, "02_analysis/rna_aliquot_choice.csv", row.names = FALSE)

clinical <- gdc_cached("00_raw_data/gdc_clinical.rds", function() {
  GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
})

head(clinical)
colnames(clinical)

sample_barcodes <- colnames(brca_data)[aliquots$keep]
patient_ids <- substr(sample_barcodes, 1, 12)

metadata <- data.frame(
  sample_barcode = sample_barcodes,
  patient_id = patient_ids
)

merged_data <- merge(metadata, clinical, by.x = "patient_id", by.y = "submitter_id")
assert_one_per_patient(merged_data$sample_barcode, "merged_data",
                       "check gdc_clinical.rds for a repeated submitter_id")

merged_data$time <- ifelse(
  is.na(merged_data$days_to_death),
  merged_data$days_to_last_follow_up,
  merged_data$days_to_death
)

merged_data$event <- ifelse(
  merged_data$vital_status == "Dead",
  1,
  0
)

counts <- assay(brca_data)

counts <- counts[, merged_data$sample_barcode]

saveRDS(merged_data, "01_processed_data/merged_data.rds")

write_manifest("01_processed_data", list(n_rna_samples = ncol(brca_data),
                                         n_patients = nrow(merged_data)),
               inputs  = c("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds",
                           "00_raw_data/gdc_clinical.rds"),
               outputs = c("02_analysis/rna_aliquot_choice.csv",
                           "01_processed_data/merged_data.rds",
                           written_since("00_raw_data/gdc_clinical.rds")),
               script  = "02_load_and_prepare_data.R",
               file    = "run_manifest_02_load_and_prepare_data.txt")
