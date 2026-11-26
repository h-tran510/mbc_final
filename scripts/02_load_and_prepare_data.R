library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

brca_data <- readRDS("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds")

colData(brca_data)

clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

head(clinical)
colnames(clinical)

sample_barcodes <- colnames(brca_data)
patient_ids <- substr(sample_barcodes, 1, 12)

metadata <- data.frame(
  sample_barcode = sample_barcodes,
  patient_id = patient_ids
)

# merged dataset with patient barcodes

merged_data <- merge(metadata, clinical, by.x = "patient_id", by.y = "submitter_id")

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

# extract expression matrix 

counts <- assay(brca_data)

# align counts with metadata

counts <- counts[, merged_data$sample_barcode]

# save into processed data

saveRDS(merged_data, "01_processed_data/merged_data.rds")




