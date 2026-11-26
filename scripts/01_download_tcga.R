library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(here)

# Load TCGA breast cancer data
brca_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       sample.type = "Primary Tumor")

dl_dir <- "00_raw_data/GDCdata"

GDCdownload(brca_query, directory = dl_dir, method = "client", files.per.chunk = 30)

brca_data <- GDCprepare(brca_query, directory = dl_dir)

saveRDS(brca_data, "00_raw_data/TCGA_BRCA_SummarizedExperiment.rds")

# downloading normal data

normal_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

GDCdownload(normal_query, directory = dl_dir, method = "client", files.per.chunk = 30)

normal_data <- GDCprepare(normal_query, directory = dl_dir)

saveRDS(normal_data, "00_raw_data/TCGA_BRCA_Normal_SummarizedExperiment.rds")

