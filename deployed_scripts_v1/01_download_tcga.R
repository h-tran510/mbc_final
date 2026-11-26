library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(here)

source("deployed_scripts_v1/00_utils.R")

dl_dir <- "00_raw_data/GDCdata"

brca_data <- gdc_cached("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds", function() {
brca_query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       sample.type = "Primary Tumor")

GDCdownload(brca_query, directory = dl_dir, method = "client", files.per.chunk = 30)

GDCprepare(brca_query, directory = dl_dir)
})

normal_data <- gdc_cached("00_raw_data/TCGA_BRCA_Normal_SummarizedExperiment.rds", function() {
normal_query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)

GDCdownload(normal_query, directory = dl_dir, method = "client", files.per.chunk = 30)

GDCprepare(normal_query, directory = dl_dir)
})
