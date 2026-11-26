library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)

source("deployed_scripts_v1/00_utils.R")

dl_dir <- "00_raw_data/GDCdata"

meth_se <- gdc_cached("00_raw_data/TCGA_BRCA_Methylation450_SE.rds", function() {
query_meth <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  sample.type   = "Primary Tumor"
)

cat("Files to download:", nrow(getResults(query_meth)), "\n")

GDCdownload(query_meth, directory = dl_dir, method = "client", files.per.chunk = 30)

GDCprepare(query_meth, directory = dl_dir)
})

cat("Methylation SE dimensions:", dim(meth_se), "\n")
cat("Assay names:", assayNames(meth_se), "\n")
