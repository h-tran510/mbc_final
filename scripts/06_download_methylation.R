library(TCGAbiolinks)
library(SummarizedExperiment)
library(sesame)

dl_dir <- "00_raw_data/GDCdata"

# Query TCGA-BRCA HM450 methylation beta values for primary tumors.
# 793 files (~1-2 GB total). Runtime depends on connection speed.

query_meth <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  data.type     = "Methylation Beta Value",
  platform      = "Illumina Human Methylation 450",
  sample.type   = "Primary Tumor"
)

cat("Files to download:", nrow(getResults(query_meth)), "\n")

GDCdownload(query_meth, directory = dl_dir, method = "client", files.per.chunk = 30)

meth_se <- GDCprepare(query_meth, directory = dl_dir)

cat("Methylation SE dimensions:", dim(meth_se), "\n")
cat("Assay names:", assayNames(meth_se), "\n")

saveRDS(meth_se, "00_raw_data/TCGA_BRCA_Methylation450_SE.rds")

cat("Saved to 00_raw_data/TCGA_BRCA_Methylation450_SE.rds\n")
