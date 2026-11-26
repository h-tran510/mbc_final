library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

if (as.integer(Sys.getenv("MBC_SUBSET_N", "0")) > 0) {
  merged_data <- readRDS("01_processed_data/merged_data.rds")
  brca_data   <- readRDS("00_raw_data/TCGA_BRCA_SummarizedExperiment.rds")
  counts      <- assay(brca_data)[, merged_data$sample_barcode]
}

keep <- rowSums(counts >= 10) >= (0.2 * ncol(counts))
counts <- counts[keep, ]

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = merged_data,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)

vsd <- vst(dds, blind = TRUE)
expr_matrix <- assay(vsd)

saveRDS(expr_matrix, "01_processed_data/expr_matrix_vst.rds")
