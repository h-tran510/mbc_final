library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)

# filter
keep <- rowSums(counts >= 10) >= (0.2 * ncol(counts))
counts <- counts[keep, ]

# normalize

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = merged_data,
  design = ~ 1   # no model yet, just normalization
)

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)

# transform for ml

vsd <- vst(dds, blind = TRUE)
expr_matrix <- assay(vsd)

saveRDS(expr_matrix, "01_processed_data/expr_matrix_vst.rds")

