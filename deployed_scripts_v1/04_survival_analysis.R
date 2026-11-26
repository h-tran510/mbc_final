library(survival)
library(SummarizedExperiment)

source("deployed_scripts_v1/00_utils.R")

if (as.integer(Sys.getenv("MBC_SUBSET_N", "0")) > 0) {
  merged_data <- readRDS("01_processed_data/merged_data.rds")
  expr_matrix <- readRDS("01_processed_data/expr_matrix_vst.rds")
}
assert_one_per_patient(merged_data$sample_barcode, "merged_data")
stopifnot("VST columns are not merged_data$sample_barcode in order -- rerun from 02" =
            identical(colnames(expr_matrix), merged_data$sample_barcode))

gene_var <- apply(expr_matrix, 1, var, na.rm = TRUE)

summary(gene_var)

var_threshold <- quantile(gene_var, 0.5, na.rm = TRUE)
expr_matrix_filtered <- expr_matrix[gene_var > var_threshold, ]

if (as.integer(Sys.getenv("MBC_SUBSET_N", "0")) > 0) {
  top <- order(gene_var[rownames(expr_matrix_filtered)], decreasing = TRUE)
  expr_matrix_filtered <- expr_matrix_filtered[head(top, 2000), ]
}

dim(expr_matrix_filtered)

surv_obj <- Surv(time = merged_data$time, event = merged_data$event)

gene_names <- rownames(expr_matrix_filtered)

results_list <- vector("list", nrow(expr_matrix_filtered))

for (i in seq_len(nrow(expr_matrix_filtered))) {
  gene_expr <- as.numeric(expr_matrix_filtered[i, ])
  model <- try(coxph(surv_obj ~ gene_expr), silent = TRUE)
  if (inherits(model, "try-error")) {
    results_list[[i]] <- data.frame(
      gene = gene_names[i],
      coef = NA,
      HR = NA,
      pvalue = NA,
      concordance = NA
    )
    next
  }
  s <- summary(model)
  coef_val <- s$coef[1, "coef"]
  hr_val   <- s$coef[1, "exp(coef)"]
  p_val    <- s$coef[1, "Pr(>|z|)"]
  c_index  <- s$concordance[1]
  if (!is.finite(coef_val) || !is.finite(hr_val) || !is.finite(p_val)) {
    results_list[[i]] <- data.frame(
      gene = gene_names[i],
      coef = NA,
      HR = NA,
      pvalue = NA,
      concordance = NA
    )
  } else {
    results_list[[i]] <- data.frame(
      gene = gene_names[i],
      coef = coef_val,
      HR = hr_val,
      pvalue = p_val,
      concordance = c_index
    )
  }
}

results <- do.call(rbind, results_list)

results <- results[complete.cases(results), ]

dim(results)
head(results)

results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

significant_genes <- results[
  results$adj_pvalue < 0.05 &
    abs(log(results$HR)) > log(1.2),
]

significant_genes <- significant_genes[order(significant_genes$adj_pvalue), ]

dim(significant_genes)
head(significant_genes)

top_genes <- significant_genes[1:min(50, nrow(significant_genes)), ]
top_genes

saveRDS(results, "02_analysis/univariate_cox_results.rds")
saveRDS(significant_genes, "02_analysis/significant_genes_cox.rds")
write.csv(results, "02_analysis/univariate_cox_results.csv", row.names = FALSE)
write.csv(significant_genes, "02_analysis/significant_genes_cox.csv", row.names = FALSE)

merged_data <- readRDS("01_processed_data/merged_data.rds")
