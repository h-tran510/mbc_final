library(survival)
library(SummarizedExperiment)

# remove genes with very low variance across samples
gene_var <- apply(expr_matrix, 1, var, na.rm = TRUE)

summary(gene_var)

var_threshold <- quantile(gene_var, 0.5, na.rm = TRUE)
expr_matrix_filtered <- expr_matrix[gene_var > var_threshold, ]

dim(expr_matrix_filtered)

# make survival object
surv_obj <- Surv(time = merged_data$time, event = merged_data$event)

# run univariate cox models
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

# remove failed and N/A rows
results <- results[complete.cases(results), ]

dim(results)
head(results)

# multiple testing correction
results$adj_pvalue <- p.adjust(results$pvalue, method = "BH")

# filter for significant genes

significant_genes <- results[
  results$adj_pvalue < 0.05 &
    abs(log(results$HR)) > log(1.2),
]

# sort by adjusted p-value
significant_genes <- significant_genes[order(significant_genes$adj_pvalue), ]

dim(significant_genes)
head(significant_genes)

# top genes
top_genes <- significant_genes[1:min(50, nrow(significant_genes)), ]
top_genes

# save outputs

saveRDS(results, "02_analysis/univariate_cox_results.rds")
saveRDS(significant_genes, "02_analysis/significant_genes_cox.rds")
write.csv(results, "02_analysis/univariate_cox_results.csv", row.names = FALSE)
write.csv(significant_genes, "02_analysis/significant_genes_cox.csv", row.names = FALSE)

merged_data <- readRDS("01_processed_data/merged_data.rds")