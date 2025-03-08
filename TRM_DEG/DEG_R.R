setwd('E:/AAA_Labwork/T cells/v2/')
srat = readRDS('forDE_TRM_cd8.rds')
#srat = readRDS('forDE_TRM_cd4.rds')
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

tissue_pairs <- list(c("LP", "L"), c("IEL", "LP"))
# Initialize list to store results
results_list <- list()

# Run DE analysis for each pair of tissues
for(i in seq_along(tissue_pairs)) {
  
  # Subset data to include only the current pair of tissues
  cds_subset <- cds[,cds@colData$tissue %in% tissue_pairs[[i]]]
  
  # Run DE analysis
  gene_fits <- fit_models(cds_subset, model_formula_str = "~tissue",cores = 4)
  fit_coefs <- coefficient_table(gene_fits)
  #saveRDS(fit_coefs,paste0(paste(tissue_pairs[[i]],collapse = "_vs_"),'_fitcoef_cd4.rds'))
  
  # Adjust condition in the filter
  terms <- fit_coefs %>% filter(term == paste0("tissue", tissue_pairs[[i]][1]))
  terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)
  
  # Store results
  results_list[[paste(tissue_pairs[[i]], collapse="_vs_")]] <- terms
  write.csv(results_list[[paste(tissue_pairs[[i]], collapse="_vs_")]],paste0(paste(tissue_pairs[[i]], collapse="_vs_"),'.csv'))
}
