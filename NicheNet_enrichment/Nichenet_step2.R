library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")
library(GenomicRanges)
library(rtracklayer)

setwd('E:/AAA_Labwork/T cells/v3')
lr_network <- readRDS('F:/CRC/Nichenet/lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('F:/CRC/Nichenet/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS('F:/CRC/Nichenet/weighted_networks_nsga2r_final.rds')
gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
# receptors =  lr_network %>% pull(to) %>% unique()
# write.table(receptors, 'receptors.csv')

# the total expressed genes in all donor all sites
background_expressed_genes <- read.table('TRM_background_genes.csv', sep = ',', header = 1, row.names = 1)                                                                                                                                                                                                                                                                                                                                                                                    

filelist = c('CD4DEG/IEL_vs_LP.csv','CD4DEG/LP_vs_L.csv','CD8DEG/IEL_vs_LP.csv','CD8DEG/LP_vs_L.csv')

for (file in filelist) {
  extracted <- sub("^(.*)\\.csv$", "\\1", file)
  extracted = sub('/','_', extracted)
  # read in the minimally expressed genes for each subset to determine what receptor should we consider
  expressed_receptors = read.table(paste0(extracted, '_1_receptors.csv'),sep = ',',skip = 1)$V2
  potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
  
  DE_table_receiver = read.table(file, sep = ',', header = 1, row.names = 1)
  geneset_oi <- DE_table_receiver %>% filter(q_value <= 0.05 & normalized_effect >= 0.5) %>% pull(gene_short_name)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes$X0,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
  sorted_table <- ligand_activities[order(ligand_activities$auroc,decreasing = TRUE), ]
  # 1 is for A in A vs B
  write.table(sorted_table, paste0(extracted,'_nichenet_1.csv'), row.names = FALSE)
  
  top_5_ligands = sorted_table$test_ligand[1:5]
  target_gene_union <- as.data.frame( matrix(NA, nrow = 15, ncol = 5))
  colnames(target_gene_union) <- top_5_ligands
  for (ligand in top_5_ligands) {
    top_15_targets =ligand_target_matrix[,ligand][order(ligand_target_matrix[,ligand],decreasing = TRUE) ][1:15]
    target_gene_union[ligand] = names(top_15_targets)
  }
  write.table(target_gene_union, paste0(extracted,'_nichenet_1_targets.csv'), row.names = FALSE, sep = '\t')
  
  expressed_receptors = read.table(paste0(extracted, '_2_receptors.csv'),sep = ',',skip = 1)$V2
  potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
  DE_table_receiver = read.table(file, sep = ',', header = 1, row.names = 1)
  geneset_oi <- DE_table_receiver %>% filter(q_value <= 0.05 & normalized_effect <= -0.5) %>% pull(gene_short_name)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes$X0,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
  sorted_table <- ligand_activities[order(ligand_activities$auroc,decreasing = TRUE), ]
  
  write.table(sorted_table, paste0(extracted,'_nichenet_2.csv'), row.names = FALSE)
  
  top_5_ligands = sorted_table$test_ligand[1:5]
  target_gene_union <- as.data.frame( matrix(NA, nrow = 15, ncol = 5))
  colnames(target_gene_union) <- top_5_ligands
  for (ligand in top_5_ligands) {
    top_15_targets =ligand_target_matrix[,ligand][order(ligand_target_matrix[,ligand],decreasing = TRUE) ][1:15]
    target_gene_union[ligand] = names(top_15_targets)
  }
  write.table(target_gene_union, paste0(extracted,'_nichenet_2_targets.csv'), row.names = FALSE, sep = '\t')
}