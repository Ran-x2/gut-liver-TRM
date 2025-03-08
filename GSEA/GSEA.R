install.packages("enrichR")
BiocManager::install("clusterProfiler")
organism = "org.Hs.eg.db"
library(enrichR)
library(organism, character.only = TRUE)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
dbs <- c("Reactome_2022","WikiPathway_2023_Human","KEGG_2021_Human","GO_Biological_Process_2023","GO_Molecular_Function_2023","GO_Cellular_Component_2023")
dbs <- c("Reactome_2022")
setwd('E:/AAA_Labwork/T cells/v2')
m_t2g_C5 <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_description, human_gene_symbol )

m_t2g_C7 <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_description, human_gene_symbol )

# CD8 IEL vs LP -------------------------------------------------------------

DE_table = as.data.frame(read.csv('CD8DEG/IEL_vs_LP.csv',row.names = 2))
prefix = 'CD8IELvsLP'
DE_table = subset(DE_table, q_value <= 0.05)
gene_list = DE_table['normalized_effect']
gene_list <- gene_list[order(gene_list$normalized_effect, decreasing = TRUE), , drop = FALSE]
rowname = row.names(gene_list)
gene_list = as.numeric(unlist(gene_list))
names(gene_list) = rowname
gene_list = gene_list[startsWith(names(gene_list),'IG')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPS')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPL')==0]
gene_list = gene_list[abs(gene_list)>0.25]
uplist = names(gene_list[gene_list>0])
downlist = names(gene_list[gene_list<0])
UPenriched <- enrichr(uplist, dbs)
DOWNenriched = enrichr(downlist, dbs)
ALLenriched = enrichr(names(gene_list), dbs)

for (DB_name in dbs) {
  sig_terms = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_UP_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
  
  sig_terms = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_DOWN_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
  
  sig_terms = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_ALL_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
}
# C7_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
# gse_df <- as.data.frame(C7_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C7_GSEAResults.csv")
# write.csv(gse_df, file = gsea_csv_filename)
# 
# C5_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g)
# gse_df <- as.data.frame(C5_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C5_GSEAResults.csv")
# write.csv(gse_df, file = gsea_csv_filename)

gse <- gseGO(geneList = gene_list,
             ont = "ALL",
             keyType = "SYMBOL", #"ENSEMBL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "BH")
gse_df <- as.data.frame(gse)
gsea_csv_filename <-  paste0(prefix,"_GSEGO_GSEAResults.csv")
write.csv(gse_df, file = gsea_csv_filename)



# CD8 LP vs L -------------------------------------------------------------
DE_table = as.data.frame(read.csv('CD8DEG/LP_vs_L.csv',row.names = 2))
prefix = 'CD8LPvsL'
DE_table = subset(DE_table, q_value <= 0.05)
gene_list = DE_table['normalized_effect']
gene_list <- gene_list[order(gene_list$normalized_effect, decreasing = TRUE), , drop = FALSE]
rowname = row.names(gene_list)
gene_list = as.numeric(unlist(gene_list))
names(gene_list) = rowname
gene_list = gene_list[startsWith(names(gene_list),'IG')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPS')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPL')==0]



# # 
# # 
# C5_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g_C5)
# gse_df_C5 <- as.data.frame(C5_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C5_GSEAResults.csv")
# write.csv(gse_df_C5, file = gsea_csv_filename)
gene_list = gene_list[abs(gene_list)>=1]
uplist = names(gene_list[gene_list>0])
downlist = names(gene_list[gene_list<0])
# UPenriched <- enrichr(uplist, dbs)
# DOWNenriched = enrichr(downlist, dbs)
# ALLenriched = enrichr(names(gene_list), dbs)
# 
# for (DB_name in dbs) {
#   sig_terms = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
#   filename <- paste0(prefix,"_",DB_name,"_UP_enrichr_significantONLY.csv")
#   write.csv(sig_terms, filename)
#   
#   sig_terms = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
#   filename <- paste0(prefix,"_",DB_name,"_DOWN_enrichr_significantONLY.csv")
#   write.csv(sig_terms, filename)
#   
#   sig_terms = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
#   filename <- paste0(prefix,"_",DB_name,"_ALL_enrichr_significantONLY.csv")
#   write.csv(sig_terms, filename)
# }

# gse <- gseGO(geneList = gene_list,
#              ont = "ALL",
#              keyType = "SYMBOL", #"ENSEMBL",
#              nPerm = 10000,
#              minGSSize = 5,
#              maxGSSize = 800,
#              pvalueCutoff = 0.05,
#              verbose = TRUE,
#              OrgDb = organism,
#              pAdjustMethod = "BH")
# 
# gse_df <- as.data.frame(gse)
# gsea_csv_filename <-  paste0(prefix,"_GSEGO_GSEAResults.csv")
# write.csv(gse_df, file = gsea_csv_filename)
C7_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g_C7)
gse_df_C7 <- as.data.frame(C7_GSEA)
gsea_csv_filename <- paste0(prefix,"_C7_GSEAResults.csv")
write.csv(gse_df_C7, file = gsea_csv_filename)

# CD4 LP vs L -------------------------------------------------------------

DE_table = as.data.frame(read.csv('CD4DEG/LP_vs_L.csv',row.names = 2))
prefix = 'CD4LPvsL'
DE_table = subset(DE_table, q_value <= 0.05)
gene_list = DE_table['normalized_effect']
gene_list <- gene_list[order(gene_list$normalized_effect, decreasing = TRUE), , drop = FALSE]
rowname = row.names(gene_list)
gene_list = as.numeric(unlist(gene_list))
names(gene_list) = rowname
gene_list = gene_list[startsWith(names(gene_list),'IG')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPS')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPL')==0]

# C7_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g_C7)
# gse_df_C7 <- as.data.frame(C7_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C7_GSEAResults.csv")
# write.csv(gse_df_C7, file = gsea_csv_filename)
# 
# 
# C5_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g_C5)
# gse_df_C5 <- as.data.frame(C5_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C5_GSEAResults.csv")
# write.csv(gse_df_C5, file = gsea_csv_filename)
gene_list = gene_list[abs(gene_list)>=1]
uplist = names(gene_list[gene_list>0])
downlist = names(gene_list[gene_list<0])
UPenriched <- enrichr(uplist, dbs)
DOWNenriched = enrichr(downlist, dbs)
ALLenriched = enrichr(names(gene_list), dbs)

for (DB_name in dbs) {
  sig_terms = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_UP_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
  
  sig_terms = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_DOWN_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
  
  sig_terms = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_ALL_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
}


gse <- gseGO(geneList = gene_list,
             ont = "ALL",
             keyType = "SYMBOL", #"ENSEMBL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,pAdjustMethod = "BH")
gse_df <- as.data.frame(gse)
gsea_csv_filename <-  paste0(prefix,"_GSEGO_GSEAResults.csv")
write.csv(gse_df, file = gsea_csv_filename)

# CD4 LP vs L -------------------------------------------------------------

DE_table = as.data.frame(read.csv('CD4DEG/IEL_vs_LP.csv',row.names = 2))
prefix = 'CD4IELvsLP'
DE_table = subset(DE_table, q_value <= 0.05)
gene_list = DE_table['normalized_effect']
gene_list <- gene_list[order(gene_list$normalized_effect, decreasing = TRUE), , drop = FALSE]
rowname = row.names(gene_list)
gene_list = as.numeric(unlist(gene_list))
names(gene_list) = rowname
gene_list = gene_list[startsWith(names(gene_list),'IG')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPS')==0]
gene_list = gene_list[startsWith(names(gene_list),'RPL')==0]

gene_list = gene_list[abs(gene_list)>0.25]
uplist = names(gene_list[gene_list>0])
downlist = names(gene_list[gene_list<0])
UPenriched <- enrichr(uplist, dbs)
DOWNenriched = enrichr(downlist, dbs)
ALLenriched = enrichr(names(gene_list), dbs)

for (DB_name in dbs) {
  sig_terms = UPenriched[[DB_name]][UPenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_UP_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
  
  sig_terms = DOWNenriched[[DB_name]][DOWNenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_DOWN_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
  
  sig_terms = ALLenriched[[DB_name]][ALLenriched[[DB_name]]$Adjusted.P.value<0.05,]
  filename <- paste0(prefix,"_",DB_name,"_ALL_enrichr_significantONLY.csv")
  write.csv(sig_terms, filename)
}

# C7_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g_C7)
# gse_df_C7 <- as.data.frame(C7_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C7_GSEAResults.csv")
# write.csv(gse_df_C7, file = gsea_csv_filename)
# 
# 
# C5_GSEA <- GSEA(gene_list, minGSSize = 5, pAdjustMethod = "BH", pvalueCutoff = 0.05, TERM2GENE = m_t2g_C5)
# gse_df_C5 <- as.data.frame(C5_GSEA)
# gsea_csv_filename <- paste0(prefix,"_C5_GSEAResults.csv")
# write.csv(gse_df_C5, file = gsea_csv_filename)

gse <- gseGO(geneList = gene_list,
             ont = "ALL",
             keyType = "SYMBOL", #"ENSEMBL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = organism,pAdjustMethod = "BH")
gse_df <- as.data.frame(gse)
gsea_csv_filename <-  paste0(prefix,"_GSEGO_GSEAResults.csv")
write.csv(gse_df, file = gsea_csv_filename)



# Select Terms ------------------------------------------------------------
# Define helper function to evaluate the "Overlap" column
# Define helper function to evaluate the "Overlap" column

enrichfiles = list.files(pattern = "Reactome.*\\.csv$", recursive = TRUE)
interested_terms = c("T cell","lymphocyte","leukocyte","antigen","chemo","CD","receptor","cytokine","effector",
                     "nterleukin","MHC","IL ","Th1","Th2","stress","Stress","signaling","Signaling","Inflam")
unwanted_terms = c("SARS","B cell","immunoglobulin","Cancer","cancer","Treg","Family","Rho")

evaldiv <- function(x) {
  eval(parse(text = x))
}

for (enrichfile in enrichfiles) {
  
  try({
    enrichdf = read.csv(enrichfile)
    # Generate output filename
    output_filename = gsub("_enrichr_significantONLY", "", enrichfile)
    output_filename = gsub("_2023_Human", "", output_filename)
    output_filename <- gsub(".csv", "_barplot.jpg", output_filename)
    
    pattern <- paste(interested_terms, collapse = "|")
    unwanted_pattern <- paste(unwanted_terms, collapse = "|")
    
    # Filter the dataframe to keep rows where the 'Description' column contains any of the interested terms
    filtered_gsea_df <- enrichdf[grepl(pattern, enrichdf$Term, ignore.case = TRUE), ]
    filtered_gsea_df <- filtered_gsea_df[!grepl(unwanted_pattern, filtered_gsea_df$Term, ignore.case = TRUE), ]
    filtered_gsea_df = filtered_gsea_df[order(filtered_gsea_df$Combined.Score,decreasing = TRUE), ]
    
    if (nrow(filtered_gsea_df) < 2) {
      filtered_gsea_df = enrichdf
    }
    
    # Clean up the 'Term' column
    filtered_gsea_df$Term <- gsub("R-HSA-\\d+", "", filtered_gsea_df$Term)
    filtered_gsea_df$Term <- trimws(filtered_gsea_df$Term)
    filtered_gsea_df$Overlap_num = sapply(filtered_gsea_df$Overlap, evaldiv)
    filtered_gsea_df$Term <- gsub("(.{30})\\s", "\\1\n", filtered_gsea_df$Term)
    
    if (nrow(filtered_gsea_df) > 5) {
      filtered_gsea_df = filtered_gsea_df[1:5, ]
    }
    
    # Set color and axis position based on whether the output file contains "UP"
    if (grepl("UP", output_filename)) {
      filtered_gsea_df$color <- 'red'
      x_position <- "bottom"
      y_position <- "left"
    } else {
      filtered_gsea_df$color <- 'blue'
      filtered_gsea_df$Combined.Score <- filtered_gsea_df$Combined.Score * -1
      
      x_position <- "top"
      y_position <- "right"
    }
    
    filtered_gsea_df$log.padj = round(-log10(filtered_gsea_df$Adjusted.P.value),3)
    
    barplot_enrichr = ggplot(filtered_gsea_df, aes(x = Combined.Score, y = reorder(Term, Combined.Score))) +
      geom_bar(aes(fill = color), stat = "identity") +  # Use 'fill' based on the color column
      # geom_text(aes(label = log.padj),               # Add overlap number labels # Position label in the middle of the bar
      #           hjust = -0.3,                           # Adjust horizontal alignment of labels
      #           size = 3, fontface  = "bold") +                             # Set label font size
      scale_fill_identity() +  # Use the exact colors from the 'color' column without mapping them
      theme_minimal() +
      labs(
        title = paste0("Bar Plot for GSEA Results - ", gsub("_", " ", gsub("_dotplot.jpg", "", output_filename))),
        x = "Enrichment",
        y = "Terms"
      ) +
      scale_x_continuous(position = x_position, labels = function(x) abs(x)) +  # Move x-axis to the top if necessary
      scale_y_discrete(position = y_position) +  # Move y-axis to the right if necessary
      theme(
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # Adjust x-axis text font size and rotation
        axis.text.y = element_text(size = 12, face = "bold"),         # Adjust y-axis text font size
        plot.title = element_text(size = 16, face = "bold"),          # Adjust plot title font size and bold
        legend.text = element_text(size = 10),                        # Adjust legend text font size
        legend.title = element_text(size = 12)                        # Adjust legend title font size
      )
    
    # Save the plot
    ggsave(output_filename, plot = barplot_enrichr, width = 9, height = 5*nrow(filtered_gsea_df)/12+1)
    
  }, silent = TRUE)
}

