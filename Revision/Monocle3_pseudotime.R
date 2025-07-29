library(SingleCellExperiment)
library(monocle3)
library(scater)
library(scran)
library(pheatmap)
library(viridisLite)
library(mgcv)

setwd('E:/AAA_Labwork/T cells/v4_revision')
LPTEMTRM_LTRM_IELTRM = readRDS('LPTRM_IELTRM_CD4_counts.rds')
expr_matrix=assay(LPTEMTRM_LTRM_IELTRM)  # if this is your raw count data
my_colors =  c('L TCRab CD4 TRM' = '#68962D','IEL TCRab CD4 TRM' = '#3a8433','LP TCRab CD4 TRM' = '#94d53f')
gene_metadata=data.frame(
  gene_short_name = rownames(expr_matrix),
  row.names = rownames(expr_matrix)
)
cell_metadata=as.data.frame(colData(LPTEMTRM_LTRM_IELTRM))

cds=new_cell_data_set(
  expression_data = expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)
reducedDims(cds)$PCA=reducedDim(LPTEMTRM_LTRM_IELTRM, "PCA")
colData(cds)$partition=factor("partition1")  # Monocle still expects a 'partition'
colData(cds)$cluster=factor(colData(cds)$tissue.celltype)

cell_names=colnames(cds)
partitions_named=setNames(rep("1", length(cell_names)), cell_names)
clusters_named=setNames(as.character(colData(cds)$tissue.celltype), cell_names)

cds@clusters$UMAP=list(
  clusters = clusters_named,
  partitions = partitions_named,
  cluster_result = NA,
  partition_result = NA
)
reducedDims(cds)$UMAP=reducedDim(LPTEMTRM_LTRM_IELTRM, "UMAP")
cds=learn_graph(cds, use_partition = FALSE)

# Set root cells by label
# Option 1: manually specify root cell
# root_cells=colnames(cds)[colData(cds)$tissue.celltype == "LP TCRab CD4 TRM"]
# cds=order_cells(cds, root_cells = root_cells)

# Option 2: interactive (or choose root principal node)
cds=order_cells(cds)

pseudotime=pseudotime(cds)
pseudotime_df=data.frame(cell = names(pseudotime), pseudotime = pseudotime)
write.csv(pseudotime_df, "LPTRM_IELTRM_CD4_monocle3_pseudotime.csv")
pseudotime = read.csv("LPTRM_IELTRM_CD4_monocle3_pseudotime.csv")$pseudotime
# Get expression matrix
#expr_mat=assay(cds, "counts")
# or normalized log (after Seurat normalization)
expr_mat=log1p(cds@assays@data$counts)

gam.pval=apply(expr_mat, 1, function(gene_expr) {
  d=data.frame(z = gene_expr, t = pseudotime)
  tmp <- gam(z ~ lo(t), data=d)
  tmp
  #p <- summary(tmp)[4][[1]][1,5]
  #p
})
saveRDS(gam.pval, 'gam_models_LPTRM_IELTRM_CD4_counts.rds')

pvals <- sapply(gam.pval, function(model) {
  tryCatch({
    summary(model)$p.pv["lo(t)"]
  }, error = function(e) NA)
})

effects <- sapply(gam.pval, function(model) {
  tryCatch({
    summary(model)[[1]]["lo(t)"]
    #summary(model)[[1]][2]
  }, error = function(e) NA)
})
adj.pvals <- p.adjust(pvals, method = "BH")
sig_genes <- names(adj.pvals)[adj.pvals < 0.05]

topgenes <- names(sort(abs(effects[sig_genes]), decreasing = TRUE))
topgenes=topgenes[!grepl("^RP|^MT-", topgenes)]  # remove RPs and MTs
topgenes=topgenes[!is.na(topgenes)]
topgenes = gsub("\\.lo\\(t\\)", "", topgenes)[1:50]
topgenes
# Subset and scale data
heatdata=expr_mat[rownames(expr_mat) %in% topgenes, order(pseudotime, na.last = NA)]
heatdata=t(scale(t(as.matrix(heatdata))))
heatdata[heatdata > 3] = 3
heatdata[heatdata < -2] = -2

# Annotation
annotation_col=data.frame(
  tissue.celltype = colData(cds)[colnames(heatdata),]$tissue.celltype,
  pseudotime = pseudotime[colnames(heatdata)]
)
annotation_col$tissue.celltype=factor(annotation_col$tissue.celltype,
                                         levels = names(my_colors))

# Heatmap color settings
pseudotime_colors=colorRampPalette(c("blue", "white", "red"))
breaks=seq(min(annotation_col$pseudotime, na.rm = TRUE),
              max(annotation_col$pseudotime, na.rm = TRUE),
              length.out = 100)

annotation_colors = list(
  tissue.celltype = my_colors,
  pseudotime = pseudotime_colors(length(breaks) - 1)
)

# Heatmap
color_gradient=viridis(100)
saveRDS(heatdata, 'LPTRM_IELTRM_CD4_monocle3_heatdata.rds')
write.table(heatdata, 'LPTRM_IELTRM_CD4_monocle3_ptime_gam.csv', sep = ',')

ptime_heatmap = pheatmap(heatdata,
                         annotation_col = annotation_col,
                         annotation_colors = annotation_colors,
                         cluster_rows = TRUE,
                         cluster_cols = FALSE,
                         show_rownames = TRUE,
                         show_colnames = FALSE,
                         border_color = NA,
                         color = color_gradient,
                         treeheight_row = 10,
                         width = 4, height = 8, fontsize = 8, fontsize_col = 12)

ggsave("LPTRM_IELTRM_CD4_monocle3_time.pdf", plot = ptime_heatmap, width = 5, height = 6, dpi = 300)
