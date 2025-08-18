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

# root_cells=colnames(cds)[colData(cds)$tissue.celltype == "L TCRab CD4 TRM"]
# cds=order_cells(cds, root_cells = root_cells)

cds=order_cells(cds)

pseudotime=pseudotime(cds,reduction_method = "UMAP")
pseudotime_df=data.frame(cell = names(pseudotime), pseudotime = pseudotime)
write.csv(pseudotime_df, "revision/monocle3_ptime/LPTRM_IELTRM_CD4_monocle3_pseudotime.csv")
# pseudotime = read.csv("LPTRM_IELTRM_CD4_monocle3_pseudotime.csv",row.names = 1)['pseudotime']

# Get expression matrix
#expr_mat=assay(cds, "counts")
# or normalized log (after Seurat normalization)
expr_mat=log1p(cds@assays@data$counts)

genes <- rownames(expr_mat)
ptime <- as.numeric(pseudotime)
k_basis <- 10

fit_one <- function(y, ptime, k = k_basis) {
  d <- data.frame(y = as.numeric(y), ptime = ptime)
  # Gaussian is fine for log-normalized expression
  mgcv::bam(y ~ s(ptime, k = k), data = d,
            method = "REML", discrete = TRUE, select = TRUE)
}
fits <- lapply(seq_len(nrow(expr_mat)), function(i) fit_one(expr_mat[i, ], ptime))
names(fits) <- genes
saveRDS(fits, 'revision/monocle3_ptime/mgcv_gam_models_LPTRM_IELTRM_CD4_counts.rds')

extract_stats <- function(m, ptime) {
  s <- summary(m)
  p <- if (!is.null(s$s.table)) as.numeric(s$s.table[1, "p-value"]) else NA_real_
  
  # deviance explained (model-level magnitude, 0..1)!
  dev <- if (!is.null(s$dev.expl)) as.numeric(s$dev.expl) else (1 - m$deviance / m$null.deviance)
  
  # amplitude of the fitted smooth on the response scale
  # evaluate at observed ptime to get the partial term contribution
  term_mat <- try(predict(m, type = "terms"), silent = TRUE)
  rng <- if (!inherits(term_mat, "try-error")) diff(range(term_mat[, 1], na.rm = TRUE)) else NA_real_
  
  c(p = p, dev_expl = dev, effect_range = rng, edf = s$edf[1])
}

stat_mat <- t(vapply(fits, extract_stats, numeric(4), ptime = ptime))
res <- as.data.frame(stat_mat)
res$padj <- p.adjust(res$p, method = "BH")


# pvals <- sapply(gam.pval, function(model) {
#   tryCatch({
#     # summary(model)$p.pv["lo(t)"]
#     summary(model)[4][[1]][1,5]
#   }, error = function(e) NA)
# })
# 
# effects <- sapply(gam.pval, function(model) {
#   tryCatch({
#     # summary(model)[[1]]["lo(t)"]
#     summary(model)[[1]][2]
#   }, error = function(e) NA)
# })
# saveRDS(gam.pval, 'revision/monocle3_ptime/gam_models_LPTRM_IELTRM_CD4_counts.rds')

sig_genes <- subset(res, is.finite(padj) & padj < 0.05)
topgenes   <- sig_genes[order(-sig_genes$dev_expl, sig_genes$padj), ]
topgenes=topgenes[!grepl("^RP|^MT-", row.names(topgenes)),]  # remove RPs and MTs
topgenes = row.names(topgenes)[1:50]
# topgenes = row.names(expr_matrix)[1:50] #testing color map!!
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
  # pseudotime = pseudotime[colnames(heatdata),]
)
rownames(annotation_col) = colnames(heatdata)
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

saveRDS(heatdata, 'revision/monocle3_ptime/LPTRM_IELTRM_CD4_monocle3_heatdata_mcgv_gam.rds')
write.table(heatdata, 'revision/monocle3_ptime/LPTRM_IELTRM_CD4_monocle3_ptime_mcgv_gam.csv', sep = ',')



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

ggsave("revision/monocle3_ptime/LPTRM_IELTRM_CD4_monocle3_time_mcgv_gam.pdf", plot = ptime_heatmap, width = 5, height = 6, dpi = 300)
