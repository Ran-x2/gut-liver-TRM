setwd('E:/AAA_Labwork/T cells/v2')
PBTeff_LTeff_CD8 = readRDS('PBTeff_LTeff_CD8.rds')
start <- slingshot(PBTeff_LTeff_CD8, clusterLabels = 'tissue.celltype', reducedDim = 'PCA')
metadata = colData(PBTeff_LTeff_CD8)

my_colors = c('PB TCRab CD8ab Teff' = '#A7E1F1','L TCRab CD8ab Teff' = '#7ABEE2','L TCRab CD8ab TRM' = '#1452a3')

#Set the pseudotime variable
t <- start$slingPseudotime_1

#Extract the gene expression matrix
Y <- assay(start)

# fit a GAM with a loess term for pseudotime
#Note: This takes a while
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})


#Select the top 50 most significant genes that change over pseudotime
gam.pval_sig = gam.pval[gam.pval<0.05]
topgenes <- c(names(sort(gam.pval_sig, decreasing = FALSE)))
topgenes[startsWith(topgenes,'RP') == TRUE] = ""
topgenes[startsWith(topgenes,'MT-') == TRUE] = ""
topgenes = topgenes[topgenes != ""]
topgenes = topgenes[1:50]
heatdata <- assay(start)[rownames(assay(start)) %in% topgenes, 
                         order(t, na.last = NA, decreasing = FALSE)]
#Scale the data per gene for visualization
heatdata <- t(scale(t(heatdata)))

#Trimm z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -2] = -2

# Create a dataframe with cell types in order of heatmap data
tissue.celltype_df <- data.frame(tissue.celltype = metadata$tissue.celltype[colnames(heatdata)])

# Ensure it's a factor with the correct levels
tissue.celltype_df$tissue.celltype <- factor(tissue.celltype_df$tissue.celltype, levels = names(my_colors))
tissue.celltype_df$pseudotime <- sort(t)
pseudotime_colors <- colorRampPalette(c("blue", "white", "red"))
breaks <- seq(from = min(tissue.celltype_df$pseudotime, na.rm = TRUE), 
              to = max(tissue.celltype_df$pseudotime, na.rm = TRUE), 
              length.out = 100)

# Define the colors for the cell types
annotation_colors = list(
  tissue.celltype = my_colors,
  pseudotime = pseudotime_colors(length(breaks) - 1)
)
# Include viridisLite package
library(viridisLite)

# Define color gradient with viridis function
color_gradient <- viridis(100)

# Generate heatmap with new color gradient
ptime_heatmap = pheatmap(heatdata,
                         annotation_col = tissue.celltype_df, 
                         annotation_colors = annotation_colors, 
                         cluster_rows = TRUE, 
                         cluster_cols = FALSE, 
                         show_rownames = TRUE, 
                         show_colnames = FALSE,
                         border_color = NA,
                         color = color_gradient,  # Apply the new color gradient
                         treeheight_row = 10,
                         width = 4, height = 8, fontsize = 8,fontsize_col = 12) 

ggsave("PBTeff_LTeff_CD8_CD8_time.pdf", plot = ptime_heatmap, width = 5, height = 7,dpi = 300)
