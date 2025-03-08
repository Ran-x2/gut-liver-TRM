library(circlize)
setwd('E:/AAA_Labwork/T cells/v2')
tissue_colors <- c("PB"="#b23429", "L"="#655e2f", "LP"="#db7843", "IEL"="#f4bf5c")
celltype_palette_cd4 <- c('PB TCM'= '#94d53f',
                          'IEL TRM'= '#3a8433',
                          'LP FOXP3+ Treg'= '#BEC765', 
                          'L TRM'= '#3a8433',
                          'IEL Mobile TRM' = '#519f38',
                          'PB FOXP3+ Treg'= '#BEC765',
                          'L FOXP3+ Treg'= '#BEC765',
                          'LP TRM'= '#3a8433',
                          'IEL FOXP3+ Treg'= '#BEC765',
                          'LP Tph'= '#6b866b',
                          'L TCM'= '#94d53f',
                          'LP Naive/TCM'= '#bff141',
                          'L Naive/TCM'= '#bff141',
                          'LP Poised TCM' = "#68b93c",
                          'LP Mobile TRM' = '#519f38',
                          'PB Naive/TCM'= '#bff141')
celltype_palette_cd8 <- c('IEL TRM'= '#132876',
                          'L MAIT'= '#09EED0',
                          'L Naive/TCM'= '#A7E1F1',
                          'L Teff'= '#4c9bd3', 
                          'L TRM'= '#132876', 
                          'LP TCM'= '#7ABEE2',
                          'LP TEM'= '#1452a3', 
                          'LP TRM'= '#132876', 
                          'PB Naive/TCM'= '#A7E1F1', 
                          'PB Teff'= '#4c9bd3',
                          'PB TEM'= '#1452a3')

# define panel function


mat1 <- readRDS('0_TCRab CD8ab_chord.rds')
mat2 <- readRDS('1_TCRab CD8ab_chord.rds')
mat3 <- readRDS('2_TCRab CD8ab_chord.rds')

panel.fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.4)
  circos.axis(h = "top", labels.cex = 0.5,  sector.index = sector.name, track.index = 2)
}

# set global parameters
circos.par(cell.padding = c(0,0,0,0), gap.degree = 2, start.degree = 90, points.overflow.warning = FALSE)
par(mar=c(5,6,4,1)+.1)

new_order <- sort(rownames(mat1))
mat1 <- mat1[new_order, new_order]
mat2 <- mat2[new_order, new_order]
mat3 <- mat3[new_order, new_order]

mat = mat1 + mat2 +mat3
rownames(mat) <- sub("TCRab ", "", rownames(mat))
colnames(mat) <- sub("TCRab ", "", colnames(mat))
rownames(mat) <- sub("CD4 ", "", rownames(mat))
colnames(mat) <- sub("CD4 ", "", colnames(mat))
rownames(mat) <- sub("CD8ab ", "", rownames(mat))
colnames(mat) <- sub("CD8ab ", "", colnames(mat))
chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = celltype_palette_cd8, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)
circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)


dev.copy(jpeg,'cd8_chord.png', width=14, height=14, units="in", res=500)
dev.off()

mat1 <- readRDS('0_TCRab CD4_chord.rds')
mat2 <- readRDS('1_TCRab CD4_chord.rds')
mat3 <- readRDS('2_TCRab CD4_chord.rds')

panel.fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2)
  circos.axis(h = "top", labels.cex = 0.5,  sector.index = sector.name, track.index = 2)
}

new_order <- sort(rownames(mat1))
mat1 <- mat1[new_order, new_order]
mat2 <- mat2[new_order, new_order]
mat3 <- mat3[new_order, new_order]

mat = mat1 + mat2 +mat3
rownames(mat) <- sub("TCRab ", "", rownames(mat))
colnames(mat) <- sub("TCRab ", "", colnames(mat))
rownames(mat) <- sub("CD4 ", "", rownames(mat))
colnames(mat) <- sub("CD4 ", "", colnames(mat))
rownames(mat) <- sub("CD8ab ", "", rownames(mat))
colnames(mat) <- sub("CD8ab ", "", colnames(mat))
chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = celltype_palette_cd4, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)
circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)


dev.copy(jpeg,'cd4_chord.png', width=14, height=14, units="in", res=500)
dev.off()
