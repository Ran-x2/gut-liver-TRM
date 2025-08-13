library(circlize)
setwd("G:/My Drive/result/publication/cellreport/revision/GLIPH2/GLIPH2_results")
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
## ===========================================================
# CD4 GLIPH2 chord diagram
mat = read.csv('CD4_GLIPH2_migration_counts.csv', row.names = 1, check.names = FALSE)
mat <- as.matrix(mat)

rownames(mat) <- sub("TCRab ", "", rownames(mat))
colnames(mat) <- sub("TCRab ", "", colnames(mat))
rownames(mat) <- sub("CD4 ", "", rownames(mat))
colnames(mat) <- sub("CD4 ", "", colnames(mat))
mat

panel.fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.4)
  circos.axis(h = "top", labels.cex = 0.5,  sector.index = sector.name, track.index = 2)
}

# set global parameters
circos.clear()
circos.par(cell.padding = c(0,0,0,0), gap.degree = 2, start.degree = 60, points.overflow.warning = FALSE)
par(mar=c(5,6,4,1)+.1)
sectors <- sort(unique(c(rownames(mat), colnames(mat))))
grid.col <- setNames(celltype_palette_cd4[sectors], sectors)
chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = grid.col, order = sectors, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)
circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)
dev.copy(jpeg,'cd4_GLIPH2_chord.png', width=14, height=14, units="in", res=500)
dev.off()
## ===========================================================
# CD8 GLIPH2 chord diagram
mat = read.csv('CD8_GLIPH2_migration_counts.csv', row.names = 1, check.names = FALSE)
mat <- as.matrix(mat)

rownames(mat) <- sub("TCRab ", "", rownames(mat))
colnames(mat) <- sub("TCRab ", "", colnames(mat))
rownames(mat) <- sub("CD8ab ", "", rownames(mat))
colnames(mat) <- sub("CD8ab ", "", colnames(mat))
mat

panel.fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.4)
  circos.axis(h = "top", labels.cex = 0.5,  sector.index = sector.name, track.index = 2)
}

# set global parameters
circos.clear()
circos.par(cell.padding = c(0,0,0,0), gap.degree = 2, start.degree = 60, points.overflow.warning = FALSE)
par(mar=c(5,6,4,1)+.1)
sectors <- sort(unique(c(rownames(mat), colnames(mat))))
grid.col <- setNames(celltype_palette_cd8[sectors], sectors)
chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = grid.col, order = sectors, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)
circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)
dev.copy(jpeg,'cd8_GLIPH2_chord.png', width=14, height=14, units="in", res=500)
dev.off()