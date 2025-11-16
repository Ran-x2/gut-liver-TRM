library(circlize)
setwd('C:/Users/andre/Documents/GitHub/gut-liver-TRM/Revision/refined_clone_chord/')
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
mat = read.csv('donor_clone_normalized_mean_TCRab CD4.csv', row.names = 1, check.names = FALSE)
mat <- as.matrix(mat)

rownames(mat) <- sub("TCRab ", "", rownames(mat))
colnames(mat) <- sub("TCRab ", "", colnames(mat))
rownames(mat) <- sub("CD4 ", "", rownames(mat))
colnames(mat) <- sub("CD4 ", "", colnames(mat))

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

# chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = grid.col, order = sectors, transparency = 0.3,keep.diagonal =FALSE, symmetric = FALSE)

stopifnot(identical(rownames(mat), colnames(mat)))

ij <- which(upper.tri(mat), arr.ind = TRUE)
df <- data.frame(
  from  = rownames(mat)[ij[,1]],
  to    = colnames(mat)[ij[,2]],
  w_from= mat[ij],                          # i -> j
  w_to  = mat[cbind(ij[,2], ij[,1])]        # j -> i
)

# drop pairs that are all zero (optional)
df <- df[df$w_from > 0 | df$w_to > 0, ]

circos.clear()
chordDiagram(
  df,
  # your styling
  annotationTrack = "grid",
  preAllocateTracks = 1,
  grid.col = grid.col,
  order = sectors,
  transparency = 0.3,
  link.overlap = FALSE,
  directional = 0                 # one ribbon; set to 2 to add arrows on both ends
)

circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)
dev.copy(jpeg,'cd4_chord_donor_averaged_total_clone_normalized.png', width=14, height=14, units="in", res=500)
dev.off()
## ===========================================================
# CD8  chord diagram
mat = read.csv('donor_clone_normalized_mean_TCRab CD8ab.csv', row.names = 1, check.names = FALSE)
mat <- as.matrix(mat)

rownames(mat) <- sub("TCRab ", "", rownames(mat))
colnames(mat) <- sub("TCRab ", "", colnames(mat))
rownames(mat) <- sub("CD8ab ", "", rownames(mat))
colnames(mat) <- sub("CD8ab ", "", colnames(mat))

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
# chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = grid.col, order = sectors, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)

stopifnot(identical(rownames(mat), colnames(mat)))

ij <- which(upper.tri(mat), arr.ind = TRUE)
df <- data.frame(
  from  = rownames(mat)[ij[,1]],
  to    = colnames(mat)[ij[,2]],
  w_from= mat[ij],                          # i -> j
  w_to  = mat[cbind(ij[,2], ij[,1])]        # j -> i
)

# drop pairs that are all zero (optional)
df <- df[df$w_from > 0 | df$w_to > 0, ]

circos.clear()
chordDiagram(
  df,
  # below is my custom styling
  annotationTrack = "grid",
  preAllocateTracks = 1,
  grid.col = grid.col,
  order = sectors,
  transparency = 0.3,
  link.overlap = FALSE,
  directional = 0                 # one ribbon; set to 2 to add arrows on both ends
)


circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)
dev.copy(jpeg, 'cd8_chord_donor_averaged_total_clone_normalized.png', width=14, height=14, units="in", res=500)
dev.off()
