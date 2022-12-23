# ------
# spliced <- Matrix::readMM('~/Downloads/Remapped_for_ParaAmb/P5TLH/velocyto/P5TLH_S_sparse_matrix.mtx')
# unspliced <- Matrix::readMM('~/Downloads/Remapped_for_ParaAmb/P5TLH/velocyto/P5TLH_U_sparse_matrix.mtx')
# barcodes_SU <- read.table('~/Downloads/Remapped_for_ParaAmb/P5TLH/velocyto/P5TLH_cellIDs.tsv')
# genes_SU <- read.table('~/Downloads/Remapped_for_ParaAmb/P5TLH/velocyto/P5TLH_geneNames.tsv')
# # fixing barcode and gene dataframes
# dash <- sub('x','-1',barcodes_SU$V1)
# new_bars <- strsplit(dash,split = "P5TLH_cellranger:")
# new_barcodes <- do.call(rbind.data.frame, new_bars)
#
# colnames(spliced) <- new_barcodes[,2]
# rownames(spliced) <- genes_SU[,1]
# colnames(unspliced) <- new_barcodes[,2]
# rownames(unspliced) <- genes_SU[,1]
#

# -------------
# SU <- Matrix::readMM('~/Downloads/Remapped_for_ParaAmb/Neuron_Ear_Dataset/Cells3/cellranger_introns_plus_exons/raw_feature_bc_matrix/matrix.mtx')
# S <- Matrix::readMM('~/Downloads/Remapped_for_ParaAmb/Neuron_Ear_Dataset/Cells3/cellranger_exons_only/raw_feature_bc_matrix/matrix.mtx')
# barcodes_SU <- read.table('~/Downloads/Remapped_for_ParaAmb/Neuron_Ear_Dataset/Cells3/cellranger_introns_plus_exons/raw_feature_bc_matrix/barcodes.tsv')
# genes_SU <- read.table('~/Downloads/Remapped_for_ParaAmb/Neuron_Ear_Dataset/Cells3/cellranger_introns_plus_exons/raw_feature_bc_matrix/features.tsv')
# -----
SU <- Matrix::readMM('~/Downloads/Remapped_for_ParaAmb/Denishenko_MouseKidney/Nuclei1/cellranger_introns_plus_exons/raw_feature_bc_matrix/matrix.mtx')
S <- Matrix::readMM('~/Downloads/Remapped_for_ParaAmb/Denishenko_MouseKidney/Nuclei1/cellranger_exons_only/raw_feature_bc_matrix/matrix.mtx')
barcodes_SU <- read.table('~/Downloads/Remapped_for_ParaAmb/Denishenko_MouseKidney/Nuclei1/cellranger_introns_plus_exons/raw_feature_bc_matrix/barcodes.tsv')
genes_SU <- read.table('~/Downloads/Remapped_for_ParaAmb/Denishenko_MouseKidney/Nuclei1/cellranger_introns_plus_exons/raw_feature_bc_matrix/features.tsv')

colnames(SU) <- barcodes_SU[,1]
rownames(SU) <- genes_SU[,1]
colnames(S) <- barcodes_SU[,1]
rownames(S) <- genes_SU[,1]
unspliced <- SU - S
spliced <- S
# -----
# creating example names
actual.cell.name <- NULL
for (i in 1:1500)
{
  actual.cell.name[[length(actual.cell.name) + 1]] <- paste("cell", i)
}
real.cell.name <- as.data.frame(matrix(unlist(actual.cell.name),nrow=length(actual.cell.name),byrow=TRUE))
colnames(real.cell.name)[1] <- c('barcodes')

amb.cell.name <- NULL
for (i in 1:1500)
{
  amb.cell.name[[length(amb.cell.name) + 1]] <- paste("amb", i)
}
ambi.cell.name <- as.data.frame(matrix(unlist(amb.cell.name),nrow=length(amb.cell.name),byrow=TRUE))
colnames(ambi.cell.name)[1] <- c('barcodes')

#----
fk <- firstKneeValue(spliced, 50)
lowerUMIs <- metadata(fk)$knee
UMI_tots <- as.data.frame(colSums(spliced))
colnames(UMI_tots)[1] = 'UMI'
orig.real.cellss <- UMI_tots %>% filter(UMI >= lowerUMIs)
orig.ambients <- UMI_tots %>% filter(UMI < lowerUMIs)

#----
# creating mean list
# mean.para.data <- NULL
# mean.ed.data <- NULL
# store data for current run
paraAmb.data <- NULL
emptyDrops.data <- NULL

for (j in 1:100) {
orig.real.cells <- cbind(cellID = rownames(orig.real.cellss), orig.real.cellss)

real.barcodes <- sample(orig.real.cells$cellID, replace = TRUE, 1500)

true.cells.spliced <- spliced[,real.barcodes]
true.cells.unspliced <- unspliced[,real.barcodes]

colnames(true.cells.spliced) <- real.cell.name$barcodes
colnames(true.cells.unspliced) <- real.cell.name$barcodes

# adding variation
down.cells.spliced <- downsampleMatrix(true.cells.spliced, prop = runif(1500, min = 0.5, max = 1), bycol = TRUE, sink = NULL)
down.cells.unspliced <- downsampleMatrix(true.cells.unspliced, prop = runif(1500, min = 0.5, max = 1), bycol = TRUE, sink = NULL)

# ambient
orig.ambient <- cbind(cellID = rownames(orig.ambients), orig.ambients)

ambient.barcodes <- sample(orig.ambient$cellID, replace = TRUE, 1500)

true.ambient.spliced <- spliced[,ambient.barcodes]
true.ambient.unspliced <- unspliced[,ambient.barcodes]

colnames(true.ambient.spliced) <- ambi.cell.name$barcodes
colnames(true.ambient.unspliced) <- ambi.cell.name$barcodes

# adding variation
down.amb.spliced <- downsampleMatrix(true.ambient.spliced, prop = runif(1500, min = 0.8, max = 1), bycol = TRUE, sink = NULL)
down.amb.unspliced <- downsampleMatrix(true.ambient.unspliced, prop = runif(1500, min = 0.25, max = 0.8), bycol = TRUE, sink = NULL)


# combining all data
all.sim.spliced <- cbind(down.cells.spliced, down.amb.spliced)
all.sim.unspliced <- cbind(down.cells.unspliced, down.amb.unspliced)

# testing emptyDrops
e.out <- emptyDrops(all.sim.spliced)
is.cell <- e.out$FDR <= 0.05
be <- as.data.frame(cbind(real.cell.name$barcodes, is.cell))
be.out <- be %>% filter(is.cell == TRUE) %>% select(V1)
colnames(be.out)[1] <- c('barcodes')
table(is.cell)

# testing paraamb
fk <- firstKneeValue(all.sim.spliced, 50)
knees <- findKnees(all.sim.spliced, fk)
plotKnees(fk, knees)
ds_comb <- combineSplicedUnsplicedData(all.sim.spliced, all.sim.unspliced)
filt_comb_data <- filterUMIData(ds_comb, 0)
cluster <- clusterGroups(filt_comb_data, knees)
bp <- barcodeInFirstCluster(cluster[1], cluster[2])


# PARA
missed.para <- anti_join(real.cell.name, bp, by = 'barcodes')

# ED
missed.ed <- anti_join(real.cell.name, be.out, by = 'barcodes')

count.missed.para <- as.numeric(nrow(missed.para))
count.missed.ed <- as.numeric(nrow(missed.ed))

total.cor.found.para <- length(actual.cell.name) - count.missed.para
total.cor.found.ed <- length(actual.cell.name) - count.missed.ed

percent.cor.para <- total.cor.found.para / length(actual.cell.name)
percent.cor.ed <- total.cor.found.ed / length(actual.cell.name)

paraAmb.data[[length(paraAmb.data) + 1]] <- percent.cor.para
emptyDrops.data[[length(emptyDrops.data) + 1]] <- percent.cor.ed

}


# calculate mean
final.para.data <- as.data.frame(matrix(unlist(paraAmb.data),nrow=length(paraAmb.data),byrow=TRUE))
final.ed.data <- as.data.frame(matrix(unlist(emptyDrops.data),nrow=length(emptyDrops.data),byrow=TRUE))


mean.para.data[[length(mean.para.data) + 1]] <- mean(final.para.data$V1)
mean.ed.data[[length(mean.ed.data) + 1]] <- mean(final.ed.data$V1)




