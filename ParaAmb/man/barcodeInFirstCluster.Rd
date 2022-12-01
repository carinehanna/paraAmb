\name{barcodeInFirstCluster}
\alias{barcodeInFirstCluster}
\title{
All Barcodes in First Cluster
}
\description{
Returns all barcodes that are found in the first cluster (actual cells)}
\usage{
barcodeInFirstCluster(kmeans_data, cluster_data)
}
\arguments{
  \item{kmeans_data}{Dataframe created from \code{kmeans}}
  \item{cluster_data}{Dataframe containing barcode and cluster value they are found in. Can be created in \code{clusterGroups}.}
}
\value{
Returns dataframe containing all barcodes found in first cluster.
}
\examples{
raw_spliced <- DropletUtils:::simCounts()
raw_unspliced <- DropletUtils:::simCounts()

spliced <- raw_spliced[30:89, 30:89]
unspliced <- raw_unspliced[30:89, 30:89]

cell_name <- NULL
for (i in 1:60)
{
  cell_name[[length(cell_name) + 1]] <- paste("cell", i)
}
colnames(spliced) <- cell_name
colnames(unspliced) <- cell_name

comb_data <- data.frame(colSums(filt_spliced), colSums(filt_unspliced))
comb_data <- log10(comb_data + 1)

kmeans_data <- kmeans(comb_data, 3)

knees <- as.data.frame(quantile(comb_data$colSums.filt_spliced., probs = c(0.75, 0.5, 0.25)))
colnames(knees) <- c('knee_vals')

cluster_data <- clusterGroups(comb_data, knees)

barcodeInFirstCluster(kmeans_data, cluster_data)
}
