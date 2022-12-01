\name{clusterGroups}
\alias{clusterGroups}
\title{Clustering of data}
\description{
Using \code{kmeans} to find clusters in data.
}
\usage{
clusterGroups(filter_data, sort_knee)
}
\arguments{
  \item{filter_data}{dataframe containing total spliced, total unspliced and barcode. Can be created from \code{filterSplicedData} or \code{combineSplicedUnsplicedData}}
  \item{sort_knee}{dataframe containg knee UMI values. Can be created from \code{findKnees}}
}
\details{
Uses \code{kmeans} and \code{fviz_cluster} to cluster the data sets. Clusters seperate barcodes containing ambient RNA from barcdoes containing actual cells}
\value{
Returns filter_data dataframe with added column "cluster", containing the cluster number the barcode is apart of}
\examples{
# Mocking up some data:
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

clusterGroups(comb_data, knees)
}
