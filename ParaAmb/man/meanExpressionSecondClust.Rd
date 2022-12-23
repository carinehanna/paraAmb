\name{meanExpressionSecondClust}
\alias{meanExpressionSecondClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean Gene Expression of Second Cluster}
\description{
Caculates the mean gene expression of the second cluster}
\usage{
meanExpressionSecondClust(kmeans_data, cluster_data)
}
\arguments{
  \item{kmeans_data}{Dataframe created from \code{kmeans}}
  \item{cluster_data}{Dataframe containing barcode and cluster value they are found in. Can be created in \code{clusterGroups}.}
}
\details{
Filters spliced and unspliced matrices by barcodes found in the second cluster. Calculates the total UMI for each gene and returns a plot of the values.}
\value{
Returns plot of gene expression between spliced and unspliced.
Returns dataframe of gene count.}
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

kmeans_data1 <- kmeans(comb_data, 3)


knees <- as.data.frame(quantile(comb_data$colSums.filt_spliced., probs = c(0.75, 0.5, 0.25)))
colnames(knees) <- c('knee_vals')

cluster_data1 <- clusterGroups(comb_data, knees)
meanExpressionSecondClust(kmeans_data1, cluster_data1)
}
