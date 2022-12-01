\name{compareKneeClusterPlot}
\alias{compareKneeClusterPlot}
\title{Compare Plots}
\description{
Compare cluster plot and knee group plot to find any correlation between them.
}
\usage{
compareKneeClusterPlot(filt_comb_data, knees, k)
}
\arguments{
  \item{filt_comb_data}{a dataframe containing total spliced and total unspliced data, can be created from \code{combineSplicedUnsplicedData}}
  \item{knees}{a dataframe containg knee UMI values, can be created from \code{findKnees}}
  \item{k}{a scatterplot containing each barcode and which knee group it belongs to, can be made from \code{kneeGroupsPlot} }
}
\value{
Returns two plots; a fviz_cluster plot and a scatter plot}
\examples{
# Mocking up some data:
set.seed(2000)
my.counts <- DropletUtils:::simCounts()

filtered_mat <- my.counts[30:100, 30:100]

cell_name <- NULL
for (i in 1:71)
{
  cell_name[[length(cell_name) + 1]] <- paste("cell", i)
}
colnames(filtered_mat) <- cell_name

umi_count <- colSums(filtered_mat)

knees <- as.data.frame(quantile(umi_count, probs = c(0.75, 0.5, 0.25)))
colnames(knees) <- c('knee_vals')


knee_groups <- filterKnees(filtered_mat, knees)

spliced <- matrix(rpois(5041,71),ncol=71)
unspliced <- matrix(rpois(5041,71),ncol=71)

comb_data <- combineSplicedUnsplicedData(spliced, unspliced)

k <- kneeGroupsPlot(comb_data, knee_groups)

compareKneeClusterPlot(comb_data, knees, k)

}
