\name{kneeGroupsPlot}
\alias{kneeGroupsPlot}
\title{
Knee Group Plot
}
\description{Scatter plot of barcodes and which group they belong in based on found knees}
\usage{
kneeGroupsPlot(ds_comb, knee_groups)
}
\arguments{
  \item{ds_comb}{dataframe containing total spliced, total unspliced and barcodes}
  \item{knee_groups}{matrix containing the knee groups that the barcodes are located in. Can be created from \code{filterKnees}}
}
\details{Using the logical matrix, barcodes are assigned a colour based on the group they are found in.
}
\value{
Returns a scatter plot}
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

kneeGroupsPlot(comb_data, knee_groups)


}

