\name{filterKnees}
\alias{filterKnees}
\title{
Grouping and Filtering of SCE data
}
\description{
Filters all barcodes to determine if the UMI count is higher than knee 1 UMI value, in between the knee 1 UMI value and the knee 2 UMI value, or in between the knee 2 UMI value and the knee 3 UMI value if a third knee exists.}
\usage{
filterKnees(data, knee_data)
}
\arguments{
  \item{data}{a dgTMatrix object.}
  \item{knee_data}{dataframe containing UMI values of all knees found from \code{findKnees} function}
}
\value{
Returns Large logic matrix with cols 'is.betweenknee1andknee2', 'is.betweenknee2andknee3', 'is.higherknee1'. Used to compare with kmeans cluster.}
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


filterKnees(filtered_mat, knees)
}
