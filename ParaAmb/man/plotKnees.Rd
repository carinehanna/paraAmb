\name{plotKnees}
\alias{plotKnees}
\title{
Plot all Knees
}
\description{
Used to ensure all knees are found properly, and visulaization of where the knees are located}
\usage{
plotKnees(first_knee, knee_data)
}
\arguments{
  \item{first_knee}{barcodeRanks dataframe of first knee}
  \item{knee_data}{dataframe of all knees and their UMI values}
}
\value{
Returns plot of all knees
}
\examples{
# Mocking up some data:
set.seed(2000)
my.counts <- DropletUtils:::simCounts()
first_knee <- DropletUtils::barcodeRanks(my.counts, lower=500, exclude.from = 50) #can edit

knees <- as.data.frame(quantile(colSums(my.counts), probs = c(.95, 0.5, 0.25)))
colnames(knees) <- c('knee_vals')

plotKnees(first_knee, knees)
}
