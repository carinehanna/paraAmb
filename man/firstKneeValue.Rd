\name{firstKneeValue}
\alias{firstKneeValue}
\title{Finding First Knee}
\description{
  Returns a barcodeRanks dataframe about the first knee that is found.
  Returns plot of cell data with location of first knee outlined.}
\usage{
firsKneeValue(data=spliced, lower_num=2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dgTMatrix object.}
  \item{lower_num}{a number to determine the lower argument of barcodeRanks.}
}
\details{
The returned plot displays the location of the first knee. It is used to visualize if the correct knee has been located. If the found knee is too high the lower_num can be reduced.
}
\value{
  \code{firstKneeValue} A barcodeRanks dataframe containing:

      "rank" : Numeric, the rank of each barcode (averaged across ties).

      "total" : Numeric, the total counts for each barcode.

      "fitted" : Numeric, the fitted value from the spline for each barcode. This is NA for points with x outside of fit.bounds.

  Also returns a plot displaying barcodeRanks first knee data
  }
\examples{
# Mocking up some data:
set.seed(2000)
my.counts <- DropletUtils:::simCounts()

firstknee <- firstKneeValue(my.counts, 500)
}
