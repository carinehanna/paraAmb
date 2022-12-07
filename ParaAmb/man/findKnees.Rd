\name{findKnees}
\alias{findKnees}
\title{Find all knees}
\description{
Finds all knees that are remaining below the first knee.
}
\usage{
findKnees(data, first_knee)
}

\arguments{
  \item{data}{a dgTMatrix object.}
  \item{first_knee}{first knee barcodeRanks dataframe found from \code{firstKneeValue} function}
}
\details{
Finds all knee values and returns the UMI value of all knees in a dataframe.
}
\value{
A dataframe output containing all knee UMI values that are found
}

\examples{
# Mocking up some data:
set.seed(2000)
my.counts <- DropletUtils:::simCounts()
first_knee <- firstKneeValue(my.counts, 500)

findKnees(my.counts, first_knee)

}
