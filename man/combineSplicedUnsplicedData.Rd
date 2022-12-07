\name{combineSplicedUnsplicedData}
\alias{combineSplicedUnsplicedData}
\title{
UMI data}
\description{
Combines and log tranforms spliced and unspliced cell data
}
\usage{
combineSplicedUnsplicedData(spliced_data = spliced, unspliced_data = unspliced, log_transform = TRUE)
}
\arguments{
  \item{spliced_data}{a dgTMatrix object, default is 'spliced'.}
  \item{unspliced_data}{a dgTMatrix object, default is, 'unspliced'.}
  \item{log_transform}{logical, if true data is log transformed, default is, 'TRUE'.}
}
\details{UMI count of each barcode of both spliced and unspliced data is calcuated and log transformed}
\value{
returns dataframe of each barcode and the log transformed UMI value}
\examples{
# Mocking up some data:
spliced <- matrix(rpois(100,10),ncol=10)
unspliced <- matrix(rpois(100,10),ncol=10)

combineSplicedUnsplicedData()
}
