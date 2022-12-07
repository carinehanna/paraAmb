\name{filterUMIData}
\alias{filterUMIData}
\title{Filter data}
\description{Filter data based on filter number to optomize clustering}
\usage{
filterUMIData(ds_comb, filter_num = 1)
}
\arguments{
  \item{ds_comb}{dataframe containing spliced and unspliced UMI data}
  \item{filter_num}{a number, any UMI value below will be removed from the dataset. Default filter value is 1 if data has been log transformed. If raw data is used a filter number of 100 is recommended.}
}
\details{Raw data can cause clustering to not work properly due to "stripes". Filtering data optomized the kmeans clustering.}
\value{
Returns dataframe, containg the filtered log transformed UMI values of each barcode.
Returns a plot of the total spliced UMI values and total unspliced UMI values.}
\examples{
# Mocking up some data:
spliced <- matrix(rpois(100,10),ncol=10)
unspliced <- matrix(rpois(100,10),ncol=10)

total_spliced <- colSums(spliced)
total_unspliced <- colSums(unspliced)

comb <-as.data.frame(cbind(total_spliced, total_unspliced))

comb_data <- log10(comb +1)

filterUMIData(comb_data, 1)
}
