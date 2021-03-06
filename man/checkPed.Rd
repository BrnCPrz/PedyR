% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{checkPed}
\alias{checkPed}
\title{Checks for errors and conformities in a pedigree file}
\usage{
checkPed(pedigree, id = 1, sire = 2, dam = 3, rmsingle = FALSE,
  verbose = FALSE)
}
\arguments{
\item{pedigree}{a data.frame with three columns: id, sire and dam}

\item{id}{integer containning the column where animal ids are located in the pedigree}

\item{sire}{integer containning the column where sire ids are located in the pedigree}

\item{dam}{integer containning the column where dam ids are located in the pedigree}
}
\description{
This function checks for common errors in genealogical structured data
}
\examples{
id = c(1,2,3,4,5,6,7,8,9,10)
sire = c(0,0,1,1,3,5,5,5,7,7)
dam = c(0,0,2,2,2,4,4,6,6,8)
ped = as.data.frame(cbind(id,sire,dam))
checkPed(ped)

}

