% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testingMisCovGeneral.R
\name{testingMisCovGeneral}
\alias{testingMisCovGeneral}
\title{testingMisCovGeneral}
\usage{
testingMisCovGeneral(trainRes, testData)
}
\arguments{
\item{trainRes}{Results results from trainging results}

\item{testData}{testing data set}
}
\value{
YHat the empirical risk of the test
}
\description{
Compute the empirical risk of cv or testing for kernel machies with non-monotone type of missing covariates
}
