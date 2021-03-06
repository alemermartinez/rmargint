% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmargint-fn.R
\name{predict.margint}
\alias{predict.margint}
\title{Fitted values for objects of class \code{margint}}
\usage{
\method{predict}{margint}(object, ...)
}
\arguments{
\item{object}{an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.}

\item{...}{additional other arguments. Currently ignored.}
}
\value{
A vector of fitted values.
}
\description{
This function returns the fitted values given the covariates of the original sample under an additive model using a classical or robust marginal integration procedure estimator computed with \code{\link{margint.cl}} or \code{\link{margint.rob}}.
}
\author{
Alejandra Mercedes Martinez \email{ale_m_martinez@hotmail.com}
}
