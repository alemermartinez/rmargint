\name{predict.margint}
\alias{predict.margint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitted values for objects of class \code{margint}
%%  ~~function to do ... ~~
}
\description{This function returns the fitted values given the covariates of the original sample under an additive model using a classical or robust marginal integration procedure estimator computed with \code{\link{margint.cl}} or \code{\link{margint.rob}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
##S3 method for class margint

predict(object, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
%%     ~~Describe \code{object} here~~
}
\item{...}{additional other arguments.}
}

\details{
This function returns the predicted values given the covariates of the original sample under an additive model using one of the three possible the classical or robust marginal integration procedures.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of fitted values.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Alejandra Martinez \email{alemartinez@unlu.edu.ar}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Additive Models }
