\name{summary.margint}
\alias{summary.margint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing fits of additive models using a marginal integration procedure
%%  ~~function to do ... ~~
}
\description{summary method for class \code{margint}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
## S3 method for class margint.

summary(object, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
%%     ~~Describe \code{object} here~~
}
\item{...}{other arguments.}
}

\details{
This function returns the estimation of the intercept and also the five-number summary and the mean of the residuals for both classical and robust estimators. For the robust estimator it also returns the estimate of the residual standard error.
%%  ~~ If necessary, more details than the description above ~~
}
%%\value{
%%\item{mu }{Estimate for the intercept.}}

\references{
%% ~put references to the literature/web site here ~
}
\author{Alejandra Martinez, Matias Salibian-Barrera.
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
