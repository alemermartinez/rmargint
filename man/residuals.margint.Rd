\name{residuals.margint}
\alias{residuals.margint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals of a marginal integration fit in an additive model
%%  ~~function to do ... ~~
}
\description{residuals method for class \code{margint}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
##S3 method for class margint

residuals(object, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.}
%%     ~~Describe \code{object} here~~
\item{...}{other arguments.}
}

\details{
This function returns the residuals of the fitted additive model using one of the three different types of classical or robust estimators based on the marginal integration procedure.
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
