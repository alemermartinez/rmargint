\name{plot.margint}
\alias{plot.margint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic plots for objects of class \code{margint}
%%  ~~function to do ... ~~
}
\description{Plot method for class \code{margint}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
##S3 method for class margint

plot(object, which=1:np, ask=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object}{an object of class \code{margint}, a result of a call to \code{\link{margint.cl}} or \code{\link{margint.rob}}.
%%     ~~Describe \code{object} here~~
}
\item{which}{vector of indices of explanatory variables for which partial residuals plots will be generated. Defaults to all available explanatory variables.}
\item{ask}{logical value. If \code{TRUE}, the graphical device will prompt before going to the next page/screen of output.}
\item{...}{additional other arguments.}
}

%\details{
%This function plots, for each additive component, the partial residuals and the estimated curve.
%}

%%\value{
%%\item{mu }{Estimate for the intercept.}}

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
