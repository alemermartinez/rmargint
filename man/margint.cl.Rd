\name{margint.cl}
\alias{margint.cl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classic marginal integration procedures for additive models
%%  ~~function to do ... ~~
}
\description{Standard marginal integration procedures for additive models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
margint.cl(Xp, yp, point=NULL, windows, epsilon=1e-06, prob=NULL, type='0', degree=NULL, qderivate=FALSE, orderkernel=2, Qmeasure=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xp}{Matrix of explanatory variables (n by p).
%%     ~~Describe \code{Xp} here~~
}
  \item{yp}{Vector of responses (missing values are allowed).
%%     ~~Describe \code{yp} here~~
}
  \item{point}{Matrix of points where predictions will be computed and returned.
%%     ~~Describe \code{point} here~~
}
  \item{windows}{Vector or a squared matrix of bandwidths for the smoothing
estimation procedure.
%%     ~~Describe \code{windows} here~~
}
  \item{epsilon}{Convergence criterion.
%%     ~~Describe \code{epsilon} here~~
}
  \item{prob}{Probabilities of observing each response (n). Defaults to ``NULL''.
%%     ~~Describe \code{prob} here~~
}
  \item{type}{Three different type of estimators can be selected: type '0' (local constant
on all the covariates), type '1' (local linear smoother on all the covariates), type 'alpha'
(local polynomial smoother only on the direction of interest).
}
  \item{degree}{Degree of the local polynomial smoother in the direction of interest
when using the estimator of type 'alpha'. Defaults to ``NULL'' for the case when using
estimators of type '0' or '1'.
%%     ~~Describe \code{degree} here~~
}
  \item{orderkernel}{Order of the kernel used in the nuisance directions when
using the estimator of type 'alpha'. Defaults to ``2''.
}
  \item{qderivate}{If TRUE, it calculates g^(q+1)/(q+1)! for each component only
for the type 'alpha' method. Defaults to ``FALSE''.
}
  \item{Qmeasure}{A matrix of points where the integration procedure ocurrs. Defaults
to ``NULL'' for calcuting the integrals over the sample.
}
}
\details{
Three types of classical marginal integration procedures for additive models, that is, considering
a squared loss function.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{mu }{Estimate for the intercept.}
\item{g.matrix }{Matrix of estimated additive components (n by p).}
\item{prediction }{Matrix of estimated additive components
for the points listed in the argument point.}
\item{mul }{A vector of size p showing in each component the estimated intercept
that considers only that direction of interest when using the type 'alpha' method.}
\item{g.derivative }{Matrix of estimated derivatives of the additive components (only
when qderivate is ``TRUE'') (n by p).}
\item{prediction.derivate }{Matrix of estimated derivatives of the additive components
for the points listed in the argument point (only when qderivate is ``TRUE'').}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
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
