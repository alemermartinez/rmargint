\name{margint.rob}
\alias{margint.rob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust marginal integration procedures for additive models
%%  ~~function to do ... ~~
}
\description{Robust marginal integration procedures for additive models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
margint.rob(Xp, yp, point=NULL, windows, prob=NULL, sigma.hat=NULL, win.sigma=NULL, epsilon=1e-06, type='0', degree=NULL, typePhi='Huber', k.h=1.345, k.t = 4.685, max.it=20, qderivate=FALSE, orderkernel=2, Qmeasure=NULL)
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
  \item{prob}{Probabilities of observing each response (n). Defaults to ``NULL''.
%%     ~~Describe \code{prob} here~~
}
 \item{sigma.hat}{Estimate of the residual standard error. If NULL we use the mad of the residuals obtained with local medians.
%%     ~~Describe \code{sigma.hat} here~~
}
  \item{win.sigma}{
}
  \item{epsilon}{Convergence criterion.
%%     ~~Describe \code{epsilon} here~~
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
  \item{typePhi}{One of either ``Tukey'' or ``Huber''.
%%     ~~Describe \code{type} here~~
}
  \item{k.h}{Tuning constant for a Huber-type loss function. Defaults to ``1.345''.%%     ~~Describe \code{k.h} here~~
}
  \item{k.t}{Tuning constant for a Tukey-type loss function. Defaults to ``4.685''.
%%     ~~Describe \code{k.t} here~~
}
  \item{max.it}{Maximum number of iterations for the algorithm.
%%     ~~Describe \code{max.it} here~~
}
  \item{qderivate}{If TRUE, it calculates g^(q+1)/(q+1)! for each component only
for the type 'alpha' method. Defaults to ``FALSE''.
}
  \item{orderkernel}{Order of the kernel used in the nuisance directions when
using the estimator of type 'alpha'. Defaults to ``2''.
}
  \item{Qmeasure}{A matrix of points where the integration procedure ocurrs. Defaults
to ``NULL'' for calcuting the integrals over the sample.
}
}
\details{
Three types of robust marginal integration procedures for additive models.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{mu }{Estimate for the intercept.}
\item{g.matrix }{Matrix of estimated additive components (n by p).}
\item{sigma.hat }{Estimate of the residual standard error.}
\item{prediction }{Matrix of estimated additive components
for the points listed in the argument point.}
\item{mul }{A vector of size q showing in each component the estimated intercept
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
