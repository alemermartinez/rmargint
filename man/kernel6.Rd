\name{kernel6}
\alias{kernel6}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Order 6 kernel
}
\description{A kernel of order 6.
}
\usage{
kernel2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A real number.
%%     ~~Describe \code{x} here~~
}
}
\value{
0 if abs(x) > 1 and ( 105/256 ) * ( 1 - x^2 ) * ( 5 - 30 * x^2 + 33 * x^4 ) otherwise.
}
\details{
A kernel L is a kernel of order 6 if it integrates 1, 
the integrals of u^j L(u) is 0 for 1 <= j < 6 (j integer) and
the integral of u^6 L(u) is different from 0.
}
\author{Alejandra Martinez, Matias Salibian-Barrera.
}