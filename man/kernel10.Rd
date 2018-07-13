\name{kernel10}
\alias{kernel10}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Order 10 kernel
}
\description{A kernel of order 10.
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
0 if abs(x) > 1 and 0.75 * ( 1 - x^2 ) * ( 315/128 - 105/32 * x^2 + 63/64 * x^4 - 3/32 * x^6 - 1/384 * x^8 ) otherwise.
}
\details{
A kernel L is a kernel of order 10 if it integrates 1, 
the integrals of u^j L(u) is 0 for 1 <= j < 10 (j integer) and
the integral of u^10 L(u) is different from 0.
}
\author{Alejandra Martinez, Matias Salibian-Barrera.
}