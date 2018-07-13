\name{kernel4}
\alias{kernel4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Order 4 kernel
}
\description{A kernel of order 4.
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
0 if abs(x) > 1 and ( 15/32 ) * ( 1 - x^2 ) * ( 3 - 7 * x^2 ) otherwise.
}
\details{
A kernel L is a kernel of order 4 if it integrates 1, 
the integrals of u^j L(u) is 0 for 1 <= j < 4 (j integer) and
the integral of u^4 L(u) is different from 0.
}
\author{Alejandra Martinez, Matias Salibian-Barrera.
}