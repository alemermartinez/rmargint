\name{kernel8}
\alias{kernel8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Order 8 kernel
}
\description{A kernel of order 8.
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
0 if abs(x) > 1 and ( 315/4096 ) * ( 1 - x^2 ) * ( 35 - 385 * x^2 + 1001 * x^4 - 715 * x^6 ) otherwise.
}
\details{
A kernel L is a kernel of order 8 if it integrates 1, 
the integrals of u^j L(u) is 0 for 1 <= j < 8 (j integer) and
the integral of u^8 L(u) is different from 0.
}
\author{Alejandra Martinez, Matias Salibian-Barrera.
}