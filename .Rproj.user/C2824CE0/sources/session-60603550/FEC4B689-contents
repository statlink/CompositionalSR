\name{Moran's I type test for compositional data using the distance correlation}
\alias{moran.adcor}
\title{
Moran's I type test for compositional data using the distance correlation
}

\description{
Moran's I type test for compositional data using the distance correlation.
}

\usage{
moran.adcor(y, a, W, R = 999)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied. If it is a vector, the test is performed for each value of \eqn{\alpha}.
}
\item{W}{
The contiguity matrix.
}
\item{R}{
The number of permutations to perform.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the distance covariance test of independence is performed.
}

\value{
Depending on the length of the argument "a", a single p-value or a vector of p-values, one for each value of \eqn{\alpha}.
}

\references{
Tsagris M. (2025). The \eqn{\alpha}--regression for compositional data: a unified framework for standard, spatially-lagged, spatial autoregressive and geographically-weighted regression models.
https://arxiv.org/pdf/2510.12663

Tsagris M. (2015). Regression analysis with compositional data containing zero values.
Chilean Journal of Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for
compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{me.aslx}, \link{me.gwar}, \link{alfa.reg}
}
}

\examples{
data(fadn)
y <- fadn[, 3:7]
coords <- fadn[, 1:2]
W <- CompositionalSR::contiguity(coords, k = 10)
moran.adcor(y, 1, W, R = 299)
}

