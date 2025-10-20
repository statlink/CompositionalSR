\name{The GWalphaR model}
\alias{gwar}
\title{
The GW\eqn{\alpha}R model
}

\description{
The GW\eqn{\alpha}R model.
}

\usage{
gwar(y, x, a, coords, h, yb = NULL, nc = 1)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1.
}
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{h}{
The bandwith value.
}
\item{yb}{
If you have already transformed the data using the \eqn{\alpha}-transformation with the same \eqn{\alpha} as given in the argument "a", put it here. Othewrise leave it NULL.
}
\item{nc}{
The number of cores to use. IF you have a multicore computer it is advisable to use more than 1. It makes the procedure faster. It is advisable to use it if you have many observations and or many variables, otherwise it will slow down th process.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the GW\eqn{\alpha}R model is applied.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{be}{
The beta coefficients.
}
\item{est}{
The fitted values.
}
}

\references{
Tsagris M. (2025). The \eqn{\alpha}--regression for compositional data: a unified framework for standard, spatially-lagged, and geographically-weighted regression models.
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
\code{\link{cv.gwar}, \link{me.gwar}, \link{alfa.slx}, \link{alfa.reg}
}
}

\examples{
data(fadn)
coords <- fadn[, 1:2]
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- gwar(y, x, a = 1, coords, h = 0.001)
}

