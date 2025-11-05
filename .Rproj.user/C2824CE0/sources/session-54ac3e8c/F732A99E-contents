\name{The alpha-MAD-SLX model}
\alias{alfamad.slx}
\title{
The \eqn{\alpha}-MAD-SLX model
}

\description{
The \eqn{\alpha}-MAD-SLX model.
}

\usage{
alfamad.slx(y, x, a, start = NULL, coords, k = 10, xnew = NULL, coordsnew, yb = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the solution exists in a closed form, since it the classical mutivariate regression.
}
\item{start}{
If you want to use initial (or starting) values, pass them here as a vector. If NULL, the coefficients from a multivariate linear model are used.
}
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{k}{
The number of nearest neighbours to consider for the contiguity matrix. For the alfa.slx2() this should be a vector. 
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{coordsnew}{
A matrix with the coordinates of the new locations. The first column is the latitude and the second is the longitude. If you do not have new data to make predictions leave this NULL.
}
\item{yb}{
If you have already transformed the data using the \eqn{\alpha}-transformation with the same \eqn{\alpha} as given in the argument "a", put it here. Othewrise leave it NULL.

This is intended to be used in the function \code{\link{cv.alfareg}} in order to speed up the process. The time difference in that function is small for small samples. But, if you have a few thousands and or a few more components, there will be bigger differences.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the mean absolute deviations 
(quantile regression) spatially lagged X (SLX) model is applied.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{be}{
The beta coefficients.
}
\item{gama}{
The gamma coefficients.
}
\item{est}{
The fitted values for xnew if xnew and coordsnew are not NULL.
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
\code{\link{alfa.mad}, \link{alfa.reg}, \link{gwar}
}
}

\examples{
data(fadn)
coords <- fadn[, 1:2]
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfamad.slx(y, x, a = 0.5, coords, k = 10)
}

