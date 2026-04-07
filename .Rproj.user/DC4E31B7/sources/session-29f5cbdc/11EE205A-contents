\name{The alpha-ESF model}
\alias{alfa.esf}
\title{
The \eqn{\alpha}-ESF model
}

\description{
The \eqn{\alpha}-ESF model.
}

\usage{
alfa.esf(y, x, a, coords, xnew = NULL, coordsnew, yb = NULL)
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
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{coordsnew}{
A matrix with the coordinates of the new locations. The first column is the latitude and the second is the longitude. If you do not have new data to make predictions leave this NULL.
}
\item{yb}{
If you have already transformed the data using the \eqn{\alpha}-transformation with the same \eqn{\alpha} as given in the argument "a", put it here. Othewrise leave it NULL.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first. Then the eigenvectors of the
kernelized distance matrix are computed and the appropriate number is selected to be included as predictors.
the selection takes place using the \eqn{\gamma}-OMP algorithm (Tsagris et al., 2022).
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
The gamma coefficients of the eigenvectors.
}
\item{ESF}{
A vector with the indices of the eigenvectors used.
}
\item{X.esf}{
A matrix with the values of the eigenvectors used.
}
\item{dev}{
The sum of the squared residuals, as produced by the function \bold{minpack.lm::nls.lm()}.
}
\item{est}{
The fitted values for xnew if xnew and coordsnew are not NULL.
}
}

\references{
Tsagris M. and Pantazis Y. (2026). The \eqn{\alpha}--regression for compositional data: a unified framework for standard, spatially-lagged, spatial autoregressive and geographically-weighted regression models.
https://arxiv.org/pdf/2510.12663

Tsagris M. (2015). Regression analysis with compositional data containing zero values.
Chilean Journal of Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for
compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Tsagris M., Papadovasilakis Z., Lakiotaki K. and Tsamardinos I. (2022).
The \eqn{\gamma}-OMP Algorithm for Feature Selection With Application to Gene Expression Data.
IEEE/ACM Transactions on Computational Biology and Bioinformatics, 19(2), 1214--1224.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cv.alfaesf}, \link{alfa.sar}, \link{alfa.slx}, \link{alfa.reg}
}
}

\examples{
data(fadn)
coords <- fadn[1:50, 1:2]
y <- fadn[1:50, 3:7]
x <- fadn[1:50, 8]
mod <- alfa.esf(y, x, a = 1, coords = coords)
}

