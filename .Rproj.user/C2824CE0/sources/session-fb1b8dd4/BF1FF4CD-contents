\name{The alpha-SAR model}
\alias{alfa.sar}
\alias{alfa.sar2}
\title{
The \eqn{\alpha}-SAR model
}

\description{
The \eqn{\alpha}-SAR model.
}

\usage{
alfa.sar(y, x, a, coords, k = 10, xnew = NULL, coordsnew, yb = NULL)
alfa.sar2(y, x, a, coords, k = 10, xnew = NULL, coordsnew, yb = NULL)
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
\item{k}{
The number of nearest neighbours to consider for the contiguity matrix.
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
The \eqn{\alpha}-transformation is applied to the compositional data first and the spatial autocorrelation (SAR) model is applied. The difference between the two is with regards to the estimation. The first function iterates between the golden ratio for the optimal \eqn{\rho}, given some initial values for the \eqn{\beta_s}, and the Levenberg-Marquardt algorithm to estimate the \eqn{\beta_s}, given a value of \eqn{\rho}. The second function performs a grid search searching for the range of good values of \eqn{\rho} and then uses that as starting value. The first option is faster, but this way, the covariance matrix cannot be estimated correctly. I suggest to the firstoption, you can test it obviously, because it is faster and I think it gives better results.
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
\item{dev}{
The sum of the squared residuals, as produced by the function \bold{minpack.lm::nls.lm()}.
}
\item{est}{
The fitted values for xnew if xnew and coordsnew are not NULL.
}
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
\code{\link{cv.alfaslx}, \link{me.aslx}, \link{gwar}, \link{alfa.reg}
}
}

\examples{
data(fadn)
coords <- fadn[, 1:2]
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfa.sar(y, x, a = 0.5, coords, k = 10)
}

