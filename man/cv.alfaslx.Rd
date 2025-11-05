\name{Spatial K-fold cross-validation for the alpha-SLX model}
\alias{cv.alfaslx}
\title{
Spatial K-fold cross-validation for the \eqn{\alpha}-SLX model
}

\description{
Spatial K-fold cross-validation for the \eqn{\alpha}-SLX model
}

\usage{
cv.alfaslx(y, x, a = seq(0.1, 1, by = 0.1), coords, k = 2:15,
nfolds = 10, size = 1000, folds = NULL)
}

\arguments{
\item{y}{
A matrix with compositional data. zero values are allowed.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied.
}
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{k}{
A vector with the nearest neighbours to consider for the contiguity matrix.
}
\item{nfolds}{
The number of folds to split the data.
}
\item{size}{
A numeric value of the specified range by which blocks are created and training/testing data are separated. This distance should be in metres. If you have big regions you should consider increasing this number. For more information see the package \emph{blockCV}.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
}
\details{
The 10-fold spatial cross-validation protocol is applied to choose the optimal values of \eqn{\alpha} and k.
}

\value{
A list including:
\item{runtime}{
The runtime required by the cross-validation.
}
\item{perf}{
A vector with the average Kullback-Leibler divergence, for every value of \eqn{\alpha}.
}
\item{opt}{
A vector with the minimum Kullback-Leibler divergence, the optimal value of \eqn{\alpha} and k.
}
}

\references{
Tsagris M. (2025). The \eqn{\alpha}--regression for compositional data: a unified framework for standard, spatially-lagged, spatial autoregressive and geographically-weighted regression models.
https://arxiv.org/pdf/2510.12663

Tsagris M. (2015). Regression analysis with compositional data containing zero values.
Chilean Journal of Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data.
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
\code{\link{alfa.slx}, \link{cv.gwar} \link{cv.alfareg}
}
}

\examples{
data(fadn)
coords <- fadn[1:60, 1:2]
y <- fadn[1:60, 3:7]
x <- fadn[1:60, 8]
mod <- cv.alfaslx(y, x, a = 0.5, coords, k = 2)
}

