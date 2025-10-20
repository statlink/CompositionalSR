\name{Prediction with the GWalphaR model}
\alias{gwar.pred}
\title{
Prediction with the GW\eqn{\alpha}R model
}

\description{
Prediction with GW\eqn{\alpha}R model.
}

\usage{
gwar.pred(y, x, a, coords, h, xnew, coordsnew)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
A vector with values for the power transformation, it has to be between -1 and 1.
}
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{h}{
A vector with bandwith values.
}
\item{xnew}{
The new data.
}
\item{coordsnew}{
A matrix with the coordinates of the new locations. The first column is the latitude and the second is the longitude.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the GW\eqn{\alpha}R model is applied and predictions are given for each observation.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{est}{
A list with the fitted values, for each combination of \eqn{\alpha} and h.
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
coords <- fadn[-c(1:10), 1:2]
y <- fadn[-c(1:10), 3:7]
x <- fadn[-c(1:10), 8]
xnew <- fadn[1:10, 8]
coordsnew <- fadn[1:10, 1:2]
mod <- gwar.pred(y, x, a = c(0.25, 0.5, 1), coords,
h = c(0.002, 0.006), xnew = xnew, coordsnew = coordsnew)
}

