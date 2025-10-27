\name{Marginal effects for the alpha-SAR model}
\alias{me.asar}
\title{
Marginal effects for the \eqn{\alpha}-SAR model
}

\description{
Marginal effects for the \eqn{\alpha}-SAR model.
}

\usage{
me.asar(be, mu, x, rho, coords, k = 10)
}

\arguments{
\item{be}{
A matrix with the beta coefficients of the \eqn{\alpha}-SAR model.
}
\item{mu}{
The fitted values of the \eqn{\alpha}-SAR model.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame. Categorical predictor variables are not suited here.
}
\item{rho}{
The spatial autoregressive coefficient.
}
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{k}{
The number of nearest neighbours to consider for the contiguity matrix.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then the \eqn{\alpha}-SLX model is applied.
}

\value{
A list including:
\item{me.dir}{
An array with the direct marginal effects of each component for each predictor variable.
}
\item{me.indir}{
An array with the indirect marginal effects of each component for each predictor variable.
}
\item{me.total}{
An array with the total marginal effects of each component for each predictor variable.
}
\item{ame.dir}{
An array with the average direct marginal effects of each component for each predictor variable.
}
\item{ame.indir}{
An array with the average indirect marginal effects of each component for each predictor variable.
}
\item{ame.total}{
An array with the aerage total marginal effects of each component for each predictor variable.
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
\code{\link{me.aslx}, \link{me.ar}, \link{me.gwar}
}
}

\examples{
data(fadn)
coords <- fadn[, 1:2]
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfa.sar(y, x, a = 0.5, coords, k = 10, xnew = x, coordsnew = coords)
me <- me.aslx(mod$be, mod$est, x, mod$rho, coords, k = 10)
}

