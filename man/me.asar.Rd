\name{Marginal effects for the alpha-SAR model}
\alias{me.aslx}
\title{
Marginal effects for the \eqn{\alpha}-SAR model
}

\description{
Marginal effects for the \eqn{\alpha}-SAR model.
}

\usage{
me.asar(be, rho, mu, x, coords, k, cov_theta = NULL) 
}

\arguments{
\item{be}{
A matrix with the beta coefficients of the \eqn{\alpha}-SAR model.
}
\item{rho}{
The spatial auto-regressive coefficient \eqn{\rho} of the \eqn{\alpha}-SAR model.
}
\item{mu}{
The fitted values of the \eqn{\alpha}-SAR model.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame. Categorical predictor variables are not suited here.
}
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{k}{
The number of nearest neighbours to consider for the contiguity matrix.
}
\item{cov_theta}{
The covariance matrix of the beta and gamma regression coefficients. If you pass this argument, then the standard error of the average marginal effects will be returned.
}
}

\details{
The marginal effects of the \eqn{\alpha}-SAR model are computed.
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
\item{se.amedir}{
An array with the standard errors of the average direct marginal effects of each component for each predictor variable. This is returned if you supply the covariance matrix cov_theta.
}
\item{se.ameindir}{
An array with the standard errors of the average indirect marginal effects of each component for each predictor variable. This is returned if you supply the covariance matrix cov_theta.
}
\item{se.ametotal}{
An array with the standard errors of the average total marginal effects of each component for each predictor variable. This is returned if you supply the covariance matrix cov_theta.
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
\code{\link{me.ar}, \link{me.aslx}, \link{me.gwar}
}
}

\examples{
data(fadn)
coords <- fadn[, 1:2]
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfa.sar(y, x, a = 0.75, coords, k = 6)
me <- me.asar(mod$be, mod$rho, mod$est, x, coords, k = 6) 
}

