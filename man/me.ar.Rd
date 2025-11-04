\name{Marginal effects for the alpha-regression model}
\alias{me.ar}
\title{
Marginal effects for the \eqn{\alpha}-regression model
}

\description{
Marginal effects for the \eqn{\alpha}-regression model.
}

\usage{
me.ar(be, mu, x, cov_be = NULL)
}

\arguments{
\item{be}{
A matrix with the beta regression coefficients of the \eqn{\alpha}-regression model.
}
\item{mu}{
The fitted values of the \eqn{\alpha}-regression.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame. Categorical predictor variables are not suited here.
}
\item{cov_be}{
The covariance matrix of the beta regression coefficients. If you pass this argument, then the standard error of the average marginal effects will be returned.
}
}

\details{
The marginal effects of the \eqn{\alpha}-regression model are computed.
}

\value{
A list including:
\item{me}{
An array with the marginal effects of each component for each predictor variable.
}
\item{ame}{
The average marginal effects of each component for each predictor variable.
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
\code{\link{me.aslx}, \link{me.gwar}, \link{alfa.reg}
}
}

\examples{
data(fadn)
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfa.reg(y, x, 0.2, xnew = x)
me <- me.ar(mod$be, mod$est, x)
}

