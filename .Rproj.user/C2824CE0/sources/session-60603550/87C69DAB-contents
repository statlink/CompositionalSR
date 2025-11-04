\name{The alpha-regression using Newton-Raphson}
\alias{alfareg.nr}
\title{
The \eqn{alpha}-regression using Newton-Raphson
}

\description{
The \eqn{alpha}-regression using Newton-Raphson.
}

\usage{
alfareg.nr(y, x, alpha = 1, beta_init = NULL, max_iter = 100,
tol = 1e-6, line_search = TRUE)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{alpha}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0.
}
\item{beta_init}{
A vector of initial parameters (optional). This is then transformed into a matrix.
}
\item{max_iter}{
The maximum number of iterations for the Newton-Raphson algorithm.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{line_search}{
Do you want to perform line search? The default value is TRUE.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then multivariate regression is applied. This involves numerical optimisation.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{iters}{
The iterations of the Newton-Raphson algorithm
}
\item{be}{
The beta coefficients.
}
\item{objective}{
The sum of the squared residuals.
}
\item{est}{
The fitted values.
}
\item{covb}{
The covariance matrix of the beta coefficients, or NULL if it is singular.
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

Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.

Aitchison J. (1986). The statistical analysis of compositional data. Chapman & Hall.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.reg}, \link{cv.alfareg}, \link{alfa.slx}
}
}

\examples{
data(fadn)
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfareg.nr(y, x, a = 0.2)
mod2 <- alfa.reg(y, x, 0.2)
}

