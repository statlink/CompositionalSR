\name{Compositional regression with compositional predictors using the alpha-transformation}
\alias{alfa.pcreg}
\title{
Compositional regression with compositional predictors using the \eqn{\alpha}-transformation
}

\description{
Compositional regression with compositional predictors using the \eqn{\alpha}-transformation.
}

\usage{
alfa.pcreg(y, x, a, k, xnew = NULL, yb = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional responses. Zero values are allowed.
}
\item{x}{
A matrix with the compositional predictors. Zero values are allowed.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it
has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the
solution exists in a closed form, since it the classical mutivariate regression.
}
\item{k}{
How many principal components to compute?
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{yb}{
If you have already transformed the data using the \eqn{\alpha}-transformation with the same
\eqn{\alpha} as given in the argument "a", put it here. Othewrise leave it NULL.
}
}

\details{
The \eqn{\alpha}-transformation is applied to both the compositional responses and predictors.
Then, principal component analysis is performed in the \eqn{\alpha}-transformed predictors and
the projected scores are used a predictors. The same value of \eqn{\alpha} is used for both the
responses and the predictors.
}

\value{
A list including:
\item{runtime}{
The time required by the regression.
}
\item{be}{
The beta coefficients.
}
\item{dev}{
The sum of the squared residuals, as produced by the function \bold{minpack.lm::nls.lm()}.
}
\item{est}{
The fitted values for xnew if xnew is not NULL.
}
}

\references{
Tsagris M. and Pantazis Y. (2026). The \eqn{\alpha}--regression for compositional data: a unified framework for standard, temporal and spatial regression models including compositional predictors.
https://arxiv.org/pdf/2510.12663

Tsagris M. (2015). Regression analysis with compositional data containing zero values.
Chilean Journal of Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for
compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Mardia K.V., Kent J.T., and Bibby J.M. (1979). Multivariate analysis. Academic press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{cv.alfapcreg}, \link{areg}
}
}

\examples{
data(fadn)
y <- fadn[, 3:7]
x <- fadn[, 8:11]
x <- x / rowSums(x)
mod <- alfa.pcreg(y, x, k = 3, 0.2)
}

