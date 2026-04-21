\name{Hypothesis test for the alpha-regression}
\alias{alfa.betest}
\title{
Hypothesis test for the \eqn{\alpha}-regression
}

\description{
Hypothesis test for the \eqn{\alpha}-regression.
}

\usage{
alfa.betest(y, x, a, index = "all", R = 499)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor
variables.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it
has to be greater than 0.
}
\item{index}{
Whose variable the coefficients to test if they are zero? Either "all", to perform a global (omnibus)
test or a number (or a vector) indicating the variable(s).
}
\item{R}{
The number of permutations to perform.
}
}

\details{
The \eqn{\alpha}-transformation is applied to the compositional data first and then multivariate
regression is applied. Then a permutation-based is applied where the test statistic is the SSE if
an omnibus test is performed or the difference in the SSE between the null and the full model in
case of a partial test.
}

\value{
The permutation based p-value.
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
\code{\link{cv.alfareg}, \link{alfareg.nr}, \link{alfa.slx}
}
}

\examples{
data(fadn)
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- areg(y, x, 1)
}

