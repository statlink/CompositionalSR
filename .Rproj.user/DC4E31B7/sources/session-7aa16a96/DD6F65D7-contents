\name{ICE plot for the alpha-regression}
\alias{ice.reg}
\title{
ICE plot for the \eqn{\alpha}-regression
}

\description{
ICE plot for the \eqn{\alpha}-regression.
}

\usage{
ice.areg(be, x, ind = 1, frac = 0.25)
}

\arguments{
\item{be}{
A numerical matrix with the estimated \eqn{\alpha}-regression coefficients.
}
\item{x}{
A numerical matrix with the predictor variables.
}
\item{ind}{
Which variable to select?.
}
\item{frac}{
Fraction of observations to use. The default value is 0.25.
}
}

\details{
This function implements the Individual Conditional Expecation plots of Goldstein et al. (2015).
See the references for more details.
}

\value{
A graph with several curves, one for each component. The horizontal axis contains the selected variable,
whereas the vertical axis contains the locally smoothed predicted compositional lines.
}

\references{
https://christophm.github.io/interpretable-ml-book/ice.html

Goldstein, A., Kapelner, A., Bleich, J. and Pitkin, E. (2015).
Peeking inside the black box: Visualizing statistical
learning with plots of individual conditional expectation.
Journal of Computational and Graphical Statistics 24(1): 44-65.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{areg}, \link{cv.alfareg}
}
}

\examples{
data(fadn)
y <- fadn[, 3:7]
x <- fadn[, 8]
mod <- alfa.reg(y, x, 0.2)
ice <- ice.areg(y, x, a = 0.2, ind = 1)
}


