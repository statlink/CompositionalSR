\name{fadn}
\alias{fadn}
\docType{data}
\title{
FADN dataset
}
\description{
A matrix with 11 columns. The first two are the locations (latitude and longitude), the next five
contain the compositional data (percentages of cultivated area of five crops), Y1.1: cereals, Y2.1: cotton,
Y3.1: tree crops, Y4.1: other annual crops and pasture and Y5.1: grapes and wine. The next four columns contain
the covariates, G1:  Human Influence Index, G2: soil pH, G3: topsoil organic carbon content and G7: erosion.
}
\usage{fadn}
\format{
A matrix with 168 rows and 11 columns.
}
\source{
Clark and Dixon (2021), available at \url{https://github.com/nick3703/Chicago-Data}.
}
\references{
Clark, N. J. and P. M. Dixon (2021). A class of spatially correlated
self-exciting statistical models. \emph{Spatial Statistics}, 43, 1--18.
}
\seealso{
\code{\link{alfa.slx}, \link{gwar}, \link{alfa.reg}
}
}
\examples{
data(fadn)
y <- fadn[, 3:7]
x <- fadn[, 8:11]
mod <- alfa.reg(y, x, a = 0.1)
}
