\name{Computation of the contiguity matrix W}
\alias{contiguity}
\title{
Computation of the contiguity matrix W
}

\description{
Computation of the contiguity matrix W.
}

\usage{
contiguity(coords, k = 10)
}

\arguments{
\item{coords}{
A matrix with the coordinates of the locations. The first column is the latitude and the second is the longitude.
}
\item{k}{
The number of nearest neighbours to consider for the contiguity matrix.
}
}

%\details{
%}

\value{
The contiguity matrix W. A square matrix with row standardised values (the elements of each row sum to 1).
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{alfa.slx}, \link{cv.alfaslx}, \link{me.aslx}
}
}

\examples{
data(fadn)
W <- contiguity(fadn[, 1:2])
}

