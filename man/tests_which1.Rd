\name{gofWhich}
\alias{gofWhich}

\title{Applicable gof tests for testing problem
}
\description{
\code{\link{gofWhich}} returns for a given copula and a dimension the applicable implemented tests.
}
\usage{
gofWhich(copula, d)
}
\arguments{
\item{copula}{
The copula to test for. Possible are the copulae \code{"normal"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
}
  \item{d}{
The dimension to search for.
}
}
\details{
Before performing a gof test on a dataset, it pays out to have a close look at the Scatterplot to receive an idea about the possible type of copula. Afterwards follows the decision about the test. The tests in this package can be used for different types of copulae functions and dimensions. This function is dedicated to help finding the applicable gof tests for the dataset.
}
\value{
A character vector which consists of the names of the tests.
}
\examples{
gofWhich("clayton", d = 2)

gofWhich("gumbel", d = 5)
}
