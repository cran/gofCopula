\name{gofTest4Copula}
\alias{gofTest4Copula}

\title{Applicable gof tests for testing problem
}
\description{
\code{\link{gofTest4Copula}} returns for a given copula and a dimension the applicable implemented tests.
}
\usage{
gofTest4Copula(copula = NULL, d = 2)
}
\arguments{
\item{copula}{
The copula to test for. Possible are the copulae \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}. For the default \code{NULL}, all tests in the package are returned. 
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
gofTest4Copula("clayton", d = 2)

gofTest4Copula("gumbel", d = 5)
}
