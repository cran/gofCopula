\name{gofCopula4Test}
\alias{gofCopula4Test}

\title{Implemented copula for a certain test
}
\description{
\code{\link{gofCopula4Test}} returns for a given test the applicable implemented copula.
}
\usage{
gofCopula4Test(test)
}
\arguments{
\item{test}{
The test to search for copula.
}
}
\details{
In case that the decision for a certain gof test was already done, it is interesting to know which copula can be used with this test. }
\value{
A character vector which consists of the names of the copula.
}
\examples{
gofCopula4Test("gofRosenblattSnB")

gofCopula4Test("gofPIOSTn")
}
