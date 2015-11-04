\name{gof}
\alias{gof}

\title{Combining function for tests
}
\description{
\code{\link{gof}} computes for a given dataset and based on the choices of the user either all tests for a given amount of copulae, performs for a given testset every test with all available copulae or computes for given copulae and tests all possible combinations.
}
\usage{
gof(x, priority = "tests", copula = NULL, tests = NULL, M = 50, MJ = 50, 
    param = 0.5, param.est = T, df = 4, df.est = T, m = 1, delta.J = 0.5, 
    nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05)
}
\arguments{
  \item{x}{
A matrix containing the residuals of the data.
}
  \item{priority}{
A character string which is either \code{"tests"} or \code{"copula"}. \code{"tests"} indicates that all implemented tests are performed for all copulae which the tests share. This are e.g. \code{"gaussian"} and \code{"clayton"}. If \code{"copula"} is chosen, the tests which are able to test for \code{"gaussian"}, \code{"t"}, \code{"frank"}, \code{"gumbel"} and \code{"clayton"} are performed. If one of the arguments \code{tests} or \code{copula} is not \code{NULL}, then \code{priority} doesn't affect the choice of the copulae and tests.
}
\item{copula}{
A character vector which indicates the copula to test for.
}
\item{tests}{
A character vector which indicates the tests to use.
}
\item{M}{
The amount of bootstrap rounds to be performed by each test. Default is 1000.
}
\item{MJ}{
Just for the test gofKernel. Size of bootstrapping sample.
}
\item{param}{
The copulae parameters to use for each test, if it shall not be estimated.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated.
}
\item{df}{
The degrees of freedom, if not meant to be estimated. Only necessary if tested for \code{"t"}-copula. 
}
\item{df.est}{
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated.
}
\item{m}{
Length of blocks. Only necessary if the test \code{gofPIOSTn} is part of \code{testset}.
}
\item{delta.J}{
Scaling parameter for the matrix of smoothing parameters. Only necessary if the test \code{gofKernel} is part of \code{testset}.
}
\item{nodes.Integration}{
Number of knots of the bivariate Gauss-Legendre quadrature. Only necessary if the test \code{gofKernel} is part of \code{testset}.
}
\item{m_b}{
The power of the statistic. Only necessary if the test \code{gofRn} is part of \code{testset}.
}
\item{zeta.m}{
The adjustment parameter. Only necessary if the test \code{gofRn} is part of \code{testset}.
}
\item{b_Rn}{
The bandwidth for the estimation of the first-order partial derivatives based on the empirical copula. Only necessary if the test \code{gofRn} is part of \code{testset}.
}
}
\details{
If a character vector is given for the argument \code{copula} and nothing for \code{tests}, then all tests are performed for which the given copulae are implemented. If \code{tests} contains a character vector of tests and \code{copula = NULL}, then will be this tests performed for all implemented copulae. If character vectors are given for \code{copula} and \code{tests}, then the tests are performed with the given copulae. If \code{tests = NULL} and \code{copula = NULL}, then the argument \code{priority} catches in and defines the procedure.

The time to compute the entire procedure is always estimated in case that \code{M} or \code{MJ} are 100 or higher.
} 
\value{
A list containing several objects of \code{class} gofCOP with the following components for each copulae
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-values and test statistics of the individual tests including the hybrid test}
}
\examples{
data = cbind(rnorm(100), rnorm(100), rnorm(100))

gof(data, priority = "tests", copula = "gaussian", 
tests = c("gofRosenblattSnB", "gofRosenblattSnC"), M = 20)
}
