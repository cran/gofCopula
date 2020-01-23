\name{gofOutputHybrid}
\alias{gofOutputHybrid}

\title{Output Hybrid gof test
}
\description{
\code{\link{gofOutputHybrid}} outputs the desired Hybrid tests from previous test results from this package for the specified testing size.
}
\usage{
gofOutputHybrid(result, tests = NULL, nsets = NULL)
}
\arguments{
\item{result}{
An object of \code{class} gofCOP.}

\item{tests}{
Individual tests which should be used in the hybrid test. Submit a vector containing the position of the individual tests as they appear in the object submitted, e.g. \code{c(1,4)} for the 1st and 4th tests. If \code{tests} is set NULL (default), all possible testing sizes are returned.}

\item{nsets}{
The desired number of tests to be included in each hybrid test. It should be an integer larger than 1 and smaller or equal than the number of tests given in \code{result}. If \code{nsets} is set NULL (default), all possible testing sizes are calculated.}
}

\details{
For details about hybrid testing please see the detail section in \code{\link{gofHybrid}}.
}

\value{
An object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{copula}{the copula tested for}
\item{margins}{the method used to estimate the margin distribution.}
\item{param.margins}{the parameters of the estimated margin distributions. Only applicable if the margins were not specified as \code{"ranks"} or \code{NULL}.}
\item{theta}{dependence parameters of the copulae}
\item{df}{the degrees of freedem of the copula. Only applicable for t-copula.}
\item{res.tests}{a matrix with the p-values and test statistics of the hybrid and the individual tests}
}
\references{
Zhang, S., Okhrin, O., Zhou, Q., and Song, P.. Goodness-of-fit Test For Specification of Semiparametric Copula Dependence Models. \emph{Journal of Econometrics, 193, 2016, pp. 215-233} \url{https://doi.org/10.1016/j.jeconom.2016.02.017} \cr \cr
}
\examples{
data(IndexReturns2D)

res1 = gof(IndexReturns2D, priority = "tests", copula = "normal", 
           tests = c("gofKendallCvM", "gofRosenblattSnC", "gofKendallKS"), M = 10)
gofOutputHybrid(res1, tests = 1, nsets = 2)
# mind the difference to the regular output
res1
}
