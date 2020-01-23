\name{gofGetHybrid}
\alias{gofGetHybrid}

\title{GetHybrid gof test
}
\description{
\code{\link{gofGetHybrid}} computes based on previous test results from this package and on p-values from your own goodness-of-fit tests the hybrid test p-values for the specified testing size.
}
\usage{
gofGetHybrid(result, p_values = NULL, nsets = NULL)
}
\arguments{
\item{result}{
An object of \code{class} gofCOP.}

\item{p_values}{
A vector containing different p-values from your own goodness-of-fit tests. If the elements are unnamed, the test results in the output will be labeled with Test_A, Test_B etc.}

\item{nsets}{
The desired number of tests to be included in each hybrid test. It should be an integer larger than 1 and smaller or equal than the number of tests given in \code{result} and \code{p_values} together. If \code{nsets} is set NULL (default), all possible testing sizes are calculated.}
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

res_2 = gof(x = IndexReturns2D, copula = "normal", 
            tests = c("gofKernel", "gofKendallCvM", "gofWhite"), M = 10)
gofGetHybrid(result = res_2, p_values = c("MyTest" = 0.3, "AnotherTest" = 0.7), nsets = 3)
}