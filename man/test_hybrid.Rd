\name{gofHybrid}
\alias{gofHybrid}

\title{Hybrid gof test
}
\description{
\code{\link{gofHybrid}} combines all tests in this package to perform the hybrid test presented in Zhang et al. (2015). The test gives the possibility to combine several single tests which is helpful since in different test scenarios are different tests the most powerful.
}
\usage{
gofHybrid(copula, x, testset = c("gofPIOSRn", "gofKernel"), margins = "ranks", 
          M = 1000, execute.times.comp = T, param = 0.5, param.est = T, df = 4, 
          df.est = T, m = 1, MJ = 100, delta.J = 0.5, nodes.Integration = 12, 
          m_b = 0.5, zeta.m = 0, b_Rn = 0.05)
}
\arguments{
\item{copula}{
The copula to test for. Possible are "gaussian", "t", "clayton", "gumbel" and "frank".
}
  \item{x}{
A matrix containing the residuals of the data. Take into account that the dimension of the matrix must fit the required dimensions from the single tests.
}
  \item{testset}{
The tests to be used. Possible are \code{"gofPIOSRn"}, \code{"gofPIOSTn"}, \code{"gofKernel"}, \code{"gofRosenblattSnB"}, \code{"gofRosenblattSnC"}, \code{"gofADChisq"}, \code{"gofADGamma"}, \code{"gofSn"}, \code{"gofKendallCvM"}, \code{"gofKendallKS"}, \code{"gofWhite"}, \code{"gofRn"}.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{M}{
Number of bootstrapping samples in the single tests.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
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
\item{MJ}{
Size of bootstrapping sample. Only necessary if the test \code{gofKernel} is part of \code{testset}.
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
In most of scenarios for goodness-of-fit tests, including the one for copula models (e.g. Genest et al. (2009)) there exists no single dominate optimal test. Zhang et al. (2015) proposed a hybrid test which performed in their simulation study more desirably compared to the applied single tests.

The p-value is a combination of the single tests in the following way:
\deqn{p_n^{hybrid} = q \cdot \min{(p_n^{(1)}, \dots, p_n^{(q)})}}{pn^(hybrid) = q x min(pn^(1), ..., pn^(q))}
where \eqn{q}{q} is the number of tests and \eqn{p_n^{(i)}}{pn^(i)} the p-value of the test \eqn{i}. It is ensured that the hybrid test is consistent as long as at least one of the tests is consistent.

The computation of the individual p-values is performed as described in the details of this tests. Note that the derivation differs.
}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-values and test statistics of the hybrid and the individual tests}
}
\references{
Zhang, S., Okhrin, O., Zhou, Q., and Song, P.. Goodness-of-fit Test For Specification of Semiparametric Copula Dependence Models. \emph{under revision in Journal of Econometrics from 15.01.2014} \url{http://sfb649.wiwi.hu-berlin.de/papers/pdf/SFB649DP2013-041.pdf} \cr \cr
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}
}
\examples{
data = cbind(rnorm(100), rnorm(100))

gofHybrid("gaussian", data, testset = c("gofRosenblattSnB", "gofRosenblattSnC"), M = 10)
}
