\name{gofHybrid}
\alias{gofHybrid}

\title{Hybrid gof test
}
\description{
\code{\link{gofHybrid}} combines all tests in this package to perform the hybrid test presented in Zhang et al. (2015). It also allows for the inclusion of customized tests. The test gives the possibility to combine several single tests which is helpful since in different test scenarios different tests are most powerful.
}
\usage{
gofHybrid(copula, x, tests = c("gofPIOSRn", "gofKernel"), customTests, 
          param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, 
          margins = "ranks", M = 1000, MJ = 100, dispstr = "ex", 
          m = 1, delta.J = 0.5, nodes.Integration = 12, 
          seed.active = NULL, processes = 1)
}
\arguments{
\item{copula}{
The copula to test for. Possible are "normal", "t", "clayton", "gumbel" and "frank".
}
  \item{x}{
A matrix containing the data. 
}
  \item{tests}{
The tests to be used. Possible are \code{"gofPIOSRn"}, \code{"gofPIOSTn"}, \code{"gofKernel"}, \code{"gofRosenblattSnB"}, \code{"gofRosenblattSnC"}, \code{"gofRosenblattChisq"}, \code{"gofRosenblattGamma"}, \code{"gofSn"}, \code{"gofKendallCvM"}, \code{"gofKendallKS"}, \code{"gofWhite"}.
}
\item{customTests}{
A character vector which indicates the customized test to use, if any. The test has to be loaded into the workspace. Currently the function containing the test has to have 2 arguments, the first one for the dataset and the second one for the copula to test for. The arguments have to be named "x" and "copula" respectively.
}
\item{param}{
The copulae parameters to use for each test, if it shall not be estimated.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated.
}
\item{df}{
The degrees of freedom, if not meant to be estimated. Only necessary if tested for \code{"t"}-copula. For the \code{"gofPIOSTn"} test the entry is limited to 60 degrees of freedom for computational reasons. 
}
\item{df.est}{
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated. For the \code{"gofPIOSTn"} test the estimate is limited to 60 degrees of freedom for computational reasons.
}
\item{margins}{
Specifies which estimation method for the margins shall be used. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively the following distributions can be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}). Input can be either one method, e.g. \code{"ranks"}, which will be used for estimation of all data sequences. Also an individual method for each margin can be specified, e.g. \code{c("ranks", "norm", "t")} for 3 data sequences. If one does not want to estimate the margins, set it to \code{NULL}.
}
\item{M}{
Number of bootstrapping samples in the single tests.
}
\item{MJ}{
Size of bootstrapping sample. Only necessary if the test \code{gofKernel} is part of \code{tests}.
}
\item{dispstr}{
A character string specifying the type of the symmetric positive definite matrix characterizing the elliptical copula. Implemented structures are "ex" for exchangeable and "un" for unstructured, see package \code{copula}.
}
\item{m}{
Length of blocks. Only necessary if the test \code{gofPIOSTn} is part of \code{tests}.
}
\item{delta.J}{
Scaling parameter for the matrix of smoothing parameters. Only necessary if the test \code{gofKernel} is part of \code{tests}.
}
\item{nodes.Integration}{
Number of knots of the bivariate Gauss-Legendre quadrature. Only necessary if the test \code{gofKernel} is part of \code{tests}.
}
\item{seed.active}{
Has to be either an integer or a vector of M+1 integers. If an integer, then the seeds for the bootstrapping procedure will be simulated. If M+1 seeds are provided, then these seeds are used in the bootstrapping procedure. Defaults to \code{NULL}, then \code{R} generates the seeds from the computer runtime. Controlling the seeds is useful for reproducibility of a simulation study to compare the power of the tests or for reproducibility of an empirical study.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors. Please see the details.}
}
\details{
In most of scenarios for goodness-of-fit tests, including the one for copula models (e.g. Genest et al. (2009)) there exists no single dominant optimal test. Zhang et al. (2015) proposed a hybrid test which performed in their simulation study more desirably compared to the applied single tests.

The p-value is a combination of the single tests in the following way:
\deqn{p_n^{hybrid} = \min(q \cdot \min{(p_n^{(1)}, \dots, p_n^{(q)})}, 1)}{pn^(hybrid) = min(q x min(pn^(1), ..., pn^(q)), 1)}
where \eqn{q}{q} is the number of tests and \eqn{p_n^{(i)}}{pn^(i)} the p-value of the test \eqn{i}. It is ensured that the hybrid test is consistent as long as at least one of the tests is consistent.

The computation of the individual p-values is performed as described in the details of this tests. Note that the derivation differs.

For small values of \code{M}, initializing the parallelisation via \code{processes} does not make sense. The registration of the parallel processes increases the computation time. Please consider to enable parallelisation just for high values of \code{M}.
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
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}
}
\examples{
data(IndexReturns2D)

gofHybrid("normal", IndexReturns2D, 
tests = c("gofRosenblattSnB", "gofRosenblattSnC"), M = 10)
}
