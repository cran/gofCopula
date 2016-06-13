\name{gofco}
\alias{gofco}

\title{Interface with copula class
}
\description{
\code{\link{gofco}} is an interface with the \code{copula} package. It reads out the information from a copula class object and hands them over to a specified gof test or set of tests.
}
\usage{
gofco(copulaobject, x, testset = c("gofPIOSRn", "gofKernel"), margins = "ranks", 
      dispstr = "ex", M = 1000, execute.times.comp = T, m = 1, MJ = 100, 
      delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05, 
      processes = 1)
}
\arguments{
\item{copulaobject}{
An object with the copula class from the copula package.
}
  \item{x}{
A matrix containing the residuals of the data. Take into account that the dimension of the matrix must fit the required dimensions from the single tests.
}
  \item{testset}{
The tests to be used. Possible are \code{"gofPIOSRn"}, \code{"gofPIOSTn"}, \code{"gofKernel"}, \code{"gofRosenblattSnB"}, \code{"gofRosenblattSnC"}, \code{"gofRosenblattChisq"}, \code{"gofRosenblattGamma"}, \code{"gofSn"}, \code{"gofKendallCvM"}, \code{"gofKendallKS"}, \code{"gofWhite"}, \code{"gofRn"}.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{dispstr}{
A character string specifying the type of the symmetric positive definite matrix characterizing the elliptical copula. Implemented structures are "ex" for exchangeable, "ar1" for AR(1), "toep" for Toeplitz, and "un" for unstructured, see package \code{copula}.
}
\item{M}{
Number of bootstrapping samples in the single tests.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
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
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors.}
}
\details{
The function reads out the arguments in the copula class object. If the dependence parameter is not specified in the object, it is estimated. In case that the object describes a "t"-copula, then the same holds for the degrees of freedom. The dimension is not extracted from the object. It is obtained from the inserted dataset.

When more than one test shall be performed, the hybrid test is computed too.
}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-values and test statistics of the hybrid and the individual tests}
}
\references{
Yan, Jun. Enjoy the joy of copulas: with a package copula. Journal of Statistical Software 21.4 (2007): 1-21.
}
\examples{
data(IndexReturns)
copObject = normalCopula()

gofco(copObject, x = IndexReturns[c(1:100),c(1:2)], testset = c("gofPIOSRn", "gofKernel"), M = 20)
}
