\name{gofCheckTime}
\alias{gofCheckTime}

\title{Combining function for tests
}
\description{
The computation of a gof test can take very long, especially when the number of bootstrap rounds is high. The function \code{\link{gofCheckTime}} computes the time which the estimation most likely take.
}
\usage{
gofCheckTime(copula, x, test, M = 1000, MJ = 10000, print.res = T, margins = "ranks", 
              dispstr = "ex", param = 0.5, param.est = T, df = 4, df.est = T, 
              m = 1, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, 
              zeta.m = 0, b_Rn = 0.05, processes = 1)
}
\arguments{
\item{copula}{
A character vector which indicates the copula to test for.
}
  \item{x}{
A matrix containing the data.
}
\item{test}{
A character vector which indicates the test to use.
}
\item{M}{
The number of bootstrapping rounds which shall be later taken in the estimation.
}
\item{MJ}{
Just for the test gofKernel. Size of bootstrapping sample.
}
\item{print.res}{
Logical which defines if the resulting time shall be printed or given as value. Default is TRUE.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{dispstr}{
A character string specifying the type of the symmetric positive definite matrix characterizing the elliptical copula. Implemented structures are "ex" for exchangeable, "ar1" for AR(1), "toep" for Toeplitz, and "un" for unstructured, see package \code{copula}.
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
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors.}
}
\details{
The function estimates the time which the entire gof test will take.
} 
\examples{
data(IndexReturns)

gofCheckTime("normal", IndexReturns[c(1:50),c(1:2)], "gofRosenblattSnC", M = 10000)
}
