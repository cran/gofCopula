\name{gofCheckTime}
\alias{gofCheckTime}

\title{Combining function for tests
}
\description{
The computation of a gof test can take very long, especially when the number of bootstrap rounds is high. The function \code{\link{gofCheckTime}} computes the time which the estimation most likely takes.
}
\usage{
gofCheckTime(copula, x, tests, customTests, param = 0.5, param.est = TRUE, df = 4, 
              df.est = TRUE, margins = "ranks", M = 1000, MJ = 100, dispstr = "ex", 
              print.res = TRUE, m = 1, delta.J = 0.5, nodes.Integration = 12, 
              seed.active = NULL, processes = 1)
}
\arguments{
\item{copula}{
A character vector which indicates the copula to test for.
}
  \item{x}{
A matrix containing the data.
}
\item{tests}{
A character vector which indicates the test to use.
}
\item{customTests}{
A character vector which indicates the customized test to use, if any.
}
\item{print.res}{
Logical which defines if the resulting time shall be printed or given as value. Default is TRUE.
}
\item{param}{
The copulae parameters to use for each test, if it shall not be estimated.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated.
}
\item{df}{
The degrees of freedom, if not meant to be estimated. Only necessary if tested for \code{"t"}-copula.  For the \code{"gofPIOSTn"} test the entry is limited to 60 degrees of freedom for computational reasons. 
}
\item{df.est}{
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated. For the \code{"gofPIOSTn"} test the estimate is limited to 60 degrees of freedom for computational reasons.
}
\item{margins}{
Specifies which estimation method for the margins shall be used. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively the following distributions can be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}). Input can be either one method, e.g. \code{"ranks"}, which will be used for estimation of all data sequences. Also an individual method for each margin can be specified, e.g. \code{c("ranks", "norm", "t")} for 3 data sequences. If one does not want to estimate the margins, set it to \code{NULL}.
}
\item{M}{
The number of bootstrapping rounds which shall be later taken in the estimation.
}
\item{MJ}{
Just for the test gofKernel. Size of bootstrapping sample.
}
\item{dispstr}{
A character string specifying the type of the symmetric positive definite matrix characterizing the elliptical copula. Implemented structures are "ex" for exchangeable and "un" for unstructured, see package \code{copula}.
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
\item{seed.active}{
Has to be either an integer or a vector of M+1 integers. If an integer, then the seeds for the bootstrapping procedure will be simulated. If M+1 seeds are provided, then these seeds are used in the bootstrapping procedure. Defaults to \code{NULL}, then \code{R} generates the seeds from the computer runtime. Controlling the seeds is useful for reproducibility of a simulation study to compare the power of the tests or for reproducibility of an empirical study.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors.}
}
\details{
The function estimates the time which the entire gof test will take.
} 
\examples{
data(IndexReturns2D)

gofCheckTime("normal", IndexReturns2D, "gofRosenblattSnC", M = 10000)
}
