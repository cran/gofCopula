% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests_gof.R
\name{gof}
\alias{gof}
\title{Combining function for tests}
\usage{
gof(
  x,
  priority = "copula",
  copula = NULL,
  tests = NULL,
  customTests = NULL,
  param = 0.5,
  param.est = TRUE,
  df = 4,
  df.est = TRUE,
  margins = "ranks",
  flip = 0,
  M = 1000,
  MJ = 100,
  dispstr = "ex",
  m = 1,
  delta.J = 0.5,
  nodes.Integration = 12,
  lower = NULL,
  upper = NULL,
  seed.active = NULL,
  processes = 1
)
}
\arguments{
\item{x}{A matrix containing the data with rows being observations and
columns being variables.}

\item{priority}{A character string which is either \code{"tests"} or
\code{"copula"}. \code{"tests"} indicates that all implemented tests are
performed for all copulae which the tests share. These are e.g.
\code{"normal"} and \code{"clayton"}. If \code{"copula"} is chosen, 
all copula are chosen and only the tests are performend which these copula 
share. If one of the arguments \code{tests} or \code{copula} is not 
\code{NULL}, then \code{priority} doesn't affect the choice of the copulae 
and tests.}

\item{copula}{A character vector which indicates the copula to test for.
Possible are \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"}, 
\code{"frank"}, \code{"joe"}, \code{"amh"}, \code{"galambos"}, 
\code{"huslerReiss"}, \code{"tawn"}, \code{"tev"}, \code{"fgm"} and 
\code{"plackett"}.}

\item{tests}{A character vector which indicates the tests to use. Possible
choices are the individual tests implemented in this package.}

\item{customTests}{A character vector which indicates the customized test to
use, if any. The test has to be loaded into the workspace. Currently the
function containing the test has to have 2 arguments, the first one for the
dataset and the second one for the copula to test for. The arguments have to
be named "x" and "copula" respectively.}

\item{param}{The copulae parameters to use for each test, if it shall not be
estimated.}

\item{param.est}{Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE}
means that \code{param} will be estimated.}

\item{df}{The degrees of freedom, if not meant to be estimated. Only
necessary if tested for \code{"t"}-copula. For the \code{"gofPIOSTn"} test
the entry is limited to 60 degrees of freedom for computational reasons.}

\item{df.est}{Indicates if \code{df} shall be estimated. Has to be either
\code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be
estimated. For the \code{"gofPIOSTn"} test the estimate is limited to 60
degrees of freedom for computational reasons.}

\item{margins}{Specifies which estimation method for the margins shall be
used. The default is \code{"ranks"}, which is the standard approach to
convert data in such a case. Alternatively the following distributions can
be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}),
\code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal
(\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
Input can be either one method, e.g. \code{"ranks"}, which will be used for
estimation of all data sequences. Also an individual method for each margin
can be specified, e.g. \code{c("ranks", "norm", "t")} for 3 data sequences.
If one does not want to estimate the margins, set it to \code{NULL}.}

\item{flip}{The vector of control parameters to flip the copula 
by 90, 180, 270 degrees clockwise. Only applicable for bivariate copula. 
Default is 0 and possible inputs are 0, 90, 180, 270 and NULL. One can either
specify one flip degree which will be applied on all copulae or choose
an individual flip for each copula in which case the input has to be a vector.}

\item{M}{The amount of bootstrap rounds to be performed by each test.
Default is 1000.}

\item{MJ}{Just for the test gofKernel. Size of bootstrapping sample.}

\item{dispstr}{A character string specifying the type of the symmetric
positive definite matrix characterizing the elliptical copula. Implemented
structures are "ex" for exchangeable and "un" for unstructured, see package
\code{copula}.}

\item{m}{Length of blocks. Only necessary if the test \code{gofPIOSTn} is
part of \code{tests}.}

\item{delta.J}{Scaling parameter for the matrix of smoothing parameters.
Only necessary if the test \code{gofKernel} is part of \code{tests}.}

\item{nodes.Integration}{Number of knots of the bivariate Gauss-Legendre
quadrature. Only necessary if the test \code{gofKernel} is part of
\code{tests}.}

\item{lower}{Lower bound for the maximum likelihood estimation of the copula
parameter. The constraint is also active in the bootstrapping procedure. The
constraint is not active when a switch to inversion of Kendall's tau is
necessary. Default \code{NULL}.}

\item{upper}{Upper bound for the maximum likelihood estimation of the copula
parameter. The constraint is also active in the bootstrapping procedure. The
constraint is not active when a switch to inversion of Kendall's tau is
necessary. Default \code{NULL}.}

\item{seed.active}{Has to be either an integer or a vector of M+1 integers.
If an integer, then the seeds for the bootstrapping procedure will be
simulated. If M+1 seeds are provided, then these seeds are used in the
bootstrapping procedure. Defaults to \code{NULL}, then \code{R} generates
the seeds from the computer runtime. Controlling the seeds is useful for
reproducibility of a simulation study to compare the power of the tests or
for reproducibility of an empirical study.}

\item{processes}{The number of parallel processes which are performed to
speed up the bootstrapping. Shouldn't be higher than the number of logical
processors. Please see the details.}
}
\value{
A list containing several objects of \code{class} gofCOP with the
following components for each copulae \item{method}{a character which
informs about the performed analysis} \item{copula}{the copula tested for}
\item{margins}{the method used to estimate the margin distribution.}
\item{param.margins}{the parameters of the estimated margin distributions.
Only applicable if the margins were not specified as \code{"ranks"} or
\code{NULL}.} \item{theta}{dependence parameters of the copulae}
\item{df}{the degrees of freedem of the copula. Only applicable for
t-copula.} \item{res.tests}{a matrix with the p-values and test statistics
of the hybrid and the individual tests}
}
\description{
\code{\link{gof}} computes for a given dataset and based on the choices of
the user different tests for different copulae. If copulae are given, all
the implemented tests for those copulae are calculated. If tests are given,
all the implemented copulae for every test are used. If both copulae and
tests are given, all possible combinations are calculated.
}
\details{
If a character vector is given for the argument \code{copula} and nothing
for \code{tests}, then all tests are performed for which the given copulae
are implemented. If \code{tests} contains a character vector of tests and
\code{copula = NULL}, then this tests will be performed for all implemented
copulae. If character vectors are given for \code{copula} and \code{tests},
then the tests are performed with the given copulae. If \code{tests = NULL}
and \code{copula = NULL}, then the argument \code{priority} catches in and
defines the procedure.

For small values of \code{M}, initializing the parallelisation via
\code{processes} does not make sense. The registration of the parallel
processes increases the computation time. Please consider to enable
parallelisation just for high values of \code{M}.

Note that this function does not display \code{warning()} messages. Due to
the large amount of tests run at once, the messages are not tracable to the
situation when they appeared. Hence they are omitted for this function.
}
\examples{

data(IndexReturns2D)

gof(IndexReturns2D, priority = "tests", copula = "normal", 
tests = c("gofRosenblattSnB", "gofRosenblattSnC"), M = 5)

}
