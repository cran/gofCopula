% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests_White.R
\name{gofWhite}
\alias{gofWhite}
\title{2 dimensional gof tests based on White's information matrix equality.}
\usage{
gofWhite(
  copula = c("normal", "t", "clayton", "gumbel", "frank", "joe"),
  x,
  param = 0.5,
  param.est = TRUE,
  df = 4,
  df.est = TRUE,
  margins = "ranks",
  flip = 0,
  M = 1000,
  lower = NULL,
  upper = NULL,
  seed.active = NULL,
  processes = 1
)
}
\arguments{
\item{copula}{The copula to test for. Possible are the copulae
\code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"}, 
\code{"frank"} and \code{"joe"}.}

\item{x}{A matrix containing the data with rows being observations and
columns being variables.}

\item{param}{The parameter to be used.}

\item{param.est}{Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE}
means that \code{param} will be estimated with a maximum likelihood
estimation.}

\item{df}{The degrees of freedom, if not meant to be estimated. Only
necessary if tested for \code{"t"}-copula.}

\item{df.est}{Indicates if \code{df} shall be estimated. Has to be either
\code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be
estimated.}

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

\item{flip}{The control parameter to flip the copula by 90, 180, 270 degrees
clockwise. Only applicable for bivariate copula. Default is 0 and possible 
inputs are 0, 90, 180, 270 and NULL.}

\item{M}{Number of bootstrap samples.}

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
An object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{copula}{the copula tested for} \item{margins}{the method used to
estimate the margin distribution.} \item{param.margins}{the parameters of
the estimated margin distributions. Only applicable if the margins were not
specified as \code{"ranks"} or \code{NULL}.} \item{theta}{dependence
parameters of the copulae} \item{df}{the degrees of freedem of the copula.
Only applicable for t-copula.} \item{res.tests}{a matrix with the p-values
and test statistics of the hybrid and the individual tests}
}
\description{
\code{\link{gofWhite}} tests a given 2 dimensional dataset for a copula with
the gof test based on White's information matrix equality. The possible
copulae are \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"}, 
\code{"frank"} and \code{"joe"}. See for reference Schepsmeier et al. (2015). 
The parameter estimation is performed with pseudo maximum likelihood method. 
In case the estimation fails, inversion of Kendall's tau is used. The margins
can be estimated by a bunch of distributions and the time which is necessary
for the estimation can be given. The approximate p-values are computed with
a parametric bootstrap, which computation can be accelerated by enabling
in-build parallel computation. The computation of the test statistic and
p-values is performed by corresponding functions from the \code{VineCopula}
package.
}
\details{
The details are obtained from Schepsmeier et al. (2015) who states that this
test uses the information matrix equality of White (1982). Under correct
model specification is the Fisher Information equivalently calculated as
minus the expected Hessian matrix or as the expected outer product of the
score function. The null hypothesis is \deqn{H_0 : \mathbf{V}(\theta) +
\mathbf{S}(\theta) = 0}{H0 : V(theta) + S(theta) = 0} where
\eqn{\mathbf{V}(\theta)}{V(theta)} is the expected Hessian matrix and
\eqn{\mathbf{S}(\theta)}{S(theta)} is the expected outer product of the
score function.

The test statistic is derived by \deqn{T_n = n(\bar{d}(\theta_n))^\top
A_{\theta_n}^{-1} \bar{d}(\theta_n)}{T_n = n(dbar(theta_n))^T
A_(theta_n)^(-1) dbar(theta_n)} with \deqn{\bar{d}(\theta_n) = \frac{1}{n}
\sum_{i=1}^n vech(\mathbf{V}_n(\theta_n|\mathbf{u}) +
\mathbf{S}_n(\theta_n|\mathbf{u})),}{dbar(theta_n) = 1/n
sum(vech(V(theta_n|u) + S(theta_n|u)), i=1, ..., n),}

\deqn{d(\theta_n) = vech(\mathbf{V}_n(\theta_n|\mathbf{u}) +
\mathbf{S}_n(\theta_n|\mathbf{u})),}{d(theta_n) = vech(V(theta_n|u) +
S(theta_n|u)),}

\deqn{A_{\theta_n} = \frac{1}{n} \sum_{i=1}^n (d(\theta_n) - D_{\theta_n}
\mathbf{V}_n(\theta_n)^{-1} \delta l(\theta_n))(d(\theta_n) - D_{\theta_n}
\mathbf{V}_n(\theta_n)^{-1} \delta l(\theta_n))^\top}{Atheta_n] = 1/n
sum(d(theta_n) - D(theta_n) V_n(theta_n)^(-1) d l(theta_n)d(theta_n) -
D(theta_n) V_n(theta_n)^(-1) d l(theta_n)^T, i=1, ..., n)} and
\deqn{D_{\theta_n} = \frac{1}{n} \sum_{i=1}^n [\delta_{\theta_k}
d_l(\theta_n)]_{l=1, \dots, \frac{p(p+1)}{2}, k=1, \dots, p}}{D(theta_n =
1/n sum(d_(theta_k) d[l](theta_n)_(l=1, ...,(p(p+1))/2, k=1, ..., p), i=1,
..., n)} where \eqn{l(\theta_n)}{l(theta_n)} represents the log likelihood
function and \eqn{p} is the length of the parameter vector
\eqn{\theta}{theta}.

The test statistic will be rejected if \deqn{T > (1 - \alpha)
(\chi^2_{p(p+1)/2})^{-1}.}{T > (1 - alpha) (chi^2_(p(p+1)/2))^(-1).}

For small values of \code{M}, initializing the parallelisation via
\code{processes} does not make sense. The registration of the parallel
processes increases the computation time. Please consider to enable
parallelisation just for high values of \code{M}.

Please note, the test gofWhite may be unstable for t-copula. Please handle
the results carefully.
}
\examples{

data(IndexReturns2D)

gofWhite("normal", IndexReturns2D, M = 10)

}
\references{
Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann,
Benedikt Graeler (2015). VineCopula: Statistical Inference of Vine Copulas.
\emph{R package version 1.4.}.
\url{https://cran.r-project.org/package=VineCopula} \cr \cr Schepsmeier, U.
and J. Stoeber (2014). Derivatives and Fisher information of bivariate
copulas. \emph{Statistical Papers, 55(2), 525-542.}
\url{https://link.springer.com/article/10.1007/s00362-013-0498-x} \cr \cr
Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular
vine copula models \emph{Computational Statistics, 28 (6), 2679-2707} \cr
\cr Schepsmeier, U. (2015). Efficient information based goodness-of-fit
tests for vine copula models with fixed margins. \emph{Journal of
Multivariate Analysis 138, 34-52.} Schepsmeier, U. (2014). A goodness-of-fit
test for regular vine copula models.
}
