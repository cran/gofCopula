\name{gofWhite}
\alias{gofWhite}

\title{2 dimensional gof tests based on White's information matrix equality.
}
\description{
\code{\link{gofWhite}} tests a given 2 dimensional dataset for a copula  with the gof test based on White's information matrix equality. The possible copulae are "normal", "t", "gumbel", "clayton" and "frank". See for reference Schepsmeier et al. (2015). The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. The approximate p-values are computed with a parametric bootstrap, which computation can be accelerated by enabling in-build parallel computation. The computation of the test statistic and p-values is performed by corresponding functions from the \code{\link{VineCopula}} package.
}
\usage{
gofWhite(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, 
          margins = "ranks", execute.times.comp = T, processes = 1)
}
\arguments{
\item{copula}{
The copula to test for. Possible are the copulae \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
}
  \item{x}{
A 2 dimensional matrix containing the residuals of the data.
}
\item{M}{
Number of bootstrap samples.
}
  \item{param}{
The parameter to be used.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated with a maximum likelihood estimation.
}
  \item{df}{
The degrees of freedom for the \code{"t"}-copula.
}
\item{df.est}{
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors. Please see the details.}
}
\details{
The details are obtained from Schepsmeier et al. (2015) who states that this test uses the information matrix equality of White (1982). Under correct model specification is the Fisher Information equivalently calculated as minus the expected Hessian matrix or as the expected outer product of the score function. The null hypothesis is 
\deqn{H_0 : \mathbf{H}(\theta) + \mathbf{S}(\theta) = 0}{H0 : H(theta) + S(theta) = 0}
where \eqn{\mathbf{H}(\theta)}{H(theta)} is the expected Hessian matrix and \eqn{\mathbf{S}(\theta)}{S(theta)} is the expected outer product of the score function.

The test statistic is derived by
\deqn{T_n = n(\bar{d}(\theta_n))^\top V_{\theta_n}^{-1} \bar{d}(\theta_n)}{T_n = n(dbar(theta_n))^T V_(theta_n)^(-1) dbar(theta_n)}
with
\deqn{\bar{d}(\theta_n) = \frac{1}{n} \sum_{i=1}^n vech(\mathbf{H}_n(\theta_n|\mathbf{u}) + \mathbf{S}_n(\theta_n|\mathbf{u})),}{dbar(theta_n) = 1/n sum(vech(H(theta_n|u) + S(theta_n|u)), i=1, ..., n),}

\deqn{d(\theta_n) = vech(\mathbf{H}_n(\theta_n|\mathbf{u}) + \mathbf{S}_n(\theta_n|\mathbf{u})),}{d(theta_n) = vech(H(theta_n|u) + S(theta_n|u)),}

\deqn{V_{\theta_n} = \frac{1}{n} \sum_{i=1}^n (d(\theta_n) - D_{\theta_n} \mathbf{H}_n(\theta_n)^{-1} \delta l(\theta_n))(d(\theta_n) - D_{\theta_n} \mathbf{H}_n(\theta_n)^{-1} \delta l(\theta_n))^\top}{Vtheta_n] = 1/n sum(d(theta_n) - D(theta_n) H_n(theta_n)^(-1) d l(theta_n)d(theta_n) - D(theta_n) H_n(theta_n)^(-1) d l(theta_n)^T, i=1, ..., n)}
and
\deqn{D_{\theta_n} = \frac{1}{n} \sum_{i=1}^n [\delta_{\theta_k} d_l(\theta_n)]_{l=1, \dots, \frac{p(p+1)}{2}, k=1, \dots, p}}{D(theta_n = 1/n sum(d_(theta_k) d[l](theta_n)_(l=1, ...,(p(p+1))/2, k=1, ..., p), i=1, ..., n)}
where \eqn{l(\theta_n)}{l(theta_n)} represents the log likelihood function and \eqn{p} is the length of the parameter vector \eqn{\theta}{theta}.

The test statistic will be rejected if
\deqn{T > (1 - \alpha) (\chi^2_{p(p+1)/2})^{-1}.}{T > (1 - alpha) (chi^2_(p(p+1)/2))^(-1).}

For small values of \code{M}, initializing the parallization via \code{processes} does not make sense. The registration of the parallel processes increases the computation time. Please consider to enable parallelization just for high values of \code{M}.
}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-value and test statistic of test}
}
\references{
Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler (2015). VineCopula: Statistical Inference of Vine Copulas. \emph{R package version 1.4.}. \url{https://cran.r-project.org/package=VineCopula} \cr \cr
Schepsmeier, U. and J. Stoeber (2014). Derivatives and Fisher information of bivariate copulas. \emph{Statistical Papers, 55(2), 525-542.} \url{http://link.springer.com/article/10.1007/s00362-013-0498-x} \cr \cr
Stoeber, J. and U. Schepsmeier (2013). Estimating standard errors in regular vine copula models \emph{Computational Statistics, 28 (6), 2679-2707} \cr \cr
Schepsmeier, U. (2015). Efficient information based goodness-of-fit tests for vine copula models with fixed margins. \emph{Journal of Multivariate Analysis 138, 34-52.}
Schepsmeier, U. (2014). A goodness-of-fit test for regular vine copula models.
}
\examples{
data(IndexReturns)

gofWhite("normal", IndexReturns[c(1:100),c(1:2)], M = 10)
}

