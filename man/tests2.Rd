\name{gofKernel}
\alias{gofKernel}

\title{2 and 3 dimensional gof test of Scaillet
}
\description{
\code{gofKernel} tests a 2 or 3 dimensional dataset with the Scaillet test for a copula. The possible copulae are "normal", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used. The approximate p-values are computed with a parametric bootstrap, which computation can be accelerated by enabling in-build parallel computation.
}
\usage{
gofKernel(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, 
          margins = "ranks", dispstr = "ex", MJ = 100, delta.J = 0.5, 
          nodes.Integration = 12, execute.times.comp = T, processes = 1)
}
\arguments{
\item{copula}{
The copula to test for. Possible are the copulae \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
}
 \item{x}{
A 2 or 3 dimensional matrix containing the residuals of the data.
}
  \item{M}{
Number of bootstrapping loops.
}
  \item{param}{
The parameter to be used.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated with a maximum likelihood estimation.
}
\item{df}{
Degrees of freedom, if not meant to be estimated. Only necessary if tested for \code{"t"}-copula.
}
\item{df.est}{
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{dispstr}{
A character string specifying the type of the symmetric positive definite matrix characterizing the elliptical copula. Implemented structures are "ex" for exchangeable and "un" for unstructured, see package \code{copula}.
}
\item{MJ}{
Size of bootstrapping sample.
}
\item{delta.J}{
Scaling parameter for the matrix of smoothing parameters.
}
\item{nodes.Integration}{
Number of knots of the bivariate Gauss-Legendre quadrature.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors. Please see the details.}
}
\details{
The Scaillet test is a kernel-based goodness-of-fit test with a fixed smoothing parameter. For the copula density \eqn{c(\mathbf{u}, \theta)}{c(u, theta)}, the corresponding kernel estimator is given by
\deqn{c_n(\mathbf{u}) = \frac{1}{n} \sum_{i=1}^n K_H[\mathbf{u} - (U_{i1}, \dots, U_{id})^{\top}], }
where \eqn{U_{ij}}{U[ij]} for \eqn{i = 1, \dots,n}{i = 1, ...,n}; \eqn{j = 1, \dots,d}{j = 1, ...,d} are the pseudo observations, \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} and \eqn{K_H(y) = K(H^{-1}y)/\det(H)}{KH(y) = K(H^(-1)y)/det(H)} for which a bivariate quadratic kernel is used, as in Scaillet (2007). The matrix of smoothing parameters is \eqn{H = 2.6073n^{-1/6} \hat{\Sigma}^{1/2}}{H = 2.6073n^{-1/6} {Sigma_hat}^{1/2}} with \eqn{\hat{\Sigma}}{Sigma_hat} the sample covariance matrix. 
The test statistic is then given by
\deqn{T = \int_{[0,1]^d} \{c_n(\mathbf{u}) - K_H * c(\mathbf{u}, \theta_n)\} \omega(\mathbf{u}) d \mathbf{u}, }{int_([0,1]^d) {c_n(u) - K_H * c(u, theta_n)} omega(u) d u,}
where \eqn{*}{*} denotes the convolution operator and \eqn{\omega}{omega} is a weight function, see Zhang et al. (2015).
The bivariate Gauss-Legendre quadrature method is used to compute the integral in the test statistic numerically, see Scaillet (2007).

The approximate p-value is computed by the formula
\deqn{\sum_{b=1}^M \mathbf{I}(|T_b| \geq |T|) / M,}{sum(|T[b]| >= |T|, b=1, .., M) / M,}

For small values of \code{M}, initializing the parallization via \code{processes} does not make sense. The registration of the parallel processes increases the computation time. Please consider to enable parallelization just for high values of \code{M}.
}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{statistic}{value of the test statistic}
\item{p.value}{the approximate p-value}
}
\references{
Zhang, S., Okhrin, O., Zhou, Q., and Song, P.. Goodness-of-fit Test For Specification of Semiparametric Copula Dependence Models. \emph{under revision in Journal of Econometrics from 15.01.2014} \url{http://sfb649.wiwi.hu-berlin.de/papers/pdf/SFB649DP2013-041.pdf} \cr \cr
Scaillet, O. (2007). Kernel based goodness-of-fit tests for copulas with fixed smoothing parameters.
\emph{Journal of Multivariate Analysis, 98:533-543}
}
\examples{
data(IndexReturns)

gofKernel("normal", IndexReturns[c(1:100),c(1:2)], M = 5, MJ = 5)
}
