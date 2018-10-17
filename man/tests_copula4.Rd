\name{gofSn}
\alias{gofSn}

\title{The Sn gof test using the empirical copula
}
\description{
\code{\link{gofSn}} performs the \code{"Sn"} gof test, described in Genest et al. (2009), for copulae and compares the empirical copula against a parametric estimate of the copula derived under the null hypothesis. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. The approximate p-values are computed with a parametric bootstrap, which computation can be accelerated by enabling in-build parallel computation. The gof statistics are computed with the function \code{\link{gofTstat}} from the package copula. It is possible to insert datasets of all dimensions above 1 and the possible copulae are "normal", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used.
}
\usage{
gofSn(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, 
      margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1)
}
\arguments{
\item{copula}{
The copula to test for. Possible are \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
}
  \item{x}{
A matrix containing the residuals of the data.
}
  \item{M}{
Number of bootstrapping loops.
}
\item{param}{
The copula parameter to use, if it shall not be estimated.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated.
}
\item{df}{
Degrees of freedom, if not meant to be estimated. Only necessary if tested for \code{"t"}-copula.
}
\item{df.est}{
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where\code{TRUE} means that it will be estimated.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{dispstr}{
A character string specifying the type of the symmetric positive definite matrix characterizing the elliptical copula. Implemented structures are "ex" for exchangeable and "un" for unstructured, see package \code{copula}.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors. Please see the details.}
}
\details{
With the pseudo observations \eqn{U_{ij}}{U[ij]} for \eqn{i = 1, \dots,n}{i = 1, ...,n}, \eqn{j = 1, \dots,d}{j = 1, ...,d} and \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} is the empirical copula given by \eqn{C_n(\mathbf{u}) = \frac{1}{n} \sum_{i = 1}^n \mathbf{I}(U_{i1} \leq u_1, \dots, U_{id} \leq u_d).}{1/n sum(U[i1] <= u_1, ..., U[id] <= u_d, i = 1, ..., n).} It shall be tested the \eqn{H_0}{H0} hypothesis:
\deqn{C \in \mathcal{C}_0}{C in Ccal0}
with \eqn{\mathcal{C}_0}{Ccal0} as the true class of copulae under \eqn{H_0}{H0}.
The test statistic \eqn{T} is then defined as

\deqn{T = n \int_{[0,1]^d} \{ C_n(\mathbf{u}) - C_{\theta_n}(\mathbf{u}) \}^2 d C_n(\mathbf{u})}{n int_{[0,1]^d} ( {Cn(u) - Cthetan(u)}^2 )dCn(u)}
with
\eqn{C_{\theta_n}(\mathbf{u})}{Cthetan(u)} the estimation of \eqn{C} under the \eqn{H_0}{H0}.

The approximate p-value is computed by the formula,

\deqn{\sum_{b=1}^M \mathbf{I}(|T_b| \geq |T|) / M,}{sum(|T[b]| >= |T|, b=1, .., M) / M,}

where \eqn{T} and \eqn{T_b}{T[b]} denote the test statistic and the bootstrapped test statistc, respectively. 

For small values of \code{M}, initializing the parallization via \code{processes} does not make sense. The registration of the parallel processes increases the computation time. Please consider to enable parallelization just for high values of \code{M}.
}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-value and test statistic of test}
}
\references{
Rosenblatt, M. (1952). Remarks on a Multivariate Transformation. \emph{The Annals of Mathematical Statistics 23, 3, 470-472}.\cr \cr
Hering, C. and Hofert, M. (2014). Goodness-of-fit tests for Archimedean copulas in high dimensions. \emph{Innovations in Quantitative Risk Management}.\cr \cr
Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan (2014). copula: Multivariate Dependence with Copulas. \emph{R package version 0.999-15.}. \url{https://cran.r-project.org/package=copula} \cr \cr
}
\examples{
data(IndexReturns)

gofSn("normal", IndexReturns[c(1:100),c(1:2)], M = 10)
}