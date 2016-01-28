\name{gofRn}
\alias{gofRn}

\title{The Rn gof test from package \pkg{copula}
}
\description{
\code{\link{gofRn}} is a wrapper for the functions \code{gofCopula}, \code{fitCopula}, \code{ellipCopula} and \code{archmCopula} from the package \pkg{copula}. It combines these functions to test a dataset for a copula directly without all the necessary intermediate steps. \code{\link{gofRn}} performs the gof test from Genest et al. (2013) for copulae and compares the empirical copula against a parametric estimate of the copula derived under the null hypothesis. The approximate p-values are computed with a fast multiplier approach. It is just possible to insert datasets of dimension 2 and the possible copulae are "gaussian", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used.
}
\usage{
gofRn(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T,
      m_b = 0.5, zeta.m = 0, b_Rn = 0.05, execute.times.comp = T)
}
\arguments{
\item{copula}{
The copula to test for. Possible are \code{"gaussian"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
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
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated.
}
\item{m_b}{
The power of the statistic.
}
\item{zeta.m}{
The adjustment parameter.
}
\item{b_Rn}{
The bandwidth for the estimation of the first-order partial derivatives based on the empirical copula.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
}
\details{
With the pseudo observations \eqn{U_{ij}}{U[ij]} for \eqn{i = 1, \dots,n}{i = 1, ...,n}; \eqn{j = 1, \dots,d}{j = 1, ...,d} and \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} is the empirical copula given by \eqn{C_n(\mathbf{u}) = \frac{1}{n} \sum_{i = 1}^n \mathbf{I}(U_{i1} \leq u_1, \dots, U_{id} \leq u_d).}{1/n sum(U[i1] <= u_1, ..., U[id] <= u_d, i = 1, ..., n).} It shall be tested the \eqn{H_0}{H0} hypothesis:
\deqn{C \in \mathcal{C}_0}{C in Ccal0}
with \eqn{\mathcal{C}_0}{Ccal0} as the true class of copulae under \eqn{H_0}{H0}.
The test statistic \eqn{T} is defined as

\deqn{T = n \int_{[0,1]^d} \{ \frac{C_n(\mathbf{u}) - C_{\theta_n}(\mathbf{u})}{[C_{\theta_n}(\mathbf{u}) \{1 - C_{\theta_n}(\mathbf{u})\} + \zeta_m]^m} \}^2 d C_n(\mathbf{u})}{n int_{[0,1]^d} [{Cn(u) - Cthetan(u)}/(Cthetan(u){1 - Cthetan(u)} + zetam)^m]^2 dCn(u)}
with 
\eqn{C_{\theta_n}(\mathbf{u})}{Cthetan(u)} an estimation of \eqn{C} under the \eqn{H_0}{H0}. \eqn{m \geq 0}{m >= 0} and \eqn{\zeta_m \geq 0}{zetam >= 0} are tuning parameters.

The approximate p-value is computed by the formula, see \pkg{copula},

\deqn{(0.5 + \sum_{b=1}^N \mathbf{I}(T_b \geq T) / (N+1),}{(0.5 + sum(T[b] >= T, b=1, .., N)) / (N+1),}

where \eqn{T} and \eqn{T_b}{T[b]} denote the test statistic and the bootstrapped test statistc, respectively. This ensures that the approximate p-value is a number strictly between 0 and 1, which is sometimes necessary for further treatments. See Pesarin (2001) for more details.

}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-value and test statistic of test}
}
\references{
Christian Genest, Wanling Huang and Jean-Marie Dufour (2013). A regularized goodness-of-fit test for copulas. \emph{Journal de la Societe Francaise de Statistique et revue de statistique appliquee 154.1: 64-77}.\cr \cr
Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan (2014). copula: Multivariate Dependence with Copulas. \emph{R package version 0.999-12.}. \url{http://CRAN.R-project.org/package=copula} \cr \cr
Pesarin, F. (2001). \emph{Multivariate Permutation Tests: With applications in Biostatistics}, Wiley
}
\examples{
data = cbind(rnorm(100), rnorm(100))

gofRn("gaussian", data, M = 20)
}
