\name{gofKendallCvM}
\alias{gofKendallCvM}

\title{2 dimensional gof test (Cramer-von Mises) based on Kendall's process
}
\description{
\code{\link{gofKendallCvM}} is a wrapper for the function \code{BiCopGofTest} from the package \pkg{VineCopula} which extends its functionality. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. It tests a given 2 dimensional dataset for a copula based on Kendall's process with the Cramer-von Mises test statistic. The possible copulae are "gaussian", "gumbel", "clayton" and "frank". See for reference Genest et al. (2009). The parameter estimation is performed with pseudo maximum likelihood method.
}
\usage{
gofKendallCvM(copula, x, param = 0, param.est = T, margins = "ranks", 
              M = 100, execute.times.comp = T)
}
\arguments{
\item{copula}{
The copula to test for. Possible are the copulae \code{"gaussian"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
}
  \item{x}{
A 2 dimensional matrix containing the residuals of the data.
}
  \item{param}{
The parameter to be used.
}
\item{param.est}{
Shall be either \code{TRUE} or \code{FALSE}. \code{TRUE} means that \code{param} will be estimated with a maximum likelihood estimation.
}
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{M}{
Number of bootstrap samples.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
}
\details{
With the pseudo observations \eqn{U_{ij}}{U[ij]} for \eqn{i = 1, \dots,n}{i = 1, ...,n}, \eqn{j = 1, \dots,d}{j = 1, ...,d} and \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} is the empirical copula given by \eqn{C_n(\mathbf{u}) = \frac{1}{n} \sum_{i = 1}^n \mathbf{I}(U_{i1} \leq u_1, \dots, U_{id} \leq u_d).}{1/n sum(U[i1] <= u_1, ..., U[id] <= u_d, i = 1, ..., n).} Let the rescaled pseudo observations be \eqn{\mathbf{V}_1 = C_n(\mathbf{U}_1), \dots, \mathbf{V}_n = C_n(\mathbf{U}_n)}{V[1] = Cn(U[1]), ..., V[n] = Cn(U[n]} and the distribution function of \eqn{\mathbf{V}}{V} shall be \eqn{K}. The estimated version is given by
\deqn{K_n(v) = \frac{1}{n} \sum_{i=1}^n \mathbf{I}(\mathbf{V}_i \leq v)}{Kn(v) = 1/n sum(I(V[i] <= v, i = 1, ..., n))}
with \eqn{v \in [0,1]^d.}{v in [0,1]^d.}
The testable \eqn{H_0^{'}}{H0^'} hypothesis is then
\deqn{K \in \mathcal{K}_0 = \{K_{\theta} : \theta \in \Theta \}}{K in K_0 = {Ktheta : theta in Theta}}
with \eqn{\Theta}{Theta} being an open subset of \eqn{R^p}{R^p} for an integer \eqn{p \geq 1}{p >= 1}, see Genest et al. (2009). The resulting Cramer-von Mises test statistic is then given by
\deqn{T = n \int_0^1 (K_n(v) - K_{\theta_n})^2 d K_{\theta_n}(v).}{T = n int_0^1 (Kn(v) - Kthetan)^2 d Kthetan(v).}

Because \eqn{H_0^{'}}{H0^'} consists of more distributions than the \eqn{H_0}{H0} is the test not necessarily consistent.

The approximate p-value is computed by the formula

\deqn{\frac{1}{M} \sum_{b=1}^M \mathbf{I}_{\{T_b \leq T\}}.}{1/M sum(T[b] >= T, b=1, ...,M).}

}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-value and test statistic of test}
}
\references{
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}\cr \cr
Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler (2015). VineCopula: Statistical Inference of Vine Copulas. \emph{R package version 1.4.}. \url{https://cran.r-project.org/package=VineCopula}
}
\examples{
data = cbind(rnorm(100), rnorm(100))

gofKendallCvM("gaussian", data, M = 1)
}
