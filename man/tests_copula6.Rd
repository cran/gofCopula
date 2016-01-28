\name{gofRosenblattSnC}
\alias{gofRosenblattSnC}

\title{The SnC test based on the Rosenblatt transformation
}
\description{
\code{\link{gofRosenblattSnC}} is a wrapper for the functions \code{gofCopula}, \code{fitCopula}, \code{ellipCopula} and \code{archmCopula} from the package \pkg{copula}. It combines these functions to test a dataset for a copula directly without all the necessary intermediate steps and extends its functionality. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. \code{\link{gofRosenblattSnC}} contains the SnC gof test from Genest (2009) for copulae and compares the empirical copula against a parametric estimate of the copula derived under the null hypothesis. The approximate p-values are obtained via parametric bootstrapping. It is possible to insert datasets of all dimensions above 1 and the possible copulae are "gaussian", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used.
}
\usage{
gofRosenblattSnC(copula, x, M = 1000, param = 0.5, 
                  param.est = T, df = 4, df.est = T, margins = "ranks", 
                  execute.times.comp = T)
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
\item{margins}{
Specifies which estimation method shall be used in case that the input data are not in the range [0,1]. The default is \code{"ranks"}, which is the standard approach to convert data in such a case. Alternatively can the following distributions be specified: \code{"beta"}, \code{"cauchy"}, Chi-squared (\code{"chisq"}), \code{"f"}, \code{"gamma"}, Log normal (\code{"lnorm"}), Normal (\code{"norm"}), \code{"t"}, \code{"weibull"}, Exponential (\code{"exp"}).
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
}
\details{
This test is based on the Rosenblatt probability integral transform which uses the mapping \eqn{\mathcal{R}: (0,1)^d \rightarrow (0,1)^d}{R : (0,1)^d -> (0,1)^d} to test the \eqn{H_0}{H0} hypothesis
\deqn{C \in \mathcal{C}_0}{C in Ccal0}
with \eqn{\mathcal{C}_0}{Ccal0} as the true class of copulae under \eqn{H_0}{H0}. Following Genest et al. (2009) ensures this transformation the decomposition of a random vector \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} with a distribution into mutually independent elements with a uniform distribution on the unit interval. The mapping provides pseudo observations \eqn{\mathbf{E}_i}{E[i]}, given by \deqn{E_1 = \mathcal{R}(U_1), \dots, E_n = \mathcal{R}(U_n).}{E_1 = R(U_1), ..., E_n = R(U_n).} The mapping is performed by assigning to every vector \eqn{\mathbf{u}}{u} for \eqn{e_1 = u_1}{e[1] = u[1]} and for \eqn{i \in \{2, \dots, d\}}{i in {2, ..., d}},
\deqn{e_i = \frac{\partial^{i-1} C(u_1, \dots, u_i, 1, \dots, 1)}{\partial u_1 \cdots \partial u_{i-1}} / \frac{\partial^{i-1} C(u_1, \dots, u_{i-1}, 1, \dots, 1)}{\partial u_1 \cdots \partial u_{i-1}}.}{e[i] = (d^(i-1) C(u[1], ..., u[i], 1, ..., 1))/(d u[1] ... d u[i-1]) / (d^(i-1) C(u[1], ..., u[i-1], 1, ..., 1))/(d u[1] ... d u[i-1]).} The resulting independence copula is given by \eqn{C_{\bot}(\mathbf{u}) = u_1 \cdot \dots \cdot u_d}{Cind(u) = u[1] x ... x u[d]}.

The test statistic \eqn{T} is then defined as

\deqn{T = n \int_{[0,1]^d} \{ D_n(\mathbf{u}) - C_{\bot}(\mathbf{u}) \}^2 d D_n(\mathbf{u})}{n int_{[0,1]^d} ( {Dn(u) - Cind(u)}^2 )dDn(u)}
with
\eqn{D_n(\mathbf{u}) = \frac{1}{n} \sum_{i = 1}^n \mathbf{I}(\mathbf{E}_i \leq \mathbf{u})}{Dn(u) = 1/n sum(E[i] <= u, i = 1, ..., n)}.

The approximate p-value is computed by the formula, see \pkg{copula},

\deqn{(0.5 + \sum_{b=1}^N \mathbf{I}(T_b \geq T) / (N+1)}{(0.5 + sum(T[b] >= T, b=1, .., N)) / (N+1)}

where \eqn{T} and \eqn{T_b}{T[b]} denote the test statistic and the bootstrapped test statistc, respectively. This ensures that the approximate p-value is a number strictly between 0 and 1, which is sometimes necessary for further treatments. See Pesarin (2001) for more details.
}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-value and test statistic of test}
}
\references{
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}\cr \cr
Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan (2014). copula: Multivariate Dependence with Copulas. \emph{R package version 0.999-12.}. \url{http://CRAN.R-project.org/package=copula} \cr \cr
Pesarin, F. (2001). \emph{Multivariate Permutation Tests: With applications in Biostatistics}, Wiley
}
\examples{
data = cbind(rnorm(100), rnorm(100), rnorm(100))

gofRosenblattSnC("gaussian", data, M = 20)
}
