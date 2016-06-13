\name{gofRosenblattChisq}
\alias{gofRosenblattChisq}

\title{Gof test using the Anderson-Darling test statistic and the chi-square distribution
}
\description{
\code{\link{gofRosenblattChisq}} contains the RosenblattChisq gof test for copulae, described in Genest (2009) and Hofert (2014), and compares the empirical copula against a parametric estimate of the copula derived under the null hypothesis. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. The approximate p-values are computed with a parametric bootstrap, which computation can be accelerated by enabling in-build parallel computation. The gof statistics are computed with the function \code{\link{gofTstat}} from the package copula. It is possible to insert datasets of all dimensions above 1 and the possible copulae are "normal", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used.
}
\usage{
gofRosenblattChisq(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, 
            df.est = T, margins = "ranks", dispstr = "ex", 
            execute.times.comp = T, processes = 1)
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
Indicates if \code{df} shall be estimated. Has to be either \code{FALSE} or \code{TRUE}, where \code{TRUE} means that it will be estimated.
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
As written in Hofert et al. (2014) computes this Anderson-Darling test statistic (supposedly) U[0,1]-distributed (under H_0) random variates via the distribution function of chi-square distribution with d degrees of freedom. The \eqn{H_0}{H0} hypothesis is 
\deqn{C \in \mathcal{C}_0}{C in Ccal0}
with \eqn{\mathcal{C}_0}{Ccal0} as the true class of copulae under \eqn{H_0}{H0}.

This test is based on the Rosenblatt probability integral transform which uses the mapping \eqn{\mathcal{R}: (0,1)^d \rightarrow (0,1)^d}{R : (0,1)^d -> (0,1)^d}. Following Genest et al. (2009) ensures this transformation the decomposition of a random vector \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} with a distribution into mutually independent elements with a uniform distribution on the unit interval. The mapping provides pseudo observations \eqn{E_i}{E[i]}, given by \deqn{E_1 = \mathcal{R}(U_1), \dots, E_n = \mathcal{R}(U_n).}{E_1 = R(U_1), ..., E_n = R(U_n).} The mapping is performed by assigning to every vector \eqn{\mathbf{u}}{u} for \eqn{e_1 = u_1}{e[1] = u[1]} and for \eqn{i \in \{2, \dots, d\}}{i in {2, ..., d}},
\deqn{e_i = \frac{\partial^{i-1} C(u_1, \dots, u_i, 1, \dots, 1)}{\partial u_1 \cdots \partial u_{i-1}} / \frac{\partial^{i-1} C(u_1, \dots, u_{i-1}, 1, \dots, 1)}{\partial u_1 \cdots \partial u_{i-1}}.}{e[i] = (d^(i-1) C(u[1], ..., u[i], 1, ..., 1))/(d u[1] ... d u[i-1]) / (d^(i-1) C(u[1], ..., u[i-1], 1, ..., 1))/(d u[1] ... d u[i-1]).}

The Anderson-Darling test statistic of the variates

\deqn{G(x_j) = \chi_d^2 \left( x_j \right)}{G(x[j]) = pchisq(x[j], df=d)}

is computed (via \code{ADGofTest::ad.test}), where \eqn{x_j = \sum_{i=1}^d (\Phi^{-1}(e_{ij}))^2}{x[j] = (Phi^{-1}(e_{1j}))^2+...+(Phi^{-1}(e_{dj}))^2}, \eqn{\Phi^{-1}}{Phi^{-1}} denotes the quantile function of the standard normal distribution function, \eqn{\chi_d^2}{pchisq(.,df=d)} denotes the distribution function of the chi-square distribution with \code{d} degrees of freedom, and \eqn{u_{ij}} is the \eqn{j}th component in the \eqn{i}th row of \eqn{\mathbf{u}}{u}.

The test statistic is then given by
\deqn{T = -n - \sum_{j=1}^n \frac{2j - 1}{n} [\ln(G(x_j)) + \ln(1 - G(x_{n+1-j}))].}{T = -n - sum((2j - 1)/n [ln(G(x[j])) + ln(1 - G(x[n+1-j]))], j = 1, ..., n).} 

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
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}\cr \cr
Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan (2014). copula: Multivariate Dependence with Copulas. \emph{R package version 0.999-12.}. \url{http://CRAN.R-project.org/package=copula} \cr \cr
}
\examples{
data(IndexReturns)

gofRosenblattChisq("normal", IndexReturns[c(1:100),c(1:2)], M = 20)
}
