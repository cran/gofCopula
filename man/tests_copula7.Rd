\name{gofADGamma}
\alias{gofADGamma}

\title{Gof test using the Anderson-Darling test statistic and the gamma distribution
}
\description{
\code{\link{gofADGamma}} is a wrapper for the functions \code{gofCopula}, \code{fitCopula}, \code{ellipCopula} and \code{archmCopula} from the package \pkg{copula}. It combines these functions to test a dataset for a copula directly without all the necessary intermediate steps and extends its functionality. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. \code{\link{gofADGamma}} contains the ADGamma gof tests for copulae, described in Genest (2009) and Hofert (2014), and compares the empirical copula against a parametric estimate of the copula derived under the null hypothesis. The approximate p-values are computed with a parametric bootstrap. It is possible to insert datasets of all dimensions above 1 and the possible copulae are "gaussian", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used.
}
\usage{
gofADGamma(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, 
            margins = "ranks", execute.times.comp = T)
}
\arguments{
\item{copula}{
The copula to test for. Possible are \code{"gaussian"}, \code{"t"}, \code{"clayton"} and \code{"gumbel"}.
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
As written in Hofert et al. (2014) computes this Anderson-Darling test statistic for (supposedly) U[0,1]-distributed (under H_0) random variates via the distribution function of the gamma distribution. The \eqn{H_0}{H0} hypothesis is 
\deqn{C \in \mathcal{C}_0}{C in Ccal0}
with \eqn{\mathcal{C}_0}{Ccal0} as the true class of copulae under \eqn{H_0}{H0}. 
The Anderson-Darling test statistic of the variates

\deqn{G(x_i) = \Gamma_d \left( x_i \right)}{G(x[i]) = pgamma(x[i], shape=d)}

is computed (via \code{ADGofTest::ad.test}), where \eqn{x_i = \sum_{j=1}^d (- \ln u_{ij})}{x[i] = -log(u_{i1})-...-log(u_{id})}, \eqn{\Gamma_d()}{pgamma( . ,shape=d)} denotes the distribution function of the gamma distribution with shape parameter d and shape parameter one (being equal to an Erlang(d) distribution function).

The test statistic is then given by
\deqn{T = -n - \sum_{i=1}^n \frac{2i - 1}{n} [\ln(G(x_i)) + \ln(1 - G(x_{n+1-i}))].}{T = -n - sum((2i - 1)/n [ln(G(x[i])) + ln(1 - G(x[n+1-i]))], i = 1, ..., n).} 

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
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}\cr \cr
Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan (2014). copula: Multivariate Dependence with Copulas. \emph{R package version 0.999-12.}. \url{http://CRAN.R-project.org/package=copula} \cr \cr
Pesarin, F. (2001). \emph{Multivariate Permutation Tests: With applications in Biostatistics}, Wiley
}
\examples{
data = cbind(rnorm(100), rnorm(100), rnorm(100))

gofADGamma("gaussian", data, M = 20)
}
