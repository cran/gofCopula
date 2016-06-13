\name{gofKendallKS}
\alias{gofKendallKS}

\title{gof test (Kolmogorov-Smirnof) based on Kendall's process
}
\description{
\code{\link{gofKendallKS}} tests a given dataset for a copula based on Kendall's process with the Kolmogorov-Smirnof test statistic. The margins can be estimated by a bunch of distributions and the time which is necessary for the estimation can be given. The possible copulae are "normal", "t", "gumbel", "clayton" and "frank". See for reference Genest et al. (2009). The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used. The approximate p-values are computed with a parametric bootstrap, which computation can be accelerated by enabling in-build parallel computation.
}
\usage{
gofKendallKS(copula, x, param = 0.5, param.est = T, df = 4, df.est = T, 
             margins = "ranks", dispstr = "ex", M = 100, 
             execute.times.comp = T, processes = 1)
}
\arguments{
\item{copula}{
The copula to test for. Possible are the copulae \code{"normal"}, \code{"t"}, \code{"clayton"}, \code{"gumbel"} and \code{"frank"}.
}
  \item{x}{
A matrix containing the residuals of the data.
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
\item{M}{
Number of bootstrap samples.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors. Please see the details.}
}
\details{
With the pseudo observations \eqn{U_{ij}}{U[ij]} for \eqn{i = 1, \dots,n}{i = 1, ...,n}, \eqn{j = 1, \dots,d}{j = 1, ...,d} and \eqn{\mathbf{u} \in [0,1]^d}{u in [0,1]^d} is the empirical copula given by \eqn{C_n(\mathbf{u}) = \frac{1}{n} \sum_{i = 1}^n \mathbf{I}(U_{i1} \leq u_1, \dots, U_{id} \leq u_d).}{1/n sum(U[i1] <= u_1, ..., U[id] <= u_d, i = 1, ..., n).} Let the rescaled pseudo observations be \eqn{V_1 = C_n(U_1), \dots, V_n = C_n(U_n)}{V[1] = Cn(U[1]), ..., V[n] = Cn(U[n]} and the distribution function of \eqn{V}{V} shall be \eqn{K}. The estimated version is given by
\deqn{K_n(v) = \frac{1}{n} \sum_{i=1}^n \mathbf{I}(V_i \leq v)}{Kn(v) = 1/n sum(I(V[i] <= v, i = 1, ..., n))}
with \eqn{v \in [0,1]^d.}{v in [0,1]^d.}
The testable \eqn{H_0^{'}}{H0^'} hypothesis is then
\deqn{K \in \mathcal{K}_0 = \{K_{\theta} : \theta \in \Theta \}}{K in K_0 = {Ktheta : theta in Theta}}
with \eqn{\Theta}{Theta} being an open subset of \eqn{R^p}{R^p} for an integer \eqn{p \geq 1}{p >= 1}, see Genest et al. (2009). The resulting Cramer-von Mises test statistic is then given by
\deqn{T = \sqrt{n} \sup_{v \in [0,1]} |K_n(v) - K_{\theta_n}| .}{T = n^{1/2} sup_{v in [0,1]} |Kn(v) - Kthetan|.}

Because \eqn{H_0^{'}}{H0^'} consists of more distributions than the \eqn{H_0}{H0} is the test not necessarily consistent.

The approximate p-value is computed by the formula

\deqn{\sum_{b=1}^M \mathbf{I}(|T_b| \geq |T|) / M,}{sum(|T[b]| >= |T|, b=1, .., M) / M,}

For small values of \code{M}, initializing the parallization via \code{processes} does not make sense. The registration of the parallel processes increases the computation time. Please consider to enable parallelization just for high values of \code{M}.

}
\value{
A object of the \code{class} gofCOP with the components
\item{method}{a character which informs about the performed analysis}
\item{erg.tests}{a matrix with the p-value and test statistic of test}
}
\references{
Christian Genest, Bruno Remillard, David Beaudoin (2009). Goodness-of-fit tests for copulas: A review and a power study. \emph{Insurance: Mathematics and Economics, Volume 44, Issue 2, April 2009, Pages 199-213, ISSN 0167-6687}. \url{http://dx.doi.org/10.1016/j.insmatheco.2007.10.005}\cr \cr
Christian Genest, Jean-Francois Quessy, Bruno Remillard (2002). Tests of serial independence based on Kendall's process. \emph{The Canadian Journal of Statistics, Volume 30, Issue 3, Sep. 2002, Pages 441-461}. \url{http://www.jstor.org/stable/3316147}\cr \cr
cp. Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler (2015). VineCopula: Statistical Inference of Vine Copulas. \emph{R package version 1.4.}. \url{https://cran.r-project.org/package=VineCopula}
}
\examples{
data(IndexReturns)

gofKendallKS("normal", IndexReturns[c(1:100),c(1:2)], M = 10)
}
