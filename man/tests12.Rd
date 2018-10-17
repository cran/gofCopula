\name{gofPIOSTn}
\alias{gofPIOSTn}

\title{2 and 3 dimensional gof test based on the in-and-out-of-sample approach
}
\description{
\code{gofPIOSTn} tests a 2 or 3 dimensional dataset with the PIOS test for a copula. The possible copulae are "normal", "t", "gumbel", "clayton" and "frank". The parameter estimation is performed with pseudo maximum likelihood method. In case the estimation fails, inversion of Kendall's tau is used. The approximate p-values are computed with a semiparametric bootstrap, which computation can be accelerated by enabling in-build parallel computation.
}
\usage{
gofPIOSTn(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, 
          margins = "ranks", dispstr = "ex", m = 1, 
          execute.times.comp = T, processes = 1)
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
\item{m}{
Length of blocks.
}
\item{execute.times.comp}{
Logical. Defines if the time which the estimation most likely takes shall be computed. It'll be just given if \code{M} is at least 100.
}
\item{processes}{
The number of parallel processes which are performed to speed up the bootstrapping. Shouldn't be higher than the number of logical processors. Please see the details.}
}
\details{
The "Tn" test is introduced in Zhang et al. (2015). It tests the \eqn{H_0}{H0} hypothesis
\deqn{H_0 : C_0 \in \mathcal{C}.}{H0 : C0 in Ccal.}
For the test are constructed blocks of length \code{m} out of the data. The test compares then the pseudo likelihood of the data in each block with the overall parameter and with the parameter by leaving out the data in the block. By this procedure can be determined if the data in the block influence the parameter estimation significantly. The test statistic is defined as
\deqn{T = \sum_{b=1}^M \sum_{k=1}^m [l\{U_k^b;\theta_n \} - l\{U_k^b;\theta_n^{-b} \}]}{T = sum(sum(l(U_k^b;theta_n ) - l(U_k^b;theta_n^(-b) ), k=1, ...,m ), b=1, ...,B)}

with the pseudo observations \eqn{U_{ij}}{U[ij]} for \eqn{i = 1, \dots,n}{i = 1, ...,n}; \eqn{j = 1, \dots,d}{j = 1, ...,d} and
\deqn{\theta_n = \arg \min_{\theta} \sum_{i=1}^n l(U_i; \theta)}{theta_n = arg max_theta sum(l(U_i; theta), i=1, ..., n)}
and
\deqn{\theta_n^{-b} = \arg \min_{\theta} \sum_{b^{'} \neq b}^M \sum_{i=1}^m l(U_i^{b^{'}}; \theta), b=1, \dots, M.}{theta_n^(-b) = arg max_theta sum(sum(l(U_i^(b^'); theta), i=1, ..., m), b^'=1, ..., M, b^' != b), b = 1, ..., M.}

The approximate p-value is computed by the formula
\deqn{\sum_{b=1}^M \mathbf{I}(|T_b| \geq |T|) / M,}{sum(|T[b]| >= |T|, b=1, .., M) / M,}

The applied estimation method is the two-step pseudo maximum likelihood approach, see Genest and Rivest (1995).

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
Genest, C., K. G. and Rivest, L.-P. (1995). A semiparametric estimation procedure of dependence
parameters in multivariate families of distributions. \emph{Biometrika, 82:534-552}
}
\examples{
data(IndexReturns)

gofPIOSTn("normal", IndexReturns[c(1:100),c(1:2)], M = 10)
}
