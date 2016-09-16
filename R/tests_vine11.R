gofKendallCvM = function(copula = c("normal", "t", "clayton", "frank", "gumbel"), x, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", M = 100, execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofKendallCvM.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofKendallCvM", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnK", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnK", estim.method = "itau", processes = processes)} else {res}
}
