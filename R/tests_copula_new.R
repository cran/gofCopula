########################################################## SnB
gofRosenblattSnB = function (copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRosenblattSnB.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattSnB", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnB", estim.method = "mpl", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnB", estim.method = "itau", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr)} else {res}
}

########################################################## SnC
gofRosenblattSnC = function (copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRosenblattSnC.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattSnC", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnC", estim.method = "mpl", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnC", estim.method = "itau", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr)} else {res}
}

########################################################## AnChisq
gofRosenblattChisq = function (copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofADChisq.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattChisq", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "AnChisq", estim.method = "mpl", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "AnChisq", estim.method = "itau", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr)} else {res}
}

########################################################## RosenblattGamma
gofRosenblattGamma = function (copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofADGamma.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattGamma", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "AnGamma", estim.method = "mpl", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "AnGamma", estim.method = "itau", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr)} else {res}
}

########################################################## Sn
gofSn = function (copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofSn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofSn", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Sn", estim.method = "mpl", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Sn", estim.method = "itau", processes = processes, param.est = param.est, df.est = df.est, dispstr=dispstr)} else {res}
}
