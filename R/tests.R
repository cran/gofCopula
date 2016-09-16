gofPIOSRn = function(copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] < 2 || dim(x)[2] > 3){stop("x must be of dimension 2 or 3")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofPIOSRn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofPIOSRn", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Rn", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Rn", estim.method = "itau", processes = processes)} else {res}
}

#################################################################################
gofPIOSTn = function(copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", m = 1, execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] < 2 || dim(x)[2] > 3){stop("x must be of dimension 2 or 3")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofPIOSTn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  n = dim(x)[1]
  B = n / m
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofPIOSTn", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  add.parameters = list(B, m, param.est, "mpl")
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Tn", estim.method = "mpl", processes = processes, add.parameters = add.parameters), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); add.parameters = list(B, m, param.est, "itau"); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Tn", estim.method = "itau", processes = processes, add.parameters = add.parameters)} else {res}
}

#################################################################################
gofKernel = function(copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", MJ = 100, delta.J = 0.5, nodes.Integration = 12, execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofKernel.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  add.parameters = list(nodes.Integration, MJ, delta.J)
  
  if (execute.times.comp == T & M >= 100 | execute.times.comp == T & MJ >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofKernel", dispstr = dispstr, M = M, MJ = MJ, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }

  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Kernel", estim.method = "mpl", processes = processes, add.parameters = add.parameters), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Kernel", estim.method = "itau", processes = processes, add.parameters = add.parameters)} else {res}
}
