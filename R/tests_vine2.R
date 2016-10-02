########################################################################
gofWhite = function(copula = c("normal", "t", "clayton", "frank", "gumbel"), x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}

  if (copula == "independence"){fam = 0
  } else if (copula == "normal" || copula == "gaussian"){fam = 1
  } else if (copula == "t"){fam = 2
  } else if (copula == "clayton"){fam = 3
  } else if (copula == "gumbel"){fam = 4
  } else if (copula == "frank"){fam = 5
  } else {stop("This copula is not implemented for gofWhite.")}
  add.parameters = list(fam)
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofWhite", M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr="ex")
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "White", estim.method = "mpl", processes = processes, add.parameters = add.parameters, param.est = param.est, df.est = df.est, dispstr="ex"), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "White", estim.method = "itau", processes = processes, add.parameters = add.parameters, param.est = param.est, df.est = df.est, dispstr="ex")} else {res}
}
