########################################################## SnB
gofRosenblattSnB = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRosenblattSnB.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattSnB", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
tryCatch(.gofRosenblattSnB(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins, dispstr = dispstr, processes = processes), error = .ErrorHandler)
}
.gofRosenblattSnB = function (copula, x, M, param, param.est, df, df.est, margins, dispstr, processes) {
  if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    res.margins = .margins(x, margins)
    param.margins = list()
    if (margins == "ranks"){
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[1]]}
    } else {
      for (i in 1:length(res.margins)) {param.margins[[i]] = res.margins[[i]][[1]]}
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[2]]}
    }
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      copula = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    } else {
    copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "SnB", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "SnB", estim.method = "itau", processes = processes)} else {res}
}

########################################################## SnC
gofRosenblattSnC = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofRosenblattSnC.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattSnC", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  tryCatch(.gofRosenblattSnC(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins, dispstr = dispstr, processes = processes), error = .ErrorHandler)
}
  .gofRosenblattSnC = function (copula, x, M, param, param.est, df, df.est, margins, dispstr, processes) {
    if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    res.margins = .margins(x, margins)
    param.margins = list()
    if (margins == "ranks"){
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[1]]}
    } else {
      for (i in 1:length(res.margins)) {param.margins[[i]] = res.margins[[i]][[1]]}
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[2]]}
    }
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      copula = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "SnC", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "SnC", estim.method = "itau", processes = processes)} else {res}
  
}

########################################################## AnChisq
gofRosenblattChisq = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofADChisq.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattChisq", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  tryCatch(.gofRosenblattChisq(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins, dispstr = dispstr, processes = processes), error = .ErrorHandler)
}
.gofRosenblattChisq = function (copula, x, M, param, param.est, df, df.est, margins, dispstr, processes) {
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    res.margins = .margins(x, margins)
    param.margins = list()
    if (margins == "ranks"){
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[1]]}
    } else {
      for (i in 1:length(res.margins)) {param.margins[[i]] = res.margins[[i]][[1]]}
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[2]]}
    }
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      copula = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "AnChisq", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "AnChisq", estim.method = "itau", processes = processes)} else {res}
}

########################################################## RosenblattGamma
gofRosenblattGamma = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofADGamma.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofRosenblattGamma", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  tryCatch(.gofRosenblattGamma(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins, dispstr = dispstr, processes = processes), error = .ErrorHandler)
}
.gofRosenblattGamma = function (copula, x, M, param, param.est, df, df.est, margins, dispstr, processes) {
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    res.margins = .margins(x, margins)
    param.margins = list()
    if (margins == "ranks"){
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[1]]}
    } else {
      for (i in 1:length(res.margins)) {param.margins[[i]] = res.margins[[i]][[1]]}
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[2]]}
    }
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      copula = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "AnGamma", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "AnGamma", estim.method = "itau", processes = processes)} else {res}
}

########################################################## Sn
gofSn = function (copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1) {
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofSn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofSn", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  tryCatch(.gofSn(copula = copula, x = x, M = M, param = param, param.est = param.est, df = df, df.est = df.est, margins = margins, dispstr = dispstr, processes = processes), error = .ErrorHandler)
}
.gofSn = function (copula, x, M, param, param.est, df, df.est, margins, dispstr, processes) {
  if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}
  
  if (any(x > 1) || any(x < 0)){
    if (margins == "ranks"){
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " of the observations.", sep = ""))
    } else {
      warning(paste("The observations aren't in [0,1]. The margins will be estimated by the ", margins, " distribution.", sep = ""))
    }
    
    res.margins = .margins(x, margins)
    param.margins = list()
    if (margins == "ranks"){
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[1]]}
    } else {
      for (i in 1:length(res.margins)) {param.margins[[i]] = res.margins[[i]][[1]]}
      for (i in 1:dim(x)[2]) {x[,i] = res.margins[[i]][[2]]}
    }
  }
  
  if ("normal" == copula || "t" == copula){
    if (param.est == T){
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate}
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      copula = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr)
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T)
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  res = try(.gofCopulapb(copula = copula, x = x, N = M, method = "Sn", estim.method = "mpl", processes = processes), silent = T)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau."); .gofCopulapb(copula = copula, x = x, N = M, method = "Sn", estim.method = "itau", processes = processes)} else {res}
}

