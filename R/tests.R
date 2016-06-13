gofPIOSRn = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] < 2 || dim(x)[2] > 3){stop("x must be of dimension 2 or 3")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofPIOSRn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  dims = dim(x)[2]
  n = dim(x)[1]
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofPIOSRn", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  
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
    if (copula == "clayton" & dims > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
    if (copula == "frank" & dims > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  
  bs.ac.c = c()
  ac.c  = .Rn(x = x, copula = copula, dims = dims)
  if (processes > 1) {
  cl = makeCluster(processes, type = "PSOCK")
  clusterEvalQ(cl, library(gofCopula))
  clusterEvalQ(cl, library(copula))
  clusterEvalQ(cl, library(foreach))
  registerDoParallel(cl)
  } else {registerDoSEQ()}
  bs.ac.c = foreach(i=1:M) %dopar% {
      xsim = rCopula(n, copula)
    .Rn(xsim, copula = copula, dims = dims)
  }
  if (processes > 1) {
  stopCluster(cl)
  }
 bs.ac.c = as.numeric(bs.ac.c)
  test = sum(abs(bs.ac.c) >= abs(ac.c))/M
 
structure(class = "gofCOP", 
          list(method = sprintf("Parametric bootstrap goodness-of-fit test (PIOS)"),
               erg.tests = matrix(c(test, ac.c, copula@parameters, if(class(try(copula@df, silent = TRUE)) == "try-error") {NULL} else {copula@df}), ncol = if(class(try(copula@df, silent = TRUE)) == "try-error") {2 + length(copula@parameters)} else {3 + length(copula@parameters)},  
                                  dimnames = list("PIOSRn", if(class(try(copula@df, silent = TRUE)) == "try-error") {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""))} else {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""), "df")}))))
}

#################################################################################
gofPIOSTn = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", m = 1, execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] < 2 || dim(x)[2] > 3){stop("x must be of dimension 2 or 3")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofPIOSTn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  dims = dim(x)[2]
  n = dim(x)[1]
  B = n / m
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofPIOSTn", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  
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
    if (copula == "clayton" & dims > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
    if (copula == "frank" & dims > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
    copula = archmCopula(copula, param = param, dim = dim(x)[2])
  }
  
  bs.ac.c = c()
  ac.c  = .Tn(x = x, copula = copula, B = B, m = m, dims = dims, param.est = param.est)
  if (processes > 1) {
  server = as.character(Sys.info()["nodename"])
  cl = makeCluster(rep("localhost",processes), type = "PSOCK")
  clusterEvalQ(cl, library(gofCopula))
  clusterEvalQ(cl, library(copula))
  clusterEvalQ(cl, library(foreach))
  registerDoParallel(cl)
  } else {registerDoSEQ()}
  bs.ac.c = foreach(i=1:M) %dopar% {
      xsim = rCopula(n, copula)
    .Tn(xsim, copula = copula, B = B, m = m, dims = dims, param.est = param.est)
  }
  if (processes > 1) {
  stopCluster(cl)
  }
  bs.ac.c = as.numeric(bs.ac.c)
  test = sum(abs(bs.ac.c) >= abs(ac.c))/M
  
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test (approximate PIOS)"),
                 erg.tests = matrix(c(test, ac.c, copula@parameters, if(class(try(copula@df, silent = TRUE)) == "try-error") {NULL} else {copula@df}), ncol = if(class(try(copula@df, silent = TRUE)) == "try-error") {2 + length(copula@parameters)} else {3 + length(copula@parameters)},  
                                    dimnames = list("PIOSTn", if(class(try(copula@df, silent = TRUE)) == "try-error") {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""))} else {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""), "df")}))))
}

#################################################################################
gofKernel = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", MJ = 100, delta.J = 0.5, nodes.Integration = 12, execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofKernel.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  dims = 2
  n = dim(x)[1]
  
  if (execute.times.comp == T & M >= 100 | execute.times.comp == T & MJ >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofKernel", dispstr = dispstr, M = M, MJ = MJ, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }

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
  
  bs.ac.c = c()
  ac.c  = .Kernel(x = x, copula = copula, dims = dims, n = n, nodes.Integration = nodes.Integration, MJ = MJ, delta.J = delta.J)
  if (processes > 1) {
  server = as.character(Sys.info()["nodename"])
  cl = makeCluster(rep("localhost",processes), type = "PSOCK")
  clusterEvalQ(cl, library(gofCopula))
  clusterEvalQ(cl, library(copula))
  clusterEvalQ(cl, library(foreach))
  registerDoParallel(cl)
  } else {registerDoSEQ()}
  bs.ac.c = foreach(i=1:M) %dopar% {
      xsim = rCopula(n, copula)
    .Kernel(xsim, copula = copula, dims = dims, n = n, nodes.Integration = nodes.Integration, MJ = MJ, delta.J = delta.J)
  }
  if (processes > 1) {
  stopCluster(cl)
  }
  bs.ac.c = as.numeric(bs.ac.c)
  test = sum(abs(bs.ac.c) >= abs(ac.c))/M

  structure(class = "gofCOP", 
            list(method = sprintf("Semiparametric bootstrap goodness-of-fit test with Scaillet test"),
                 erg.tests = matrix(c(test, ac.c, copula@parameters, if(class(try(copula@df, silent = TRUE)) == "try-error") {NULL} else {copula@df}), ncol = if(class(try(copula@df, silent = TRUE)) == "try-error") {2 + length(copula@parameters)} else {3 + length(copula@parameters)},  
                                    dimnames = list("Kernel", if(class(try(copula@df, silent = TRUE)) == "try-error") {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""))} else {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""), "df")}))))
}
