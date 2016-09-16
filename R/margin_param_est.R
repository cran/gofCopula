.margins.param.est = function(copula, margins, x, param, param.est, df, df.est, dispstr) {
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
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T); estim.method = "mpl"
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate}; estim.method = "itau"
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      copula = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr); estim.method = "mpl"
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr); estim.method = "mpl"
    }
    if (class(copula) == "indepCopula") {stop("The parameter estimation is at boundary and an independence copula was returned.")}
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
            "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T); estim.method = "mpl"
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate}; estim.method = "itau"
    }
    if (copula == "clayton" & dim(x)[2] > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
    if (copula == "frank" & dim(x)[2] > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
    copula = archmCopula(copula, param = param, dim = dim(x)[2]); estim.method = "mpl"
    if (class(copula) == "indepCopula") {stop("The parameter estimation is at boundary and an independence copula was returned.")}
  }
  return(list(copula, x, estim.method))
}