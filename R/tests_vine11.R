gofKendallCvM = function(copula, x, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", dispstr = "ex", M = 100, execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  dims = dim(x)[2]
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == F){stop("This copula is not implemented for gofKendallCvM.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofKendallCvM", dispstr = dispstr, M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
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
  
  if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  if (df.est == T){df.fixed = F} else if (df.est == F){df.fixed = T}

#### test
if ("normal" == copula || "t" == copula){
  if (param.est == T){
    param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "mpl")@estimate, silent = T); estim.method = "mpl"
    if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed, dispstr = dispstr), data = x, method = "itau")@estimate; estim.method = "itau"}
  }
  if (copula == "t" & df.fixed == F & param.est == T){
    df = tail(param, n=1)
    cop = ellipCopula(copula, param = param[-length(param)], dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr); estim.method = "mpl"
  } else {
    cop = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T, dispstr = dispstr); estim.method = "mpl"
  }
} else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
            "amh" == copula || "joe" == copula){
  if (param.est == T){
    param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T); estim.method = "mpl"
    if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate; estim.method = "itau"}
  }
  if (copula == "clayton" & dims > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the clayton copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
  if (copula == "frank" & dims > 2 & param < 0) {stop("The dependence parameter is negative for the dataset. For the frank copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")}
  cop = archmCopula(copula, param = param, dim = dim(x)[2]); estim.method = "mpl"
}

n = dim(x)[1]
Cn = C.n(x,x)
Kn = C.n(as.matrix(Cn), as.matrix(Cn))
if (is.element(copula, c("clayton", "frank", "gumbel"))){
  tcop = seq(1, n)/n
  texe = sapply(Kn[-n], function(x,y) length(which(x <= y)), tcop)
 Scvm1 = sum(texe[-n]^2 * (apply(matrix(tcop[-1], nrow = (n-1), ncol = dims),1,pCopula,cop) - apply(matrix(tcop[-n], nrow = (n-1), ncol = dims),1,pCopula,cop)))/n
 Scvm2 = sum(texe[-n] * (apply(matrix(tcop[-1], nrow = (n-1), ncol = dims),1,pCopula,cop)^2 - apply(matrix(tcop[-n], nrow = (n-1), ncol = dims),1,pCopula,cop)^2))
  SnK = n/3 + Scvm1 + Scvm2 
} else if (is.element(copula, c("normal", "t"))) {
  Csample = rCopula(n*2, cop)
  Csample.n = C.n(Csample,Csample)
  Ksample = C.n(as.matrix(Csample.n), as.matrix(Csample.n))
  SnK = n * mean((Kn - Ksample)^2)
}

Scvm = c()
if (processes > 1) {
cl = makeCluster(processes, type = "PSOCK")
clusterEvalQ(cl, library(copula))
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)
} else {registerDoSEQ()}
Scvm = foreach(i=1:M) %dopar% {
  repeat {
    Csampleb = rCopula(n, cop)
    copb <- try(fitCopula(cop, Csampleb, method = estim.method, 
                          estimate.variance = FALSE)@copula, silent = T)
    if (class(copb) != "try-error"){break}
  }
  Csample.nb = C.n(Csampleb,Csampleb)
  Ksampleb = C.n(as.matrix(Csample.nb), as.matrix(Csample.nb))
  
  if (is.element(copula, c("clayton", "frank", "gumbel"))){
    tcop = seq(1, n)/n
    texe = sapply(Ksampleb[-n], function(x,y) length(which(x <= y)), tcop)
   Scvm1 = sum(texe[-n]^2 * (apply(matrix(tcop[-1], nrow = (n-1), ncol = dims),1,pCopula,copb) - apply(matrix(tcop[-n], nrow = (n-1), ncol = dims),1,pCopula,copb)))/n
   Scvm2 = sum(texe[-n] * (apply(matrix(tcop[-1], nrow = (n-1), ncol = dims),1,pCopula,copb)^2 - apply(matrix(tcop[-n], nrow = (n-1), ncol = dims),1,pCopula,copb)^2))
    n/3 + Scvm1 + Scvm2 
  } else if (is.element(copula, c("normal", "t"))) {
    Csample = rCopula(n*2, cop)
    Csample.n = C.n(Csample,Csample)
    Ksample = C.n(as.matrix(Csample.n), as.matrix(Csample.n))
    n * mean((Ksampleb - Ksample)^2)
  }
}
if (processes > 1) {
stopCluster(cl)
}
Scvm = as.numeric(Scvm)
pvalue = sum(abs(Scvm) >= abs(SnK))/M

#### end test

    structure(class = "gofCOP", 
              list(method = sprintf("Parametric bootstrap goodness-of-fit test (Cramer-von Mises) based on Kendall's process"), 
                   erg.tests = matrix(c(pvalue, SnK, cop@parameters, if(class(try(cop@df, silent = TRUE)) == "try-error") {NULL} else {cop@df}), ncol = if(class(try(cop@df, silent = TRUE)) == "try-error") {2 + length(cop@parameters)} else {3 + length(cop@parameters)}, 
                                      dimnames = list("KendallCvM", if(class(try(cop@df, silent = TRUE)) == "try-error") {c("p.value", "test statistic", paste("rho.", 1:length(cop@parameters), sep=""))} else {c("p.value", "test statistic", paste("rho.", 1:length(cop@parameters), sep=""), "df")}))))
}
