########################################################################
gofWhite = function(copula, x, M = 1000, param = 0.5, param.est = T, df = 4, df.est = T, margins = "ranks", execute.times.comp = T, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  
  if (execute.times.comp == T & M >= 100){
    times.comp = gofCheckTime(copula = copula, x=x, test = "gofWhite", M = M, print.res = F, processes = processes)
    print(.get.time(times.comp))
  }
  
  if (copula == "gaussian"){ warning("Please note that the old (pre 0.1-3) term 'gaussian' was replaced with 'normal'."); copula = "normal"}
  if (copula == "independence"){fam = 0
  } else if (copula == "normal"){fam = 1
  } else if (copula == "t"){fam = 2
  } else if (copula == "clayton"){fam = 3
  } else if (copula == "gumbel"){fam = 4
  } else if (copula == "frank"){fam = 5
  } else {stop("This copula is not implemented for gofWhite.")}
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
      param = try(fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "mpl")@estimate, silent = T); estim.method = "mpl"
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(ellipCopula(copula, dim = dim(x)[2], df = df, df.fixed = df.fixed), data = x, method = "itau")@estimate; estim.method = "itau"}
    }
    if (copula == "t" & df.fixed == F & param.est == T){
      df = tail(param, n=1)
      param = param[-length(param)]
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T); estim.method = "mpl"
    } else {
      copula = ellipCopula(copula, param = param, dim = dim(x)[2], df = df, df.fixed = T); estim.method = "mpl"
    }
  } else if("clayton" == copula || "frank" == copula || "gumbel" == copula || 
              "amh" == copula || "joe" == copula){
    if (param.est == T){
      param = try(fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "mpl")@estimate, silent = T); estim.method = "mpl"
      if (class(param) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameter failed. The estimation was performed with inversion of Kendall's Tau."); param = fitCopula(archmCopula(copula, dim = dim(x)[2]), data = x, method = "itau")@estimate; estim.method = "itau"}
    }
    copula = archmCopula(copula, param = param, dim = dim(x)[2]); estim.method = "mpl"
  }
    
  Tres = BiCopGofTest(u1 = x[,1], u2 = x[,2], family = fam, par = param, par2 = df, method = "white", B = 0)
  C.th.n = copula
  
  if (processes > 1) {
    cl = makeCluster(processes, type = "PSOCK")
    clusterEvalQ(cl, library(copula))
    clusterEvalQ(cl, library(foreach))
    clusterEvalQ(cl, library(gofCopula))
    clusterEvalQ(cl, "BiCopDeriv")
    clusterEvalQ(cl, "BiCopDeriv2")
    clusterEvalQ(cl, "BiCopPDF")
    clusterEvalQ(cl, "ginv")
    registerDoParallel(cl)
  } else {registerDoSEQ()}
  T0 = suppressWarnings(foreach(i=1:M) %dopar% {
    repeat {
      Uhat = rCopula(nrow(x), C.th.n)
      C.th.n. <- try(fitCopula(copula, Uhat, method = estim.method, 
                               estimate.variance = FALSE)@copula, silent = T)
      if (class(C.th.n.) != "try-error"){break}
    }
    T0. = BiCopGofTest(u1 = Uhat[,1], u2 = Uhat[,2], family = fam, par = C.th.n@parameters, par2 = if(class(C.th.n) == "tCopula"){C.th.n@df}else{0}, method = "white", B = 0)$statistic
    T0.
  })
  if (processes > 1) {
    stopCluster(cl)
  }
  
  structure(class = "gofCOP", 
            list(method = sprintf("Goodness-of-fit test based on White's information equality matrix"), 
                 erg.tests = matrix(c(mean(as.numeric(T0) >= as.numeric(Tres$statistic)), Tres$statistic, copula@parameters, if(class(try(copula@df, silent = TRUE)) == "try-error") {NULL} else {copula@df}), ncol = if(class(try(copula@df, silent = TRUE)) == "try-error") {3} else {4}, 
                                    dimnames = list("White", if(class(try(copula@df, silent = TRUE)) == "try-error") {c("p.value", "test statistic", "parameters")} else {c("p.value", "test statistic", "parameters", "df")}))))
}
