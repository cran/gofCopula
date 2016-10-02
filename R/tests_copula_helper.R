.gofCopulapb = function (copula, x, M = 1000, method = eval(formals(.Tstats)$Tstat), estim.method = eval(formals(fitCopula)$method), processes, param.est, df.est, dispstr, ...) {

  bs.ac.c = c()

  if (is.element(method, c("SnB", "SnC", "AnChisq", "AnGamma"))){
      x = do.call(cCopula, c(list(x, copula = copula)))
  }
  ac.c <- if (method == "Sn") {
    .Tstats(x, Tstat = method, copula = copula)
  } else if (method == "Tn" || method == "Kernel" || method == "White") {
    add.parameters = list(...)$add.parameters
    .Tstats(x, Tstat = method, copula = copula, add.parameters = add.parameters)
  } else {
    .Tstats(x, Tstat = method, copula = copula)
  } 
  
  if (processes > 1) {
  cl = makeCluster(processes, type = "PSOCK")
  clusterEvalQ(cl, library(copula))
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(gofCopula))
  if (method == "White") {
      clusterEvalQ(cl, "BiCopDeriv")
      clusterEvalQ(cl, "BiCopDeriv2")
      clusterEvalQ(cl, "BiCopPDF")
      clusterEvalQ(cl, "ginv")
  }
  registerDoParallel(cl)
  } else {registerDoSEQ()}
  bs.ac.c = foreach(i=1:M) %dopar% {
    repeat {
      xsim = rCopula(nrow(x), copula)
      copula.sim = try(.margins.param.est(copula=as.character(substr(class(copula),1,nchar(class(copula))-6)), margins="ranks", x=xsim, 
                        param=if(class(copula) != "tCopula"){copula@parameters}else{copula@parameters[-length(copula@parameters)]}, 
                        param.est=param.est, 
                        df=if(class(copula) != "tCopula"){4}else{copula@parameters[length(copula@parameters)]}, 
                        df.est=df.est, 
                        dispstr=dispstr)[[1]], silent = T)
      if (class(copula.sim) != "try-error"){break}
    }
    if (is.element(method, c("SnB", "SnC", "AnChisq", "AnGamma"))){
        xsim = do.call(cCopula, c(list(xsim, copula = copula.sim)))
    }
    if (method == "Sn") {
      .Tstats(xsim, Tstat = method, copula = copula.sim)
    } else if (method == "Tn" || method == "Kernel" || method == "White") {
      .Tstats(xsim, Tstat = method, copula = copula.sim, add.parameters = add.parameters)
    } else {
      .Tstats(xsim, Tstat = method, copula = copula.sim)
    }
  }
  if (processes > 1) {stopCluster(cl)}
  
  ac.c = as.numeric(ac.c)
  bs.ac.c = as.numeric(bs.ac.c)
  test = sum(abs(bs.ac.c) >= abs(ac.c))/M
  
  switch(method,
         SnB = {matrix_names = "RosenblattSnB"},
         SnC = {matrix_names = "RosenblattSnC"},
         AnChisq = {matrix_names = "RosenblattChisq"},
         AnGamma = {matrix_names = "RosenblattGamma"},
         Sn = {matrix_names = "Sn"},
         Rn = {matrix_names = "PIOSRn"},
         Tn = {matrix_names = "PIOSTn"},
         Kernel = {matrix_names = "Kernel"},
         SnK = {matrix_names = "KendallCvM"},
         TnK = {matrix_names = "KendallKS"},
         White = {matrix_names = "White"})
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test with %s test and %s copula", 
                                                   matrix_names, as.character(substr(class(copula),1,nchar(class(copula))-6))),
                 erg.tests = matrix(c(test, ac.c, copula@parameters), ncol = 2 + length(copula@parameters),  
                                    dimnames = list(matrix_names, if(class(copula) != "tCopula") {c("p.value", 
                                        "test statistic", paste("rho.", 1:length(copula@parameters), sep=""))} else {
                                            c("p.value", "test statistic", paste("rho.", 1:(length(copula@parameters)-1), 
                                                                                 sep=""), "df")}))))
}

