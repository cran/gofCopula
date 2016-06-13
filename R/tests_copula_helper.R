.gofCopulapb = function (copula, x, N = 1000, method = eval(formals(gofTstat)$method), estim.method = eval(formals(fitCopula)$method), processes, ...) 
{

  C.th.n <- fitCopula(copula, x, method = estim.method, 
                      estimate.variance = FALSE, ...)@copula

  Tstat <- if (method == "Sn") 
    gofTstat(x, method = method, copula = C.th.n)
  else gofTstat(x, method = method)
  
  if (processes > 1) {
  cl = makeCluster(processes, type = "PSOCK")
  clusterEvalQ(cl, library(copula))
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(gofCopula))
  registerDoParallel(cl)
  } else {registerDoSEQ()}
  T0 = foreach(i=1:N) %dopar% {
    repeat {
        Uhat = rCopula(nrow(x), C.th.n)
      C.th.n. <- try(fitCopula(copula, Uhat, method = estim.method, 
                               estimate.variance = FALSE, ...)@copula, silent = T)
      if (class(C.th.n.) != "try-error"){break}
    }
    u. = Uhat
    T0. <- if (method == "Sn") 
      gofTstat(u., method = method, copula = C.th.n.)
    else gofTstat(u., method = method)
    T0.
  }
  if (processes > 1) {stopCluster(cl)}
  
  
  switch(method,
         SnB = {matrix_names = "RosenblattSnB"},
         SnC = {matrix_names = "RosenblattSnC"},
         AnChisq = {matrix_names = "RosenblattChisq"},
         AnGamma = {matrix_names = "RosenblattGamma"},
         Sn = {matrix_names = "Sn"})
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test with %s test", 
                                                   matrix_names),
                 erg.tests = matrix(c(sum(abs(as.numeric(T0)) >= abs(Tstat))/N, Tstat, copula@parameters, if(class(try(copula@df, silent = TRUE)) == "try-error") {NULL} else {copula@df}), ncol = if(class(try(copula@df, silent = TRUE)) == "try-error") {2 + length(copula@parameters)} else {3 + length(copula@parameters)},  
                                    dimnames = list(matrix_names, if(class(try(copula@df, silent = TRUE)) == "try-error") {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""))} else {c("p.value", "test statistic", paste("rho.", 1:length(copula@parameters), sep=""), "df")}))))
  
}

