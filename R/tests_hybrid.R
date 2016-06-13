gofHybrid = function(copula, x, testset = c("gofPIOSRn", "gofKernel"), margins = "ranks", dispstr = "ex", M = 1000, execute.times.comp = T, param = 0.5, param.est = T, df = 4, df.est = T, m = 1, MJ = 100, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05, processes = 1){
  if(any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofRn"))}) == T)){stop("The test gofRn was removed due to inconsistencies with the remaining tests.")}
  if(any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofPIOSRn", "gofPIOSTn", "gofKernel", "gofRosenblattSnB", 
                                                                 "gofRosenblattSnC", "gofRosenblattChisq", "gofRosenblattGamma", 
                                                                 "gofSn", "gofKendallCvM", "gofKendallKS",
                                                                 "gofWhite", "gofRn", "gofADChisq", "gofADGamma"))}) == F)==T){stop("At least one of the tests in 'testset' is not implemented")}
  if (dim(x)[2] == 3 & any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofKernel", "gofWhite"))}) == T)){
    stop("At least one test in the testset can't handle dimension 3. Please use the tests gofRosenblattSnB, 
          gofRosenblattSnC, gofRosenblattChisq, gofRosenblattGamma", "gofKendallCvM", "gofKendallKS", "gofPIOSRn", "gofPIOSTn", "gofSn or gofRn.")
  }
  if (dim(x)[2] > 3 & any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofPIOSRn", "gofPIOSTn", "gofKernel", "gofWhite"))}) == T)){
    stop("At least one test in the testset can't handle dimensions greater 3. Please use the tests gofRosenblattSnB, 
          gofRosenblattSnC, gofRosenblattChisq, gofRosenblattGamma", "gofKendallCvM", "gofKendallKS", "gofSn or gofRn.")
  }
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  
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
  
  if (execute.times.comp == T & M >= 100){
  times.comp = gofCheckTime(copula = copula, x=x, test = testset, M = M, dispstr = dispstr, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, m = m, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, processes = processes)
  print(.get.time(times.comp))
  }
  
  pres1 = c()
  pres = c()
  tres1 = c()
  tres = c()
  parares1 = c()
  parares = c()
  dfres1 = c()
  dfres = c()
  test_names = c()
  if (is.element("gofSn", testset) == T){
    cop = try(gofSn(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofSn")
  }
  if (is.element("gofPIOSRn", testset) == T){
    cop = try(gofPIOSRn(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofPIOSRn")
  }
  if (is.element("gofPIOSTn", testset) == T){
    cop = try(gofPIOSTn(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, m = m, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofPIOSTn")
  }
  if (is.element("gofKernel", testset) == T){
    cop = try(gofKernel(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofKernel")
  }
  if (is.element("gofKendallCvM", testset) == T){
    cop = try(gofKendallCvM(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofKendallCvM")
  }
  if (is.element("gofKendallKS", testset) == T){
    cop = try(gofKendallKS(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofKendallKS")
  }
  if (is.element("gofWhite", testset) == T){
    cop = try(gofWhite(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofWhite")
  }
  if (is.element("gofRosenblattSnB", testset) == T){
    cop = try(gofRosenblattSnB(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofRosenblattSnB")
  }
  if (is.element("gofRosenblattSnC", testset) == T){
    cop = try(gofRosenblattSnC(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofRosenblattSnC")
  }
  if (is.element("gofRosenblattChisq", testset) == T){
    cop = try(gofRosenblattChisq(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofRosenblattChisq")
  }
  if (is.element("gofADChisq", testset) == T){
    cop = try(gofADChisq(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofRosenblattChisq")
  }
  if (is.element("gofRosenblattGamma", testset) == T){
    cop = try(gofRosenblattGamma(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofRosenblattGamma")
  }
  if (is.element("gofADGamma", testset) == T){
    cop = try(gofADGamma(copula = copula, x = x, margins = margins, M = M, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes))
    if (class(cop) == "try-error"){pres = tres = parares = dfres = NA} else {
      pres = cop$erg.tests[1]
      tres = cop$erg.tests[2]
      parares = if (copula == "t") {cop$erg.tests[-c(1,2,length(cop$erg.tests))]} else {cop$erg.tests[-c(1,2)]}
      dfres = if (copula == "t") {cop$erg.tests[length(cop$erg.tests)]}
    }
    parares1 = c(parares1, parares)
    dfres1 = c(dfres1, dfres)
    pres1 = c(pres1, pres)
    tres1 = c(tres1, tres)
    test_names = c(test_names,"gofRosenblattGamma")
  }
  
  
  if (length(testset) > 1) {
    which_comb =list()
    for(i in 1:(2^length(testset))){
      which_comb[[i]]=which(as.integer(intToBits(i)) == 1)
    }
    comb_exist = which_comb[which(unlist(lapply(which_comb, length)) > 1)]
    
    for (i in 1:length(comb_exist)){
      pres1 = c(pres1, min(length(pres1[comb_exist[[i]]]) * min(pres1[comb_exist[[i]]]), 1))
      tres1 = c(tres1, NaN)
    }
    hybrid_comb_names = paste("hybrid(", lapply(comb_exist, paste, collapse = ", "), ")", sep="")
    matrix_names = c(substring(test_names, 4), hybrid_comb_names)
    structure(class = "gofCOP",
              list(method = sprintf("Hybrid Goodness-of-fit test with the testset = c(%s) and the %s copula.", 
                                    toString(testset), copula), 
                   erg.tests = matrix(if (copula == "t") {c(pres1, tres1, rep(parares, each=length(matrix_names)), rep(dfres, each=length(matrix_names)))} else {c(pres1, tres1, rep(parares, each=length(matrix_names)))}, nrow = length(pres1), 
                                      dimnames = list(matrix_names, if (copula == "t") {c("p.value", "test statistic", paste("rho.", 1:length(parares), sep=""), "df")} else {c("p.value", "test statistic", paste("rho.", 1:length(parares), sep=""))}))))
  } else {
  matrix_names = substring(testset, 4)
  structure(class = "gofCOP",
            list(method = sprintf("Goodness-of-fit test with the test = c(%s) and the %s copula.", 
                                  toString(testset), copula), 
                 erg.tests = matrix(if (copula == "t") {c(pres1, tres1, rep(parares, each=length(matrix_names)), rep(dfres, each=length(matrix_names)))} else {c(pres1, tres1, rep(parares, each=length(matrix_names)))}, nrow = length(pres1), 
                                    dimnames = list(matrix_names, if (copula == "t") {c("p.value", "test statistic", paste("rho.", 1:length(parares), sep=""), "df")} else {c("p.value", "test statistic", paste("rho.", 1:length(parares), sep=""))}))))
  }
}
