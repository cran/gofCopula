gofHybrid = function(copula, x, testset = c("gofPIOSRn", "gofKernel"), margins = "ranks", dispstr = "ex", M = 1000, execute.times.comp = T, param = 0.5, param.est = T, df = 4, df.est = T, m = 1, MJ = 100, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05, processes = 1){
  if(any(apply(as.matrix(testset), 2, function(x){is.element(x,c("gofRn"))}) == T)){stop("The test gofRn was removed due to inconsistencies with the remaining tests.")}
  if(any(!is.element(testset, gofWhich(copula, dim(x)[2])))) {stop("At least one of the tests in 'testset' is not implemented or can not handle a dataset of this dimension.")}
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
  
  res_list = mapply(function(k, Ms) tryCatch(doCall(.fcn = k, copula = copula, x = x, margins = margins, M = Ms, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, processes = processes)$erg.tests, error = function(e) warning(e)), testset, M, SIMPLIFY = FALSE)
  if (any(lapply(res_list, class) != "matrix")) {res_list = res_list[-which(lapply(res_list, class) != "matrix")]}
  res = do.call(rbind, res_list)
    
  
  if (length(testset) > 1) {
    which_comb =list()
    for(i in 1:(2^NROW(res))){
      which_comb[[i]]=which(as.integer(intToBits(i)) == 1)
    }
    comb_exist = which_comb[which(unlist(lapply(which_comb, length)) > 1)]
    
    pres = c()
    for (i in 1:length(comb_exist)){
      pres = c(pres, min(length(res[comb_exist[[i]],1]) * min(res[comb_exist[[i]],1]), 1))
    }
    hybrid_comb_names = paste("hybrid(", lapply(comb_exist, paste, collapse = ", "), ")", sep="")
    matrix_names = matrix(c(pres, rep(NA, length(pres)), rep(res[1,-c(1,2)], each = length(pres))), byrow = FALSE, nrow = length(pres))
    rownames(matrix_names) = hybrid_comb_names
    res1 = rbind(res, matrix_names)
  } else {
      res1 = res
  }
  structure(class = "gofCOP", 
            list(method = sprintf("Parametric bootstrap goodness-of-fit test with hybrid test and %s copula",
                                copula),
                 erg.tests = res1))
}
