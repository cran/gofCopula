gof = function(x, priority = "tests", copula = NULL, tests = NULL, margins = "ranks", dispstr = "ex", M = 50, MJ = 50, param = 0.5, param.est = T, df = 4, df.est = T, m = 1, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05, processes = 1){
  if (is.matrix(x) == F){stop("x must be a matrix")}
  if(any(lapply(tests, function(x){is.element(x,c("gofPIOSRn", "gofPIOSTn", "gofKernel", "gofRosenblattSnB", 
                                                  "gofRosenblattSnC", "gofADChisq", "gofADGamma", 
                                                  "gofSn", "gofKendallCvM", "gofKendallKS",
                                                  "gofWhite", "gofRn"))}) == F)==T){stop("At least one of the tests in 'testset' is not implemented")}
  if(any(lapply(copula, function(x){is.element(x,c("normal", "gaussian", "t", "clayton", "frank", "gumbel"))}) == F)==T){stop("At least one of the copulae is not implemented")}
  if (is.element(priority, c("tests", "copula")) == F){
    stop("Please insert a valid character string for the argument priority. It shall be either tests or copula.")
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
  
  if (!is.null(tests)){
    if (!is.null(copula)){
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, dispstr = dispstr, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, dispstr = dispstr, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
    } else {
      copula_list = lapply(tests, gofWhichCopula)
      copula = Reduce(intersect, copula_list)
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, dispstr = dispstr, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, dispstr = dispstr, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
    }
  } else {
    if (!is.null(copula)){
      tests_list = lapply(copula, gofWhich, d = dim(x)[2])
      tests = Reduce(intersect, tests_list)
      tests = tests[-which(tests == "gofHybrid")]
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, dispstr = dispstr, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, dispstr = dispstr, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
    } else if (priority == "tests") {
      tests = gofWhich("normal", dim(x)[2])
      tests = tests[-which(tests == "gofHybrid")]
      copula_list = lapply(tests, gofWhichCopula)
      copula = Reduce(intersect, copula_list)
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, dispstr = dispstr, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, dispstr = dispstr, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
    } else if (priority == "copula"){
      copula = c("normal", "t", "frank", "gumbel", "clayton")
      tests_list = lapply(copula, gofWhich, d = dim(x)[2])
      tests = Reduce(intersect, tests_list)
      tests = tests[-which(tests == "gofHybrid")]
      if (M >= 100 | MJ >= 100){
      times.comp = lapply(copula, gofCheckTime, x=x, test = tests, dispstr = dispstr, M = M, MJ = MJ, print.res = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
      print(.get.time(Reduce("+", times.comp)))
      }
      lapply(copula, gofHybrid, x=x, testset = tests, dispstr = dispstr, M = M, MJ = MJ, execute.times.comp = F, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)
    }
  }
}