gofHybrid = function(copula, x, tests = c("gofPIOSRn", "gofKernel"), customTests = NULL, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, MJ = 100, dispstr = "ex", m = 1, delta.J = 0.5, nodes.Integration = 12, seed.active = NULL, processes = 1){
  if (is.matrix(x) == FALSE){stop("x must be a matrix")}
  if (!is.null(tests) & any(is.element(tests, c("gofRn")))){stop("The test gofRn was removed due to inconsistencies with the remaining tests.")}
  if (any(!is.element(tests, gofTest4Copula(copula, dim(x)[2]))) || !is.null(customTests) & any(!is.element(customTests, ls(".GlobalEnv")))) {stop("At least one of the tests in 'tests' is not implemented, cannot handle a dataset of this dimension or at least one of the tests in 'customTests' does not match any function in the global workspace. Please check if it is correctly spelled in the function call.")}
  if (any(!vapply(customTests, function(x) all(names(formals(x)) %in% c("x", "copula")), TRUE))) {stop("At least one function in 'customTest' does not follow the requirements for the arguments. The first argument for the dataset has to be called 'x', the second one for the copula has to be called 'copula'.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  if (any(!is.element(margins, c("ranks", "beta", "cauchy", "chisq", "f", "gamma", "lnorm", "norm", "t", "weibull","exp")))) {stop(paste("At least one of the distributions in `margins' is not implemented. Please amend and run the function again. \n It has to be either of `ranks', `beta', `cauchy', `chisq', `f', `gamma', `lnorm', `norm', `t', `weibull', `exp'."))}
  if (!is.numeric(processes)) {stop("The argument 'processes' has to be a numeric.")}
  if (processes %% 1 != 0 | processes < 1) {stop("The argument 'processes' has to be a positiv integer.")}
  if (!is.numeric(M)) {stop("The argument 'M' has to be a numeric.")}
  if (M %% 1 != 0 | M < 0) {stop("The argument 'M' has to be a positiv integer.")}
  if (!is.numeric(param)) {stop("The argument 'param' has to be a numeric.")}
  if (!is.numeric(df)) {stop("The argument 'df' has to be a numeric.")}
  if (!is.numeric(delta.J)) {stop("The argument 'delta.J' has to be a numeric.")}
  if (delta.J <= 0) {stop("The argument 'delta.J' has to be larger 0.")}
  if (!is.numeric(nodes.Integration)) {stop("The argument 'nodes.Integration' has to be a numeric.")}
  if (nodes.Integration %% 1 != 0 | nodes.Integration < 0) {stop("The argument 'nodes.Integration' has to be a positiv integer.")}
  if (!is.numeric(MJ)) {stop("The argument 'MJ' has to be a numeric.")}
  if (MJ %% 1 != 0 | MJ < 0) {stop("The argument 'MJ' has to be a positiv integer.")}
  if (!is.numeric(m)) {stop("The argument 'm' has to be a numeric.")}
  n = dim(x)[1]
  if (n %% m != 0 | m < 1) {stop("The length of the blocks, 'm', has to be larger 1 and a divisor of the length of the data sequence.")}
  if (!inherits(param.est, "logical")) {stop("The argument 'param.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!inherits(df.est, "logical")) {stop("The argument 'df.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!is.null(seed.active) & length(seed.active) != 1 & length(seed.active) != (M+1)) {stop("The seed has to be an integer or a vector of M+1 seeds.")}
  if (!is.null(seed.active) & length(seed.active) == 1) {set.seed(seed.active); RNGkind(sample.kind = "default"); seed.active = sample(x = 2147483647, size = M+1)}
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")}
  
  param.margins = NULL
  if (!is.null(margins)){
    print(paste("The margins will be estimated as: ", paste0(margins, collapse = ", "), sep = ""))
    
    res.margins = .margins(x, margins)
    param.margins = list()
    if(length(margins) == 1) {margins.dummy = rep(margins, dim(x)[2])} else {margins.dummy = margins}
    for (i in 1:length(margins.dummy)) {
      if (margins.dummy[i] == "ranks"){
        x[,i] = res.margins[[i]][[1]]
      } else {
        param.margins[[i]] = res.margins[[i]][[1]]
        x[,i] = res.margins[[i]][[2]]
      }
    }
  }
  
  
  res_list = list()
  if (!is.null(tests)) {
    res_list = mapply(function(k, Ms) {print(paste0("Test ", k, " is running"));a = tryCatch(doCall(.fcn = k, copula = copula, x = x, margins = NULL, M = Ms, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, processes = processes, seed.active = seed.active), error = function(e) warning(e)); a}, tests, M, SIMPLIFY = FALSE)
    if (any(unlist(lapply(res_list, function(x) !any(is.element(class(x), "gofCOP")))))) {res_list = res_list[-which(lapply(res_list, function(x) !any(is.element(class(x), "gofCOP"))) == TRUE)]}
  }
  
  res_list2 = list()
  if (!is.null(customTests)) {
    res_list2 = mapply(function(k, Ms) {print(paste0("Test ", k, " is running"));a = tryCatch(gofCustomTest(copula = copula, x = x, customTest = k, margins = NULL, M = Ms, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, processes = processes, seed.active = seed.active), error = function(e) warning(e)); a}, customTests, M, SIMPLIFY = FALSE)
    if (any(unlist(lapply(res_list2, function(x) !any(is.element(class(x), "gofCOP")))))) {res_list2 = res_list2[-which(lapply(res_list2, function(x) !any(is.element(class(x), "gofCOP"))) == TRUE)]}
  } 
  res = do.call(rbind, c(lapply(res_list, function(x) x[[1]]$res.tests), lapply(res_list2, function(x) x[[1]]$res.tests)))
  resTheta = do.call(rbind, c(lapply(res_list, function(x) x[[1]]$theta), lapply(res_list2, function(x) x[[1]]$theta)))
  resDf = do.call(rbind, c(lapply(res_list, function(x) x[[1]]$df), lapply(res_list2, function(x) x[[1]]$df)))
    
  
  if (length(tests) > 1) {
    which_comb =list()
    for(i in 1:(2^NROW(res))){
      which_comb[[i]] = which(as.integer(intToBits(i)) == 1)
    }
    comb_exist = which_comb[which(unlist(lapply(which_comb, length)) > 1)]
    
    pres = c()
    for (i in 1:length(comb_exist)){
      pres = c(pres, min(length(res[comb_exist[[i]], 1]) * min(res[comb_exist[[i]], 1]), 1))
    }
    hybrid_comb_names = paste("hybrid(", lapply(comb_exist, paste, collapse = ", "), ")", sep="")
    matrix_names = matrix(c(pres, rep(NA, length(pres)), rep(res[1,-c(1,2)], each = length(pres))), byrow = FALSE, nrow = length(pres))
    rownames(matrix_names) = hybrid_comb_names
    res1 = rbind(res, matrix_names)
  } else {
      res1 = res
  }
  res = structure(class = "gofCOP", 
                  list(
            list(method = sprintf("Parametric bootstrap goodness-of-fit test with hybrid test and %s copula",
                                copula),
                 copula = copula, 
                 margins = margins, 
                 param.margins = param.margins, 
                 theta = unique(unname(resTheta)),
                 df = unique(unname(resDf)),
                 res.tests = res1)))
  names(res) = copula
  res
}
