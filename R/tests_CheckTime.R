gofCheckTime = function(copula, x, tests = NULL, customTests = NULL, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, MJ = 100, dispstr = "ex", print.res = TRUE, m = 1, delta.J = 0.5, nodes.Integration = 12, seed.active = NULL, processes = 1){
  if (is.matrix(x) == FALSE){stop("x must be a matrix")}
  if (is.null(tests) & is.null(customTests)) {stop("Please provide either 'tests' or 'customTests'.")}
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
  if (!is.null(seed.active) & length(seed.active) == 1) {set.seed(seed.active); RNGkind(sample.kind = "default"); seed.active = sample(x=2147483647, size = M+1)}
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")}
  
  print("An estimate of the computational time is under derivation.")
  lasted.time = c()
  N = c(2,5,10,15)
  NJ = c(2,5,10,15)
  if (!is.null(tests)) {
  for (j in 1:length(tests)){
    times.comp = c()
    if (tests[j] == "gofKernel"){
      for (i in N){
        for (ii in NJ){
          times.comp = rbind(times.comp, system.time(invisible(capture.output(suppressWarnings(gofHybrid(copula = copula, x = x, dispstr = dispstr, tests = tests[j], M = i, MJ = ii, margins = margins, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, seed.active = seed.active, processes = processes)))))[3])
        }
      }
      times.comp = cbind(sort(rep(N, length(N))), rep(NJ, length(NJ)), times.comp)
      if(processes > 1){
        times.lm = lm(times.comp[, 3] ~ times.comp[,1] + times.comp[,2])
        lasted.time[j] = round(times.lm$coefficients[1] + times.lm$coefficients[2] * M + times.lm$coefficients[3] * MJ)
      }else{
        times.lm = lm(times.comp[, 3] ~ times.comp[,1] + times.comp[,2] - 1)
        lasted.time[j] = round(times.lm$coefficients[1] * M + times.lm$coefficients[2] * MJ)
      }
      
      
      if (lasted.time[j] < 0) {print(paste0("Derivation time could not be computed for ", tests[j])); lasted.time[j] = NA}
    } else {
      for (i in N){
        times.comp = rbind(times.comp, system.time(invisible(capture.output(suppressWarnings(gofHybrid(copula = copula, x = x, dispstr = dispstr, tests = tests[j], M = i, margins = margins, param = param, param.est = param.est, df = df, df.est = df.est, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, seed.active = seed.active, processes = processes)))))[3])
      }
      times.comp = cbind(N, times.comp)
      if(processes > 1){
        times.lm = lm(times.comp[, 2] ~ times.comp[, 1])
        lasted.time[j] = round(times.lm$coefficients[1] + times.lm$coefficients[2] * M)
      }else{
        times.lm = lm(times.comp[, 2] ~ times.comp[, 1] - 1)
        lasted.time[j] = round(times.lm$coefficients[1] * M)
      }

      if (lasted.time[j] < 0) {print(paste0("Derivation time could not be computed for ", tests[j])); lasted.time[j] = NA}
    }
  }
  }
  
  lasted.time2 = c()
  if (!is.null(customTests)) {
    for (j in 1:length(customTests)){
      times.comp = c()
      for (i in N){
        times.comp = rbind(times.comp, system.time(invisible(capture.output(suppressWarnings(gofHybrid(copula = copula, x = x, dispstr = dispstr, tests = NULL, customTests = customTests[j], M = i, margins = margins, param = param, param.est = param.est, df = df, df.est = df.est, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, seed.active = seed.active, processes = processes)))))[3])
      }
      times.comp = cbind(N, times.comp)
      if(processes > 1){
        times.lm = lm(times.comp[, 2] ~ times.comp[, 1])
        lasted.time2[j] = round(times.lm$coefficients[1] + times.lm$coefficients[2] * M)
      }else{
        times.lm = lm(times.comp[, 2] ~ times.comp[, 1] - 1)
        lasted.time2[j] = round(times.lm$coefficients[1] * M)
      }
      if (lasted.time2[j] < 0) {print(paste0("Derivation time could not be computed for ", customTests[j])); lasted.time2[j] = NA}
    }
  }
  if (print.res == TRUE){
    .get.time(sum(lasted.time, lasted.time2))
  } else {
    store.time = sum(lasted.time, lasted.time2)
  }
}
