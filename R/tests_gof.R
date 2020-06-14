gof <- function(x, priority = "tests", copula = NULL, tests = NULL, customTests = NULL, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, MJ = 100, dispstr = "ex", m = 1, delta.J = 0.5, nodes.Integration = 12, lower = NULL, upper = NULL, seed.active = NULL, processes = 1) {
  if (is.matrix(x) == FALSE) {
    stop("x must be a matrix")
  }
  if (is.null(copula)) {
    tests.available <- gofTest4Copula(copula, dim(x)[2])
  } else {
    tests.available <- Reduce(intersect, lapply(copula, function(y) gofTest4Copula(y, dim(x)[2])))
  }
  if (!is.null(tests) & any(!is.element(tests, tests.available)) || !is.null(customTests) & any(!is.element(customTests, ls(".GlobalEnv")))) {
    stop("At least one of the tests in 'tests' is not implemented, cannot handle a dataset of this dimension or at least one of the tests in 'customTests' does not match any function in the global workspace. Please check if it is correctly spelled in the function call.")
  }
  if (any(lapply(copula, function(x) {
    is.element(x, c("normal", "gaussian", "t", "clayton", "frank", "gumbel"))
  }) == FALSE) == TRUE) {
    stop("At least one of the copulae is not implemented")
  }
  if (is.element(priority, c("tests", "copula")) == FALSE) {
    stop("Please insert a valid character string for the argument priority. It shall be either 'tests' or 'copula'.")
  }
  if (!is.element(dispstr, c("ex", "un"))) {
    stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")
  }
  if (!is.numeric(processes)) {
    stop("The argument 'processes' has to be a numeric.")
  }
  if (processes %% 1 != 0 | processes < 1) {
    stop("The argument 'processes' has to be a positive integer.")
  }
  if (!is.numeric(M)) {
    stop("The argument 'M' has to be a numeric.")
  }
  if (M %% 1 != 0 | M < 0) {
    stop("The argument 'M' has to be a positive integer.")
  }
  if (!is.numeric(param)) {
    stop("The argument 'param' has to be a numeric.")
  }
  if (!is.numeric(df)) {
    stop("The argument 'df' has to be a numeric.")
  }
  if (!is.numeric(delta.J)) {
    stop("The argument 'delta.J' has to be a numeric.")
  }
  if (delta.J <= 0) {
    stop("The argument 'delta.J' has to be larger 0.")
  }
  if (!is.numeric(nodes.Integration)) {
    stop("The argument 'nodes.Integration' has to be a numeric.")
  }
  if (nodes.Integration %% 1 != 0 | nodes.Integration < 0) {
    stop("The argument 'nodes.Integration' has to be a positive integer.")
  }
  if (!is.numeric(MJ)) {
    stop("The argument 'MJ' has to be a numeric.")
  }
  if (MJ %% 1 != 0 | MJ < 0) {
    stop("The argument 'MJ' has to be a positive integer.")
  }
  if (!is.numeric(m)) {
    stop("The argument 'm' has to be a numeric.")
  }
  n <- dim(x)[1]
  if (n %% m != 0 | m < 1) {
    stop("The length of the blocks, 'm', has to be larger 1 and a divisor of the length of the data sequence.")
  }
  if (!inherits(param.est, "logical")) {
    stop("The argument 'param.est' has to be either 'TRUE' or 'FALSE'.")
  }
  if ((is.null(copula) || is.element("gumbel", copula)) & param.est == FALSE & param <= 1) {
    warning("When copula is 'gumbel', 'param' has to be larger 1. Because 'param.est' was set to 'FALSE', 'param' will be set to 1.5 as default value for 'gumbel' copula in the individual tests.")
  }
  if (!inherits(df.est, "logical")) {
    stop("The argument 'df.est' has to be either 'TRUE' or 'FALSE'.")
  }
  if (!is.null(seed.active) & length(seed.active) != 1 & length(seed.active) != (M + 1)) {
    stop("The seed has to be an integer or a vector of M+1 seeds.")
  }
  if (!is.null(seed.active) & length(seed.active) == 1) {
    set.seed(seed.active)
    RNGsetting <- RNGkind()
    RNGkind(sample.kind = "default")
    on.exit(RNGkind(sample.kind = RNGsetting[3]))
    seed.active <- sample(x = 2147483647, size = M + 1)
  }
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {
    stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")
  }
  if (any(!vapply(customTests, function(x) all(names(formals(x)) %in% c("x", "copula")), TRUE))) {
    stop("At least one function in 'customTest' does not follow the requirements for the arguments. The first argument for the dataset has to be called 'x', the second one for the copula has to be called 'copula'.")
  }

  if (!is.null(tests)) {
    if (!is.null(copula)) {
      res <- sapply(copula, function(cop) {
        cat(paste(cop, "copula\n"))
        tmp <- suppressWarnings(.gofHybrid(x = x, copula = cop, tests = tests, customTests = customTests, margins = margins, dispstr = dispstr, M = M, MJ = MJ, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, lower = lower, upper = upper, seed.active = seed.active, processes = processes))
        cat("\n")
        return(tmp)
      }, USE.NAMES = FALSE)
      return(structure(
        class = "gofCOP",
        res
      ))
    } else {
      copula_list <- lapply(tests, gofCopula4Test)
      copula <- Reduce(intersect, copula_list)
      res <- sapply(copula, function(cop) {
        cat(paste(cop, "copula\n"))
        tmp <- suppressWarnings(.gofHybrid(x = x, copula = cop, tests = tests, customTests = customTests, margins = margins, dispstr = dispstr, M = M, MJ = MJ, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, lower = lower, upper = upper, seed.active = seed.active, processes = processes))
        cat("\n")
        return(tmp)
      }, USE.NAMES = FALSE)
      return(structure(
        class = "gofCOP",
        res
      ))
    }
  } else {
    if (!is.null(copula)) {
      tests_list <- lapply(copula, gofTest4Copula, d = dim(x)[2])
      tests <- Reduce(intersect, tests_list)
      tests <- tests[-c(which(tests == "gofHybrid"), which(tests == "gofCustomTest"))]
      res <- sapply(copula, function(cop) {
        cat(paste(cop, "copula\n"))
        tmp <- suppressWarnings(.gofHybrid(x = x, copula = cop, tests = tests, customTests = customTests, dispstr = dispstr, margins = margins, M = M, MJ = MJ, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, lower = lower, upper = upper, seed.active = seed.active, processes = processes))
        cat("\n")
        return(tmp)
      }, USE.NAMES = FALSE)
      return(structure(
        class = "gofCOP",
        res
      ))
    } else if (priority == "tests") {
      tests <- gofTest4Copula("normal", dim(x)[2])
      tests <- tests[-c(which(tests == "gofHybrid"), which(tests == "gofCustomTest"))]
      copula_list <- lapply(tests, gofCopula4Test)
      copula <- Reduce(intersect, copula_list)
      res <- sapply(copula, function(cop) {
        cat(paste(cop, "copula\n"))
        tmp <- suppressWarnings(.gofHybrid(x = x, copula = cop, tests = tests, customTests = customTests, dispstr = dispstr, margins = margins, M = M, MJ = MJ, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, lower = lower, upper = upper, seed.active = seed.active, processes = processes))
        cat("\n")
        return(tmp)
      }, USE.NAMES = FALSE)
      return(structure(
        class = "gofCOP",
        res
      ))
    } else if (priority == "copula") {
      copula <- c("normal", "t", "clayton", "gumbel", "frank")
      tests_list <- lapply(copula, gofTest4Copula, d = dim(x)[2])
      tests <- Reduce(intersect, tests_list)
      tests <- tests[-c(which(tests == "gofHybrid"), which(tests == "gofCustomTest"))]
      res <- sapply(copula, function(cop) {
        cat(paste(cop, "copula\n"))
        tmp <- suppressWarnings(.gofHybrid(x = x, copula = cop, tests = tests, customTests = customTests, dispstr = dispstr, margins = margins, M = M, MJ = MJ, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m = m, lower = lower, upper = upper, seed.active = seed.active, processes = processes))
        cat("\n")
        return(tmp)
      }, USE.NAMES = FALSE)
      return(structure(
        class = "gofCOP",
        res
      ))
    }
  }
}
