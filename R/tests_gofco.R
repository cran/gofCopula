gofco <- function(copulaobject, x, tests = c("gofPIOSRn", "gofKernel"), customTests = NULL, margins = "ranks",
                  M = 1000, MJ = 100, dispstr = "ex",
                  m = 1, delta.J = 0.5, nodes.Integration = 12,
                  lower = NULL, upper = NULL,
                  seed.active = NULL, processes = 1) {
  if (is.matrix(x) == FALSE) {
    stop("x must be a matrix")
  }
  if (!is.null(tests) & any(!is.element(tests, gofTest4Copula(as.character(substr(class(copulaobject), 1, nchar(class(copulaobject)) - 6)), dim(x)[2]))) || !is.null(customTests) & any(!is.element(customTests, ls(".GlobalEnv")))) {
    stop("At least one of the tests in 'tests' is not implemented, cannot handle a dataset of this dimension or at least one of the tests in 'customTests' does not match any function in the global workspace. Please check if it is correctly spelled in the function call.")
  }
  if (any(!vapply(customTests, function(x) all(names(formals(x)) %in% c("x", "copula")), TRUE))) {
    stop("At least one function in 'customTest' does not follow the requirements for the arguments. The first argument for the dataset has to be called 'x', the second one for the copula has to be called 'copula'.")
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
  switch(class(copulaobject),
    normalCopula = {
      copula <- "normal"
      param <- copulaobject@parameters
      param.est <- if (is.na(copulaobject@parameters)) {
        TRUE
      } else {
        FALSE
      }
      df <- 4
      df.est <- TRUE
      dispstr <- copulaobject@dispstr
    },
    tCopula = {
      copula <- "t"
      param <- copulaobject@parameters[-length(copulaobject@parameters)]
      param.est <- if (is.na(copulaobject@parameters[1])) {
        TRUE
      } else {
        FALSE
      }
      df <- copulaobject@parameters[length(copulaobject@parameters)]
      df.est <- if (copulaobject@df.fixed == FALSE) {
        TRUE
      } else if (copulaobject@df.fixed == TRUE) {
        FALSE
      }
      dispstr <- copulaobject@dispstr
    },
    claytonCopula = {
      copula <- "clayton"
      param <- copulaobject@parameters
      param.est <- if (is.na(copulaobject@parameters)) {
        TRUE
      } else {
        FALSE
      }
      df <- 4
      df.est <- TRUE
    },
    frankCopula = {
      copula <- "frank"
      param <- copulaobject@parameters
      param.est <- if (is.na(copulaobject@parameters)) {
        TRUE
      } else {
        FALSE
      }
      df <- 4
      df.est <- TRUE
    },
    gumbelCopula = {
      copula <- "gumbel"
      param <- copulaobject@parameters
      param.est <- if (is.na(copulaobject@parameters)) {
        TRUE
      } else {
        FALSE
      }
      df <- 4
      df.est <- TRUE
    },
    stop("The class of the object is not supported.")
  )
  res.f <- .gofHybrid(copula = copula, x = x, tests = tests, customTests = customTests, margins = margins, dispstr = dispstr, M = M, param = param, param.est = param.est, df = df, df.est = df.est, m = m, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, lower = lower, upper = upper, seed.active = seed.active, processes = processes)
  return(res.f)
}
