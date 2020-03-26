gofKendallCvM <- function(copula = c("normal", "t", "clayton", "gumbel", "frank"), x, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, dispstr = "ex", lower = NULL, upper = NULL, seed.active = NULL, processes = 1) {
  if (is.matrix(x) == FALSE) {
    stop("x must be a matrix")
  }
  if (length(copula) > 1) {
    stop("'copula' has to be a vector of length 1. Please select only one copula.")
  }
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == FALSE) {
    stop("This copula is not implemented for gofKendallCvM.")
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
  if (!inherits(param.est, "logical")) {
    stop("The argument 'param.est' has to be either 'TRUE' or 'FALSE'.")
  }
  if (copula == "gumbel" & param.est == FALSE & param <= 1) {
    param <- 1.5
    warning("When copula is 'gumbel', 'param' has to be larger 1. Because 'param.est' was set to 'FALSE', 'param' was set to 1.5 as default value for 'gumbel' copula.")
  }
  if (!inherits(df.est, "logical")) {
    stop("The argument 'df.est' has to be either 'TRUE' or 'FALSE'.")
  }
  if (!is.null(seed.active) & length(seed.active) != 1 & length(seed.active) != (M + 1)) {
    stop("The seed has to be an integer or a vector of M+1 seeds.")
  }
  if (!is.null(seed.active) & length(seed.active) == 1) {
    set.seed(seed.active)
    RNGkind(sample.kind = "default")
    seed.active <- sample(x = 2147483647, size = M + 1)
  }
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {
    stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")
  }

  erg <- .margins.param.est(copula = copula, margins = margins, x = x, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr, lower = lower, upper = upper)
  res <- try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnK", estim.method = "mpl", processes = processes, param.est = param.est, df.est = erg[[5]], dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active, lower = lower, upper = upper), silent = TRUE)
  if (inherits(res, "try-error")) {
    warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau. Therefore df.est was set to FALSE for the bootstrapping.")
    res.f <- .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "SnK", estim.method = "itau", processes = processes, param.est = param.est, df.est = FALSE, dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active, lower = lower, upper = upper)
    return(res.f)
  } else {
    return(res)
  }
}
