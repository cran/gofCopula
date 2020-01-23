gofKendallKS = function(copula = c("normal", "t", "clayton", "gumbel", "frank"), x, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, dispstr = "ex", seed.active = NULL, processes = 1){
  if (is.matrix(x) == FALSE){stop("x must be a matrix")}
  dims = dim(x)[2]
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == FALSE){stop("This copula is not implemented for gofKendallKS.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}  
  if (!is.numeric(processes)) {stop("The argument 'processes' has to be a numeric.")}
  if (processes %% 1 != 0 | processes < 1) {stop("The argument 'processes' has to be a positiv integer.")}
  if (!is.numeric(M)) {stop("The argument 'M' has to be a numeric.")}
  if (M %% 1 != 0 | M < 0) {stop("The argument 'M' has to be a positiv integer.")}
  if (!is.numeric(param)) {stop("The argument 'param' has to be a numeric.")}
  if (!is.numeric(df)) {stop("The argument 'df' has to be a numeric.")}
  if (!inherits(param.est, "logical")) {stop("The argument 'param.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!inherits(df.est, "logical")) {stop("The argument 'df.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!is.null(seed.active) & length(seed.active) != 1 & length(seed.active) != (M+1)) {stop("The seed has to be an integer or a vector of M+1 seeds.")}
  if (!is.null(seed.active) & length(seed.active) == 1) {set.seed(seed.active); RNGkind(sample.kind = "default"); seed.active = sample(x = 2147483647, size = M+1)}
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")}

  erg = .margins.param.est(copula = copula, margins = margins, x = x, param = param, param.est = param.est, df = df, df.est = df.est, dispstr = dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "TnK", estim.method = "mpl", processes = processes, param.est = param.est, df.est = erg[[5]], dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active), silent = TRUE)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau. Therefore df.est was set to FALSE for the bootstrapping."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "TnK", estim.method = "itau", processes = processes, param.est = param.est, df.est = FALSE, dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active)} else {res}
}
