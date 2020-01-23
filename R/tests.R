gofPIOSRn = function(copula = c("normal", "t", "clayton", "gumbel", "frank"), x, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, dispstr = "ex", seed.active = NULL, processes = 1){
  if (is.matrix(x) == FALSE){stop("x must be a matrix")}
  if (dim(x)[2] < 2 || dim(x)[2] > 3){stop("x must be of dimension 2 or 3")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == FALSE){stop("This copula is not implemented for gofPIOSRn.")}
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
  if (!is.null(seed.active) & length(seed.active) == 1) {set.seed(seed.active); RNGkind(sample.kind = "default"); seed.active = sample(x=2147483647, size = M+1)}
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")}
  
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Rn", estim.method = "mpl", processes = processes, param.est = param.est, df.est = erg[[5]], dispstr=dispstr, param.margins=erg[[4]], margins = margins, seed.active = seed.active), silent = TRUE)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau. Therefore df.est was set to FALSE for the bootstrapping."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Rn", estim.method = "itau", processes = processes, param.est = param.est, df.est = FALSE, dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active)} else {res}
}

#################################################################################
gofPIOSTn = function(copula = c("normal", "t", "clayton", "gumbel", "frank"), x, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, dispstr = "ex", m = 1, seed.active = NULL, processes = 1){
  if (is.matrix(x) == FALSE){stop("x must be a matrix")}
  if (dim(x)[2] < 2 || dim(x)[2] > 3){stop("x must be of dimension 2 or 3")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == FALSE){stop("This copula is not implemented for gofPIOSTn.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
  if (!is.numeric(processes)) {stop("The argument 'processes' has to be a numeric.")}
  if (processes %% 1 != 0 | processes < 1) {stop("The argument 'processes' has to be a positiv integer.")}
  if (!is.numeric(M)) {stop("The argument 'M' has to be a numeric.")}
  if (M %% 1 != 0 | M < 0) {stop("The argument 'M' has to be a positiv integer.")}
  if (!is.numeric(param)) {stop("The argument 'param' has to be a numeric.")}
  if (!is.numeric(df)) {stop("The argument 'df' has to be a numeric.")}
  if (!is.numeric(m)) {stop("The argument 'm' has to be a numeric.")}
  n = dim(x)[1]
  if (n %% m != 0 | m < 1) {stop("The length of the blocks, 'm', has to be larger 1 and a divisor of the length of the data sequence.")}
  if (!inherits(param.est, "logical")) {stop("The argument 'param.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!inherits(df.est, "logical")) {stop("The argument 'df.est' has to be either 'TRUE' or 'FALSE'.")}
  B = n / m
  if (!is.null(seed.active) & length(seed.active) != 1 & length(seed.active) != (M+1)) {stop("The seed has to be an integer or a vector of M+1 seeds.")}
  if (!is.null(seed.active) & length(seed.active) == 1) {set.seed(seed.active); RNGkind(sample.kind = "default"); seed.active = sample(x=2147483647, size = M+1)}
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")}
  
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  add.parameters = list(B, m, param.est, "mpl")
  if (copula == "t") {erg[[1]]@parameters[2] = min(erg[[1]]@parameters[2], 60)}
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Tn", estim.method = "mpl", processes = processes, add.parameters = add.parameters, param.est = param.est, df.est = erg[[5]], dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active), silent = TRUE)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau. Therefore df.est was set to FALSE for the bootstrapping."); add.parameters = list(B, m, param.est, "itau"); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Tn", estim.method = "itau", processes = processes, add.parameters = add.parameters, param.est = param.est, df.est = FALSE, dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active)} else {res}
}

#################################################################################
gofKernel = function(copula = c("normal", "t", "clayton", "gumbel", "frank"), x, param = 0.5, param.est = TRUE, df = 4, df.est = TRUE, margins = "ranks", M = 1000, MJ = 100, dispstr = "ex", delta.J = 0.5, nodes.Integration = 12, seed.active = NULL, processes = 1){
  if (is.matrix(x) == FALSE){stop("x must be a matrix")}
  if (dim(x)[2] != 2){stop("x must be of dimension 2")}
  if (is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel")) == FALSE){stop("This copula is not implemented for gofKernel.")}
  if (!is.element(dispstr, c("ex", "un"))) {stop("dispstr has to be either 'ex' or 'un'. See documentation for more information.")}
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
  if (!inherits(param.est, "logical")) {stop("The argument 'param.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!inherits(df.est, "logical")) {stop("The argument 'df.est' has to be either 'TRUE' or 'FALSE'.")}
  if (!is.null(seed.active) & length(seed.active) != 1 & length(seed.active) != (M+1)) {stop("The seed has to be an integer or a vector of M+1 seeds.")}
  if (!is.null(seed.active) & length(seed.active) == 1) {set.seed(seed.active); RNGkind(sample.kind = "default"); seed.active = sample(x=2147483647, size = M+1)}
  if (!is.null(seed.active) & all(!vapply(seed.active, function(x) x %% 1 == 0, TRUE))) {stop("All seeds have to be whole numbers. Please check seed.active for non-whole numbers.")}
  add.parameters = list(nodes.Integration, MJ, delta.J)
  
  erg = .margins.param.est(copula=copula, margins=margins, x=x, param=param, param.est=param.est, df=df, df.est=df.est, dispstr=dispstr)
  res = try(.gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Kernel", estim.method = "mpl", processes = processes, add.parameters = add.parameters, param.est = param.est, df.est = erg[[5]], dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active), silent = TRUE)
  if (class(res) == "try-error"){warning("Pseudo Maximum Likelihood estimation of the parameters while the bootstrapping procedure failed. The estimation was performed with inversion of Kendall's Tau. Therefore df.est was set to FALSE for the bootstrapping."); .gofCopulapb(copula = erg[[1]], x = erg[[2]], M = M, method = "Kernel", estim.method = "itau", processes = processes, add.parameters = add.parameters, param.est = param.est, df.est = FALSE, dispstr = dispstr, param.margins = erg[[4]], margins = margins, seed.active = seed.active)} else {res}
}
