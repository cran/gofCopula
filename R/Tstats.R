.Tstats <- function(x, Tstat, copula, ...) {
  if (is.element(Tstat, eval(formals(gofTstat)$method))) {
    res.f <- gofTstat(u = x, method = Tstat, copula = copula)
    return(res.f)
  } else if (Tstat == "Rn") {
    res.f <- .Rn(x = x, copula = copula, dims = copula@dimension)
    return(res.f)
  } else if (Tstat == "Tn") {
    add.parameters <- list(...)$add.parameters
    res.f <- .Tn(x = x, copula = copula, B = add.parameters[[1]], m = add.parameters[[2]], dims = copula@dimension, param.est = add.parameters[[3]], estim.method = add.parameters[[4]])
    return(res.f)
  } else if (Tstat == "Kernel") {
    add.parameters <- list(...)$add.parameters
    res.f <- .Kernel(x = x, copula = copula, dims = copula@dimension, n = nrow(x), nodes.Integration = add.parameters[[1]], MJ = add.parameters[[2]], delta.J = add.parameters[[3]])
    return(res.f)
  } else if (Tstat == "SnK") {
    cop.compare <- list(...)$cop.compare
    res.f <- .SnK(x = x, cop = cop.compare)
    return(res.f)
  } else if (Tstat == "TnK") {
    cop.compare <- list(...)$cop.compare
    res.f <- .TnK(x = x, cop = cop.compare)
    return(res.f)
  } else if (Tstat == "White") {
    add.parameters <- list(...)$add.parameters
    BiCopGofTest(u1 = x[, 1], u2 = x[, 2], family = add.parameters[[1]], par = if (inherits(copula, "tCopula")) {
      copula@parameters[-length(copula@parameters)]
    } else {
      copula@parameters
    }, par2 = if (inherits(copula, "tCopula")) {
      copula@parameters[length(copula@parameters)]
    } else {
      0
    }, method = "white", B = 0)$statistic
  } else if (!is.null(Tstat)) {
    res.f <- doCall(Tstat, x = x, copula = copula)
    return(res.f)
  }
}
