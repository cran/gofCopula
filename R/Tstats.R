.Tstats = function(x, Tstat = c(eval(formals(gofTstat)$method), "Rn", "Tn", "Kernel", "SnK", "TnK", "White"), copula, ...) {
  if (is.element(Tstat, eval(formals(gofTstat)$method))) {
    gofTstat(u = x, method = Tstat, copula = copula)
  } else if (Tstat == "Rn") {
    .Rn(x = x, copula = copula, dims = copula@dimension)
  } else if (Tstat == "Tn") {
    add.parameters = list(...)$add.parameters
    .Tn(x = x, copula = copula, B = add.parameters[[1]], m = add.parameters[[2]], dims = copula@dimension, param.est = add.parameters[[3]], estim.method = add.parameters[[4]])
  } else if (Tstat == "Kernel") {
    add.parameters = list(...)$add.parameters
    .Kernel(x = x, copula = copula, dims = copula@dimension, n = nrow(x), nodes.Integration = add.parameters[[1]], MJ = add.parameters[[2]], delta.J = add.parameters[[3]])
  } else if (Tstat == "SnK") {
    .SnK(x = x, cop = copula)
  } else if (Tstat == "TnK") {
    .TnK(x = x, cop = copula)
  } else if (Tstat == "White") {
      add.parameters = list(...)$add.parameters
      BiCopGofTest(u1 = x[,1], u2 = x[,2], family = add.parameters[[1]], par = copula@parameters, par2 = if(class(copula) == "tCopula"){copula@df}else{0}, method = "white", B = 0)$statistic
  }
}
