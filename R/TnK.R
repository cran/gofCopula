.TnK = function(x, cop) {
  n = dim(x)[1]
  dims = dim(x)[2]
  Cn = C.n(x,x)
  Kn = C.n(as.matrix(Cn), as.matrix(Cn))
  if (is.element(class(cop), c("claytonCopula", "frankCopula", "gumbelCopula"))){
    tcop1 = seq(0, n)/n
    tcop2 = (seq(0, n) + 1)/n
    Tn1 = abs((sapply(Kn, function(x,y) length(which(x <= y)), tcop1[-(n + 1)]) - apply(matrix(tcop1[-(n + 1)], nrow = n, ncol = dims),1,pCopula,cop)))
    Tn2 = abs((sapply(Kn, function(x,y) length(which(x <= y)), tcop1[-(n + 1)]) - apply(matrix(tcop2[-(n + 1)], nrow = n, ncol = dims),1,pCopula,cop)))
    TnK = sqrt(n) * max(Tn1, Tn2)
  } else if (is.element(class(cop), c("normalCopula", "tCopula"))) {
    Csample = rCopula(n*2, cop)
    Csample.n = C.n(Csample,Csample)
    Ksample = C.n(as.matrix(Csample.n), as.matrix(Csample.n))
    TnK = sqrt(n) * mean(abs(Kn - Ksample))
  }
  return(TnK)
}
