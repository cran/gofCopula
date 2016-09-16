.SnK = function(x, cop) {
  n = dim(x)[1]
  dims = dim(x)[2]
  Cn = C.n(x,x)
  Kn = C.n(as.matrix(Cn), as.matrix(Cn))
  if (is.element(class(cop), c("claytonCopula", "frankCopula", "gumbelCopula"))){
    tcop = seq(1, n)/n
    texe = sapply(Kn[-n], function(x,y) length(which(x <= y)), tcop)
    Scvm1 = sum(texe[-n]^2 * (apply(matrix(tcop[-1], nrow = (n-1), ncol = dims),1,pCopula,cop) - apply(matrix(tcop[-n], nrow = (n-1), ncol = dims),1,pCopula,cop)))/n
    Scvm2 = sum(texe[-n] * (apply(matrix(tcop[-1], nrow = (n-1), ncol = dims),1,pCopula,cop)^2 - apply(matrix(tcop[-n], nrow = (n-1), ncol = dims),1,pCopula,cop)^2))
    SnK = n/3 + Scvm1 + Scvm2 
  } else if (is.element(class(cop), c("normalCopula", "tCopula"))) {
    Csample = rCopula(n*2, cop)
    Csample.n = C.n(Csample,Csample)
    Ksample = C.n(as.matrix(Csample.n), as.matrix(Csample.n))
    SnK = n * mean((Kn - Ksample)^2)
  }
  return(SnK)
}
