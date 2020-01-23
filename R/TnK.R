.TnK = function(x, cop) {
  n = dim(x)[1]
  dims = dim(x)[2]
  Cn = F.n(x, x)
  Kn = F.n(as.matrix(seq(1, n)/n), as.matrix(Cn))
  Csample.n = cop
  Ksample = F.n(as.matrix(seq(1, n)/n), as.matrix(Csample.n))
  TnK = sqrt(n) * max(abs(c(Kn[-n] - Ksample[-n], Kn[-n] - Ksample[-1])))
  return(TnK)
}
