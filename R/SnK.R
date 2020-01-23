.SnK = function(x, cop) {
  n = dim(x)[1]
  dims = dim(x)[2]
  Cn = F.n(x,x)
  Kn = F.n(as.matrix(seq(1, n)/n), as.matrix(Cn))
  tcop = seq(1, n)/n
  Csample.n = cop
  Scvm1 = sum(Kn[-n]^2 * (F.n(as.matrix(tcop[-1]), as.matrix(Csample.n)) - F.n(as.matrix(tcop[-n]), as.matrix(Csample.n))))#/n
  Scvm2 = sum(Kn[-n] * (F.n(as.matrix(tcop[-1]), as.matrix(Csample.n))^2 - F.n(as.matrix(tcop[-n]), as.matrix(Csample.n))^2))
  SnK = n/3 + n * Scvm1 - n * Scvm2 
  return(SnK)
}
