gofCheckTime = function(copula, x, test, M, MJ, print.res = T, margins = "ranks", dispstr = "ex", param = 0.5, param.est = T, df = 4, df.est = T, m = 1, delta.J = 0.5, nodes.Integration = 12, m_b = 0.5, zeta.m = 0, b_Rn = 0.05, processes = 1){
  lasted.time = c()
  N = c(2,5,10,15)
  NJ = c(2,5,10,15)
  for (j in 1:length(test)){
    times.comp = c()
    if (test[j] == "gofKernel"){
      options(show.error.messages = FALSE)
      for (i in N){
        for (ii in NJ){
          times.comp = rbind(times.comp, system.time(suppressWarnings(gofHybrid(copula = copula, x = x, dispstr = dispstr, testset = test[j], M = i, MJ = ii, execute.times.comp = F, margins = margins, param = param, param.est = param.est, df = df, df.est = df.est, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)))[3])
        }
      }
      times.comp = cbind(sort(rep(N, length(N))), rep(NJ, length(NJ)), times.comp)
      times.comp.mult = times.comp[, 2]*times.comp[, 1]
      times.lm = lm(times.comp[, 3] ~ times.comp.mult)
      lasted.time[j] = round(times.lm$coefficients[1] + times.lm$coefficients[2] * MJ * M)
      options(show.error.messages = TRUE)
    } else {
      options(show.error.messages = FALSE)
      for (i in N){
        times.comp = rbind(times.comp, system.time(suppressWarnings(gofHybrid(copula = copula, x = x, dispstr = dispstr, testset = test[j], M = i, execute.times.comp = F, margins = margins, param = param, param.est = param.est, df = df, df.est = df.est, MJ = MJ, delta.J = delta.J, nodes.Integration = nodes.Integration, m_b = m_b, zeta.m = zeta.m, b_Rn = b_Rn, m = m, processes = processes)))[3])
      }
      times.comp = cbind(N, times.comp)
      times.lm = lm(times.comp[, 2] ~ times.comp[, 1])
      lasted.time[j] = round(times.lm$coefficients[1] + times.lm$coefficients[2] * M)
      options(show.error.messages = TRUE)
    }
  }
  if (print.res == T){
    .get.time(sum(lasted.time))
  } else {
    store.time = sum(lasted.time)
  }
}
