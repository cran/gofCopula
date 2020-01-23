.Rn = function(x, copula, dims) {
  switch(dims-1, {
  switch(class(copula),
         claytonCopula = {
           Rn.ac.c = -sum(.V.clay.12(copula@parameters, x[, 1], x[, 2])^2) / sum(.S.clay.12(copula@parameters, x[, 1], x[, 2])) - 1
           Rn.ac.c
         },
         gumbelCopula = {
           Rn.ac.g = -sum(.V.gumb.12(copula@parameters, x[, 1], x[, 2])^2) / sum(.S.gumb.12(copula@parameters, x[, 1], x[, 2])) - 1
           Rn.ac.g
         },
         frankCopula = {
           Rn.ac.f = -sum(.V.fran.12(copula@parameters, x[, 1], x[, 2])^2) / sum(.S.fran.12(copula@parameters, x[, 1], x[, 2])) - 1
           Rn.ac.f
         },
         tCopula = {
           psn.sample = qt(x, df = tail(copula@parameters, n = 1))
           Rn.t = -sum(apply(psn.sample, 1, FUN = .V.t, rho = copula@parameters[-length(copula@parameters)], nu = tail(copula@parameters, n = 1))^2)/sum(apply(psn.sample, 1, FUN = .S.t, rho = copula@parameters[-length(copula@parameters)], nu = tail(copula@parameters, n = 1))) - 1
           Rn.t
         },
         normalCopula = {
           sig = matrix(c(1, copula@parameters, copula@parameters, 1), ncol = 2, byrow = TRUE)
           sig.inv = rbind(c(sig[1, 1], -sig[1, 2]), c(-sig[2, 1], sig[2, 2]))/det(sig)
           psn.sample = qnorm(x)
           Rn.g = -sum(apply(psn.sample, 1, FUN = .V.ga, sig.inv = sig.inv, dims = dims)[2,]^2)/sum(apply(psn.sample, 1, FUN = .S.ga.12, sig = sig)) - 1
           Rn.g
         })
  },
  {
    switch(class(copula),
           claytonCopula = {
             Rn.ac.c = -sum(.clay.V.123(copula@parameters, x)^2)/sum(.clay.S.123(copula@parameters, x)) - 1
             Rn.ac.c
           },
           gumbelCopula = {
             Rn.ac.g = -sum(.gumb.V.123(copula@parameters, x)^2)/sum(.gumb.S.123(copula@parameters, x)) - 1
             Rn.ac.g
           },
           frankCopula = {
             Rn.ac.f = -sum(.fran.V.123(copula@parameters, x)^2)/sum(.fran.S.123(copula@parameters, x)) - 1
             Rn.ac.f
           },
           tCopula = {
             stop("The t-Copula is not implemented for dimension 3 and the test gofPIOSRn.")
           },
           normalCopula = {
             sig = getSigma(copula)
             sig.inv = solve(sig) 
             sig.inv = (sig.inv + t(sig.inv)) / 2
             psn.sample = qnorm(x)
             V = t(apply(psn.sample, 1, FUN = .V.ga, sig.inv = sig.inv, dims = dims)[c(2, 3, 6),])
             V = matrix(colSums(cbind(V[, 1]^2, V[, 1] * V[, 2], V[, 1] * V[, 3], V[, 1] * V[, 2], V[, 2]^2, V[, 2] * V[, 3], V[, 1] * V[, 3], V[, 2] * V[, 3], V[, 3]^2)), ncol = dims)
             S1212 = apply(psn.sample, 1, FUN = .S.ga.12.12, sig = sig)
             S1213 = apply(psn.sample, 1, FUN = .S.ga.12.13, sig = sig)
             S1223 = apply(psn.sample, 1, FUN = .S.ga.12.23, sig = sig)
             S1313 = apply(psn.sample, 1, FUN = .S.ga.13.13, sig = sig)
             S1323 = apply(psn.sample, 1, FUN = .S.ga.13.23, sig = sig)
             S2323 = apply(psn.sample, 1, FUN = .S.ga.23.23, sig = sig)
             S = matrix(colSums(cbind(S1212, S1213, S1223, S1213, S1313, S1323, S1223, S1323, S2323)), ncol = dims)
             Rn.g = sum(diag(-solve(S) %*% V)) - 3
             Rn.g
           })
  })
}

.Tn = function(x, copula, B, m, dims, param.est, estim.method) {
  switch(dims - 1, {
  switch(class(copula),
         claytonCopula = {
           l.ac.c = log(.clay.12.density(copula@parameters, x))
           l.ac.c.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.ac.c.b = c(l.ac.c.b, .clay.12.density(fitCopula(claytonCopula(dim = dims), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, x[((b - 1) * m + 1):(b * m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.c.b = c(l.ac.c.b, .clay.12.density(copula@parameters, x[((b - 1) * m + 1):(b * m),]))
             }
           }
           l.ac.c.b = l.ac.c.b[-1]
           stat = sum(l.ac.c - log(l.ac.c.b)) - 1
           stat
         },
         gumbelCopula = {
           l.ac.g = log(.gumb.12.density(copula@parameters, x))
           l.ac.g.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.ac.g.b = c(l.ac.g.b, .gumb.12.density(fitCopula(gumbelCopula(dim = dims), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, x[((b - 1) * m + 1):(b * m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.g.b = c(l.ac.g.b, .gumb.12.density(copula@parameters, x[((b - 1) * m + 1):(b * m),]))
             }
           }
           
           l.ac.g.b = l.ac.g.b[-1]
           stat = sum(l.ac.g - log(l.ac.g.b)) - 1
           stat
         },
         frankCopula = {
           l.ac.f = log(.fran.12.density(copula@parameters, x))
           l.ac.f.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.ac.f.b = c(l.ac.f.b, .fran.12.density(fitCopula(frankCopula(dim = dims), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, x[((b - 1) * m + 1):(b * m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.f.b = c(l.ac.f.b, .fran.12.density(copula@parameters, x[((b - 1) * m + 1):(b * m),]))
             }
           }
           l.ac.f.b = l.ac.f.b[-1]
           stat = sum(l.ac.f - log(l.ac.f.b)) - 1
           stat
         },
         tCopula = {
           psn.sample = qt(x, df = tail(copula@parameters, n = 1))
           l.t = log(.t.12.dens(copula@parameters[-length(copula@parameters)], psn.sample, nu = tail(copula@parameters, n = 1)))
           l.t.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.t.b = c(l.t.b, .t.12.dens(fitCopula(tCopula(dim = dims, df = tail(copula@parameters, n = 1), df.fixed = TRUE), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, psn.sample[((b - 1) * m + 1):(b * m),], nu = tail(copula@parameters, n = 1)))
             }
           } else {
             for(b in 1:B){
               l.t.b = c(l.t.b, .t.12.dens(copula@parameters[-length(copula@parameters)], psn.sample[((b - 1) * m + 1):(b * m),], nu = tail(copula@parameters, n = 1)))
             }
           }
           l.t.b = l.t.b[-1]
           stat = sum(l.t - log(l.t.b)) - 1
           stat
         },
         normalCopula = {
           psn.sample = qnorm(x)
           l.ac.ga.b = as.vector(sapply(1:B, FUN = .opt.ga, psn.sample = psn.sample, m = m, dims = dims))
           stat = sum(log(dCopula(x, normalCopula(copula@parameters, dims))) - log(l.ac.ga.b)) - 1 
           stat
         })
  },
{
  switch(class(copula),
         claytonCopula = {
           l.ac.c = log(.clay.123.density(copula@parameters, x))
           l.ac.c.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.ac.c.b = c(l.ac.c.b, .clay.123.density(fitCopula(claytonCopula(dim = dims), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, x[((b - 1) * m + 1):(b * m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.c.b = c(l.ac.c.b, .clay.123.density(copula@parameters, x[((b - 1) * m + 1):(b * m),]))
             }
           }
           l.ac.c.b = l.ac.c.b[-1]
           stat = sum(l.ac.c - log(l.ac.c.b)) - 1
           stat
         },
         gumbelCopula = {
           l.ac.g = log(.gumb.123.density(copula@parameters, x))
           l.ac.g.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.ac.g.b = c(l.ac.g.b, .gumb.123.density(fitCopula(gumbelCopula(dim = dims), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, x[((b - 1) * m + 1):(b * m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.g.b = c(l.ac.g.b, .gumb.123.density(copula@parameters, x[((b - 1) * m + 1):(b * m),]))
             }
           }
           
           l.ac.g.b = l.ac.g.b[-1]
           stat = sum(l.ac.g - log(l.ac.g.b)) - 1
           stat
         },
         frankCopula = {
           l.ac.f = log(.fran.123.density(copula@parameters, x))
           l.ac.f.b = 0
           if (param.est == TRUE){
             for(b in 1:B){
               l.ac.f.b = c(l.ac.f.b, .fran.123.density(fitCopula(frankCopula(dim = dims), data = x[-(((b - 1) * m + 1):(b * m)),], method = estim.method, estimate.variance = FALSE)@estimate, x[((b - 1) * m + 1):(b * m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.f.b = c(l.ac.f.b, .fran.123.density(copula@parameters, x[((b - 1) * m + 1):(b * m),]))
             }
           }
           l.ac.f.b = l.ac.f.b[-1]
           stat = sum(l.ac.f - log(l.ac.f.b)) - 1
           stat
         },
         tCopula = {
           stop("The t-Copula is not implemented for dimension 3 and the test gofPIOSTn.")
         },
         normalCopula = {
           psn.sample = qnorm(x)
           l.ac.ga.b = as.vector(sapply(1:B, FUN = .opt.ga, psn.sample = psn.sample, m = m, dims = dims))
           stat = sum(log(dCopula(x, normalCopula(copula@parameters, dims, dispstr = copula@dispstr))) - log(l.ac.ga.b)) - 3 
           stat
         })
})
}

.Kernel = function(x, copula, dims, n, nodes.Integration, MJ, delta.J) {
  switch(class(copula),
         claytonCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, claytonCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         },
         gumbelCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, gumbelCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         },
         frankCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, frankCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_f = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_f = c(Jn_f, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_f)
           stat
         },
         tCopula = {
           psn.sample = qt(x, df = tail(copula@parameters, n = 1))
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, tCopula(copula@parameters[-length(copula@parameters)], df = tail(copula@parameters, n = 1), dim = dims, dispstr = copula@dispstr))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         },
         normalCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, normalCopula(copula@parameters, dim = dims, dispstr = copula@dispstr))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         })
  
}