.gofCopulapb <- function(copula, x, M, method, estim.method, processes, param.est, df.est, dispstr, param.margins, margins, seed.active, lower, upper, ...) {
  # required to avoid error on global variable definition
  cnt <- NULL

  if (!is.null(seed.active)) {
    set.seed(seed.active[1], kind = "default")
  }
  bs.ac.c <- list()

  if (is.element(method, c("SnB", "SnC", "AnChisq", "AnGamma"))) {
    x <- do.call(cCopula, c(list(x, copula = copula)))
  }
  ac.c <- if (method == "Sn") {
    if (inherits(copula, "tCopula")) {
      copula@parameters[length(copula@parameters)] <- as.integer(copula@parameters[length(copula@parameters)])
    }
    .Tstats(x, Tstat = method, copula = copula)
  } else if (method == "Tn") {
    if (inherits(copula, "tCopula")) {
      copula@parameters[2] <- min(copula@parameters[2], 60)
    }
    add.parameters <- list(...)$add.parameters
    .Tstats(x, Tstat = method, copula = copula, add.parameters = add.parameters)
  } else if (method == "Kernel" || method == "White") {
    add.parameters <- list(...)$add.parameters
    .Tstats(x, Tstat = method, copula = copula, add.parameters = add.parameters)
  } else {
    cop.compare <- rCopula(10000, copula)
    cop.compare.n <- F.n(cop.compare, cop.compare)
    .Tstats(x, Tstat = method, copula = copula, cop.compare = cop.compare.n)
  }

  if (processes > 1) {
    cl <- makeCluster(processes)
    clusterEvalQ(cl, library(copula))
    clusterEvalQ(cl, library(foreach))
    clusterEvalQ(cl, library(progress))
    clusterEvalQ(cl, library(gofCopula))
    if (method == "White") {
      clusterEvalQ(cl, "BiCopDeriv")
      clusterEvalQ(cl, "BiCopDeriv2")
      clusterEvalQ(cl, "BiCopPDF")
      clusterEvalQ(cl, "ginv")
    }
    if (!is.element(method, c("SnB", "SnC", "AnChisq", "AnGamma", "Sn", "Tn", "Kernel", "White", "SnK", "TnK", "Rn"))) {
      clusterExport(cl, paste(method))
    }
    registerDoSNOW(cl)
  } else {
    registerDoSEQ()
  }

  pb <- progress_bar$new(total = M, format = ":dummy [:bar] :percent | time left: :eta", force = TRUE, show_after = 0)
  progBar_dummy <- rep("Progress:", M)
  progBar <- function(n) {
    pb$tick(tokens = list(dummy = progBar_dummy[n]))
  }
  opts <- list(progress = progBar)

  bs.ac.c <- foreach(cnt = seq_len(M), .options.snow = opts) %dopar% {
    progBar(cnt)
    seed.counter <- 0
    repeat {
      if (!is.null(seed.active)) {
        new.seed <- seed.active[cnt + 1] + (M + 1) * seed.counter
        while (is.element(new.seed, seed.active)) {
          seed.counter <- seed.counter + 1
          new.seed <- seed.active[cnt + 1] + (M + 1) * seed.counter
        }
        seed.active <- unique(c(seed.active, new.seed))
        set.seed(new.seed, kind = "default")
      }
      xsim <- rCopula(nrow(x), copula)
      # margins entry has to be NULL, since the data are in [0,1] already
      copula.sim <- try(.margins.param.est(
        copula = as.character(substr(class(copula), 1, nchar(class(copula)) - 6)), margins = NULL, x = xsim,
        param = if (!inherits(copula, "tCopula")) {
          copula@parameters
        } else {
          copula@parameters[-length(copula@parameters)]
        },
        param.est = param.est,
        df = if (!inherits(copula, "tCopula")) {
          4
        } else {
          copula@parameters[length(copula@parameters)]
        },
        df.est = df.est,
        dispstr = dispstr,
        lower = lower,
        upper = upper
      )[[1]], silent = TRUE)

      if (!inherits(copula.sim, "try-error")) {
        break
      }
      seed.counter <- seed.counter + 1
    }
    if (is.element(method, c("SnB", "SnC", "AnChisq", "AnGamma"))) {
      xsim <- do.call(cCopula, c(list(xsim, copula = copula.sim)))
    }
    if (method == "Sn") {
      if (inherits(copula.sim, "tCopula")) {
        copula.sim@parameters[length(copula.sim@parameters)] <- as.integer(copula.sim@parameters[length(copula.sim@parameters)])
      }
      .Tstats(xsim, Tstat = method, copula = copula.sim)
    } else if (method == "Tn") {
      if (inherits(copula, "tCopula")) {
        copula.sim@parameters[2] <- min(copula.sim@parameters[2], 60)
      }
      .Tstats(xsim, Tstat = method, copula = copula.sim, add.parameters = add.parameters)
    } else if (method == "Kernel" || method == "White") {
      .Tstats(xsim, Tstat = method, copula = copula.sim, add.parameters = add.parameters)
    } else {
      .Tstats(xsim, Tstat = method, copula = copula.sim, cop.compare = cop.compare.n)
    }
  }

  bs.ac.c <- unlist(bs.ac.c)

  if (processes > 1) {
    stopCluster(cl)
  }

  ac.c <- as.numeric(ac.c)
  bs.ac.c <- as.numeric(bs.ac.c)
  test <- sum(abs(bs.ac.c) >= abs(ac.c)) / M

  switch(method,
    SnB = {
      matrix_names <- "RosenblattSnB"
    },
    SnC = {
      matrix_names <- "RosenblattSnC"
    },
    AnChisq = {
      matrix_names <- "RosenblattChisq"
    },
    AnGamma = {
      matrix_names <- "RosenblattGamma"
    },
    Sn = {
      matrix_names <- "Sn"
    },
    Rn = {
      matrix_names <- "PIOSRn"
    },
    Tn = {
      matrix_names <- "PIOSTn"
    },
    Kernel = {
      matrix_names <- "Kernel"
    },
    SnK = {
      matrix_names <- "KendallCvM"
    },
    TnK = {
      matrix_names <- "KendallKS"
    },
    White = {
      matrix_names <- "White"
    },
    {
      matrix_names <- method
    }
  )
  res <- structure(
    class = "gofCOP",
    list(
      list(
        method = sprintf(
          "Parametric bootstrap goodness-of-fit test with %s test and %s copula",
          matrix_names, as.character(substr(class(copula), 1, nchar(class(copula)) - 6))
        ),
        copula = as.character(substr(class(copula), 1, nchar(class(copula)) - 6)),
        margins = margins,
        param.margins = param.margins,
        theta = if (!inherits(copula, "tCopula")) {
          copula@parameters
        } else {
          copula@parameters[-length(copula@parameters)]
        },
        df = if (!inherits(copula, "tCopula")) {
          NULL
        } else {
          copula@parameters[length(copula@parameters)]
        },
        res.tests = matrix(c(test, ac.c),
          ncol = 2,
          dimnames = list(matrix_names, c("p.value", "test statistic"))
        )
      )
    )
  )
  names(res) <- substr(class(copula), 1, nchar(class(copula)) - 6)
  return(res)
}
