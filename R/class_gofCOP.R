print.gofCOP <- function(x, ...) {
  for (j in seq_along(x)) {
    cat(rep("-", getOption("width")), sep = "")
    cat("\n")
    cat(strwrap(x[[j]]$method), sep = "\n")
    cat("\n")
    out <- character()
    if (!is.null(x[[j]]$statistic)) {
      out <- c(out, paste("statistic =", format(signif(x[[j]]$statistic, 2))))
    }
    if (!is.null(x[[j]]$p.value)) {
      fp <- format.pval(x[[j]]$p.value, digits = 2)
      out <- c(out, paste("p-value =", fp))
    }
    if (!is.null(x[[j]]$theta)) {
      for (i in seq_along(x[[j]]$theta)) {
        out <- c(out, paste0("theta.", i, " = ", x[[j]]$theta[i]))
      }
    }
    if (!is.null(x[[j]]$df)) {
      out <- c(out, paste("df =", x[[j]]$df))
    }
    cat("Parameters:")
    cat("\n")
    cat(paste(out, collapse = "\n"))
    cat("\n")
    cat("\n")

    if (!is.null(x[[j]]$df)) {
      if (length(x[[j]]$df) > 1 & !any(x[[j]]$df == 60)) {
        cat("For the Sn test the df have to be an integer and for this reason were transformed.")
        cat("\n")
        cat("\n")
      }
      if (any(x[[j]]$df == 60)) {
        cat("For the PIOSTn test the estimated df were too high and for computational reasons were fixed to 60.")
        cat("\n")
        cat("\n")
      }
    }

    if (!is.null(x[[j]]$res.tests)) {
      cat("Tests results:")
      cat("\n")
      print(x[[j]]$res.tests, max = 50)
      if (dim(x[[j]]$res.tests)[1] != 1) {
        cat("\n")
        cat("Please use the functions gofGetHybrid() and gofOutputHybrid() for display of subsets of Hybrid tests. To access the results, please obtain them from the structure of the gofCOP object.")
      }
      if (any(is.na(x[[j]]$res.tests[, 1]))) {
        cat("\n")
        cat(bold("At least one p-value was returned as 'NA'. Typically this is due to the copula parameter being estimated close to the boundaries. Consider adjusting the possible parameter space via 'lower' and 'upper'."))
      }
    }
    cat("\n")
    cat("\n")
    if (is.element("White", rownames(x[[j]]$res.tests)) & x[[j]]$copula == "t") {
      cat("The test gofWhite may be unstable for t-copula. Please handle the results carefully.")
      cat("\n")
      cat("\n")
    }
  }
  invisible(x)
}

plot.gofCOP <- function(x, copula = NULL, hybrid = NULL, point.bg = "white", point.col = "black", point.cex = .7,
                        jitter.val = 0.1, pal = "xmen", bean.b.o = .2, inf.method = "hdi", theme = 2, ...) {
  if (inherits(x, "gofCOP") & dim(x[[1]]$res.tests)[1] == 1) {
    stop("Plotting p-values is only supported for various tests. Please consider gofHybrid() and gof().")
  }
  if (!is.element(inf.method, c("hdi", "iqr", "sd", "se"))) {
    stop("Please consider for the argument 'inf.method' either 'hdi', 'iqr', 'sd' or 'se'.")
  }

  # selection of copulae
  if (!is.null(copula)) {
    if (!all(is.element(copula, names(x)))) {
      stop("Please select only copulae included in your object of class gofCOP.")
    }
    cops_max <- c("normal", "t", "clayton", "frank", "gumbel")
    cops_out <- cops_max[!is.element(cops_max, copula)]
    for (i in seq_along(cops_out)) {
      eval(parse(text = paste("x$", cops_out[i], "= NULL", sep = "")))
    }
  }

  numb_tests <- sum(!grepl("hybrid", rownames(x[[1]]$res.tests)))
  if (is.null(hybrid)) {
    hybrid <- 1:numb_tests
  }

  for (i in seq_along(hybrid)) {
    if (!is.numeric(hybrid[i])) {
      stop("Please select for 'hybrid' a vector of integers defining which hybrid testing sizes are desired.")
    } else {
      if (floor(hybrid[i]) != hybrid[i]) {
        stop("Please select for 'hybrid' a vector of integers defining which hybrid testing sizes are desired.")
      } else {
        if ((hybrid[i] < 1) | (hybrid[i] > numb_tests)) {
          stop("The elements of the 'hybrid' argument should not be smaller than 1 or larger than the number of single tests.")
        }
      }
    }
  }

  hybrid <- sort(unique(hybrid))

  cops <- pvalues <- hybrids <- c()
  for (j in seq_along(x)) {
    numb_tests <- sum(!grepl("hybrid", rownames(x[[j]]$res.tests)))
    cops_tmp <- rep(x[[j]]$copula, each = numb_tests)
    pvalues_tmp <- c(gofOutputHybrid(x, nsets = 1)[[j]]$res.tests[, 1])
    hybrids_tmp <- rep(c("hyb1"), numb_tests)
    hybs_tmp <- c(c("hyb1"), paste("hyb", 2:numb_tests, sep = ""))
    for (i in seq_len(numb_tests)[-1]) {
      current.comb <- gofOutputHybrid(x, nsets = i)[[j]]$res.tests[-c(1:numb_tests), 1]
      cops_tmp <- c(cops_tmp, rep(x[[j]]$copula, each = length(current.comb)))
      pvalues_tmp <- c(pvalues_tmp, current.comb)
      hybrids_tmp <- c(hybrids_tmp, rep(hybs_tmp[i], length(current.comb) * length(x[[j]]$copula)))
    }
    cops <- c(cops, cops_tmp)
    pvalues <- c(pvalues, pvalues_tmp)
    hybrids <- c(hybrids, hybrids_tmp)
  }

  toplot <- as.data.frame(cbind(cops, pvalues, hybrids))
  toplot[, 2] <- as.numeric(as.character(toplot[, 2]))
  colnames(toplot) <- c("copulae", "pvalues", "hybrids")

  # selection of hybrids
  desired_hybs <- paste("hyb", hybrid, sep = "")
  toplot <- toplot[is.element(toplot$hybrids, desired_hybs), ]

  cops <- as.character(unique(toplot$copulae))

  # binary data case with iqr
  if (inf.method == "iqr") {
    for (i in seq_along(cops)) {
      for (j in seq_along(hybrid)) {
        index <- which((toplot$copulae == cops[i]) & (!is.na(toplot$pvalues)) & (as.character(toplot$hybrids) == paste("hyb", as.character(hybrid[j]), sep = "")))
        if ((length(setdiff(toplot[index, 2], c(0, 1))) == 0) & (length(index) > 0)) {
          toplot[index[1], 2] <- toplot[index[1], 2] + 0.000000000000001
        }
      }
    }
  }

  # sd/se case
  if (is.element(inf.method, c("sd", "se"))) {
    cnt <- c()
    for (i in seq_along(cops)) {
      for (j in seq_along(hybrid)) {
        index <- which((as.character(toplot$hybrids) == paste("hyb", as.character(hybrid[j]), sep = "")) & (!is.na(toplot$pvalues)) & (toplot$copulae == cops[i]))
        cur <- toplot[index, ]
        if (nrow(cur) == 1) {
          toplot <- toplot[-which(as.character(toplot$hybrids) == paste("hyb", as.character(hybrid[j]), sep = "")), ]
          cnt <- c(cnt, j)
        }
      }
    }
    cnt <- sort(unique(cnt))
    if (length(cnt) == 1) {
      warning(paste(paste("'hyb", as.character(hybrid[cnt]), "'", sep = ""), "is removed from plotting as the chosen 'inf.method' is not applicable."))
    }
    if (length(cnt) > 1) {
      warning(paste("The hybrids", paste(as.character(hybrid[cnt]), collapse = ","), "are removed from plotting as the chosen 'inf.method' is not applicable."))
    }
    if (length(cnt) > 0) {
      hybrid <- hybrid[-cnt]
    }
  }

  # arguments given in the ellipsis
  args <- list(...)

  # base arguments which are not to be touched by the user
  tabu_args <- c("formula", "data", "sortx", "ylim", "xaxt", "xlim", "ylab")
  if (any(is.element(names(args), tabu_args))) {
    stop("Please do not alter the arguments 'formula', 'data', 'sortx', 'ylim', 'xaxt', 'xlim' or 'ylab'.")
  }

  if ((length(x) > 2) & (length(hybrid) > 2)) {
    adj <- 0.75
  } else {
    adj <- 0.5
  }

  suppressWarnings(pirateplot(
    formula = pvalues ~ copulae + hybrids, data = toplot, theme = theme, sortx = "sequentiell",
    ylim = c(0, 1), xaxt = "n", xlim = c(adj, ((length(x) + 1) * (length(hybrid)) - adj)),
    ylab = "p-value", point.bg = point.bg, point.col = point.col, point.cex = point.cex,
    jitter.val = jitter.val, pal = pal, bean.b.o = bean.b.o, inf.method = inf.method, ...
  ))

  # vertical lines and axis
  if (length(hybrid) > 1) {
    abline(v = seq(length(x) + 1, (length(x) + 1) * (length(hybrid) - 1), (length(x) + 1)), lty = 3)
  }
  axis(side = 1, at = 1:((length(x) + 1) * length(hybrid) - 1), labels = c(rep(c(names(x), NA), (length(hybrid) - 1)), names(x)), las = 2)

  # hyb texts
  ats <- cumsum(c((length(x) + 1) / 2, rep(length(x) + 1, length(hybrid) - 1)))
  for (i in 1:length(hybrid)) {
    mtext(text = desired_hybs[i], side = 3, at = ats[i])
  }
}
