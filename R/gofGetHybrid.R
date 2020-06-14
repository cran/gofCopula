gofGetHybrid <- function(result, p_values = NULL, nsets = NULL) {
  if (!inherits(result, "gofCOP")) {
    stop("Please input an object of class 'gofCOP'. Such an object will be returned by functions of this package. If you input an object obtained from 'gof()', then input the result for all copula.")
  }
  if (length(nsets) > 1) {
    stop("'nsets' has to be a single integer entry.")
  }

  res_list <- list()
  for (j in seq_along(result)) {
    ## checking setup
    if (length(result[[j]]) == 7) {
      tmp <- rownames(result[[j]]$res.tests)
      index <- which(startsWith(tmp, "hybrid"))
      if (length(index) > 0) {
        res_length <- length(tmp[-index])
      } else {
        res_length <- length(tmp)
      }
    } else {
      res_length <- 1
    }

    num_tests <- length(p_values) + res_length
    if (num_tests <= 1) {
      stop("The input should contain information of at least two different tests.")
    }
    if (!is.null(nsets)) {
      if (nsets < 1 | nsets > num_tests) {
        stop("Please set nsets larger or equal than 1 and smaller or equal than the number of single tests. Otherwise hybrid testing is not meaningful.")
      }
    }

    ## getting combinations
    which_comb <- list()
    for (i in seq_len(2^num_tests)) {
      which_comb[[i]] <- which(as.integer(intToBits(i)) == 1)
    }
    comb_exist <- which_comb[which(unlist(lapply(which_comb, length)) > 1)]

    ## building results
    # names and p-values of single tests
    s_res_names <- rownames(result[[j]]$res.tests)
    s_res_p <- result[[j]]$res.tests[, 1]
    index <- which(startsWith(s_res_names, "hybrid"))
    if (length(index) > 0) {
      s_res_names <- s_res_names[-index]
      s_res_p <- s_res_p[-index]
    }
    s_p_names <- if (!is.null(p_values)) {
      if (is.null(names(p_values))) {
        paste0("Test_", LETTERS[1:length(p_values)])
      } else {
        names(p_values)
      }
    } else {
      NULL
    }
    s_names <- c(s_res_names, s_p_names)
    p <- c(s_res_p, p_values)

    # combinations
    if (!is.null(nsets)) {
      if (nsets > 1) {
        index <- which(lapply(comb_exist, FUN = function(x) {
          length(x) == nsets
        }) == TRUE)
        comb_wanted <- comb_exist[index]
        for (i in seq_along(comb_wanted)) {
          p <- c(p, min(nsets * min(p[comb_wanted[[i]]]), 1))
        }
        comb_names <- paste("hybrid(", lapply(comb_wanted, paste, collapse = ", "), ")", sep = "")
        res <- matrix(p, ncol = 1, dimnames = list(c(s_names, comb_names), "p.value"))
      } else {
        res <- matrix(p, ncol = 1, dimnames = list(c(s_names), "p.value"))
      }
    } else {
      for (i in seq_along(comb_exist)) {
        p <- c(p, min(length(comb_exist[[i]]) * min(p[comb_exist[[i]]]), 1))
      }
      comb_names <- paste("hybrid(", lapply(comb_exist, paste, collapse = ", "), ")", sep = "")
      res <- matrix(p, ncol = 1, dimnames = list(c(s_names, comb_names), "p.value"))
    }
    res_list[[j]] <- list(
      method = "Hybrid test p-values for given single tests.",
      copula = result[[j]]$copula,
      margins = result[[j]]$margins,
      param.margins = result[[j]]$param.margins,
      theta = result[[j]]$theta,
      df = result[[j]]$df,
      res.tests = res
    )
  }
  names(res_list) <- names(result)

  # output
  return(structure(
    class = "gofCOP",
    res_list
  ))
}
