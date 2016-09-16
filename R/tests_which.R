gofWhich = function(copula, d){
    if (!is.numeric(d)) {stop("d must be a numeric.")}
  if (d <= 1){stop("The dimension d must be 2 or even higher.")}
  cops = lapply(ls(pos = "package:gofCopula"), function(x) {try(eval(formals(x)$copula), silent = TRUE)})
  cops_pos = sapply(cops, function(x) {any(is.element(x, copula))})
  if (!any(cops_pos)) {stop("This copula is for no test implemented.")}
  res = which(sapply(ls(pos = "package:gofCopula")[cops_pos], .gofCheckDim) >= d)
  c("gofHybrid", names(res))
}

gofWhichCopula = function(test){
    if (is.element(test, c("gofHybrid"))){
        print("The available copula depend on the used tests.")
    } else {
        res = try(eval(formals(eval(parse(text = test)))$copula), silent = TRUE)
        if (class(res) == "try-error") {
            print("The test is not implemented.")
        } else {
            res
        }
    }
}

.gofCheckDim = function(test) {
    if (is.element(test, c("gofRosenblattSnB", "gofRosenblattSnC", "gofRosenblattChisq", "gofRosenblattGamma", "gofSn", "gofKendallCvM", "gofKendallKS"))) {
        Inf
    } else if (is.element(test, c("gofKernel", "gofWhite"))) {
        2
    } else if (is.element(test, c("gofPIOSRn", "gofPIOSTn"))) {
        3
    } else {
        NULL
    }
}

