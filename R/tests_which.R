gofTest4Copula = function(copula = NULL, d = 2){
  if (length(copula) > 1) {stop("The argument 'copula' has to be of length 1.")}
  if (is.null(copula)) {copula = eval(formals(gofRosenblattSnB)$copula)}
  if (!is.numeric(d)) {stop("d must be a numeric.")}
  if (d <= 1){stop("The dimension d must be 2 or even higher.")}
  cops = lapply(ls(pos = "package:gofCopula"), function(x) {try(eval(formals(x)$copula), silent = TRUE)})
  cops_pos = sapply(cops, function(x) {any(is.element(x, copula))})
  if (!any(cops_pos)) {stop("This copula is for no test implemented.")}
  res = which(sapply(ls(pos = "package:gofCopula")[cops_pos], .gofCheckDim) >= d)
  c("gofHybrid", names(res))
}

gofCopula4Test = function(test){
  if (!is.character(test)) {stop("The argument 'test' has to be a character.")}
  if (length(test) > 1) {stop("The argument 'test' has to be of length 1.")}
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
    if (is.element(test, c("gofRosenblattSnB", "gofRosenblattSnC", "gofRosenblattChisq", "gofRosenblattGamma", "gofSn", "gofKendallCvM", "gofKendallKS", "gofCustomTest"))) {
        Inf
    } else if (is.element(test, c("gofKernel", "gofWhite"))) {
        2
    } else if (is.element(test, c("gofPIOSRn", "gofPIOSTn"))) {
        3
    } else {
        NULL
    }
}




gofWhich = function(copula = NULL, d = 2){
  warning("The function 'gofWhich' was renamed to 'gofTest4Copula'. Please use 'gofTest4Copula'.")
  gofTest4Copula(copula = copula, d = d)
}
gofWhichCopula = function(test){
  warning("The function 'gofWhichCopula' was renamed to 'gofCopula4Test'. Please use 'gofCopula4Test'.")
  gofCopula4Test(test = test)
}
