gofWhich = function(copula, d){
  if (d <= 1){stop("The dimension d must be 2 or even higher.")}
  if(is.element(copula, c("independence", "normal", "gaussian", "t", "clayton", "frank", "gumbel", "joe", 
                          "amh", "survival clayton", "survival gumbel", "survival joe", "90 clayton", 
                          "90 gumbel", "90 joe", "270 clayton", "270 gumbel", "270 joe", 
                          "bb1", "bb6", "bb7", "bb8", "survival bb1", "survival bb6", 
                          "survival bb7", "survival bb8", "90 bb1", "90 bb6", "90 bb7", 
                          "90 bb8", "270 bb1", "270 bb6", "270 bb7", "270 bb8")) == F){stop("This copula is for no test implemented.")}
  
  which1 = which2 = which3 = c()
  if (d >= 2 & is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel"))){
    which1 = c("gofRosenblattSnB", "gofRosenblattSnC", "gofRosenblattChisq", "gofRosenblattGamma", "gofSn", "gofKendallCvM", "gofKendallKS")
  }
  if (d == 2 & is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel"))){
    which2 = c("gofKernel", "gofWhite")
  }
  if (d == 2 || d == 3 & is.element(copula, c("normal", "gaussian", "t", "clayton", "frank", "gumbel"))){
    which3= c("gofPIOSRn", "gofPIOSTn")
  }
  res = c("gofHybrid", which1, which2, which3)
  res
}

gofWhichCopula = function(test){
  if (is.element(test, c("gofRosenblattSnB", "gofRosenblattSnC", "gofRosenblattChisq", "gofRosenblattGamma", "gofSn", "gofRn", "gofKendallCvM", "gofKendallKS", "gofPIOSRn", "gofPIOSTn", "gofKernel", "gofWhite"))){
    res = c("normal", "t", "clayton", "frank", "gumbel")#, "amh", "joe")
  } else if (is.element(test, c("gofHybrid"))){
    print("The available copula depend on the used tests.")
  } else {
    print("The test is not implemented.")
  }
  res
}