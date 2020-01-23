gofOutputHybrid = function(result, tests = NULL, nsets = NULL) {
  
  if (class(result) != "gofCOP") {stop("Please input an object of class 'gofCOP'. Such an object will be returned by functions of this package. If you input an object obtained from 'gof()', then input the result for one copula only.")}
  if (!any(grepl("hybrid", rownames(result[[1]]$res.tests)))) {stop("Please input an object containing hybrid test results.")}
  numb_tests = sum(!grepl("hybrid", rownames(result[[1]]$res.tests)))
  if (isTRUE(nsets > numb_tests)) {stop("The number of tests to obtain the hybrid test for cannot be larger than the number of tests available.")}
  if (!all(is.element(tests, 1:numb_tests))) {stop("At least one of the 'tests' is not stored in the object.")}
  
  res_list = list()
  for (j in 1:length(result)) {

    which_comb =list()
    for(i in 1:(2^numb_tests)){
      which_comb[[i]]=which(as.integer(intToBits(i)) == 1)
    }
    comb_exist = which_comb[which(unlist(lapply(which_comb, length)) > 1)]
    
    which_tests = rep(TRUE, length(comb_exist))
    if (!is.null(tests)) {
      which_tests = sapply(comb_exist, function(x) any(is.element(x,tests)))
    }
    which_nsets = rep(TRUE, length(comb_exist))
    if (!is.null(nsets)) {
      which_nsets = sapply(comb_exist, function(x) (length(x) == nsets))
    }
    which_tests_nsets = which_nsets & which_tests
    
    hybrid.results = result[[j]]$res.tests[-c(1:numb_tests),]
    new.res.tests = result[[j]]$res.tests[1:numb_tests,]
    names_new.res.tests = rownames(new.res.tests)
    if (length(which_tests_nsets) == 1) {
      if (which_tests_nsets == TRUE) {
        new.res.tests = rbind(new.res.tests, hybrid.results)
        rownames(new.res.tests) = c(names_new.res.tests, rownames(result[[j]]$res.tests)[-c(1:numb_tests)])
      }
    } else {
      new.res.tests = rbind(new.res.tests, hybrid.results[which_tests_nsets,])
      rownames(new.res.tests) = c(names_new.res.tests, rownames(hybrid.results)[which_tests_nsets])
    }
    res_list[[j]] = list(method = result[[j]]$method,
                         copula = result[[j]]$copula, 
                         margins = result[[j]]$margins, 
                         param.margins = result[[j]]$param.margins, 
                         theta = result[[j]]$theta,
                         df = result[[j]]$df,
                         res.tests = new.res.tests)
  }
  names(res_list) = names(result)
  
  structure(class = "gofCOP", 
            res_list)
}
