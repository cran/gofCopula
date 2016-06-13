.ErrorHandler = function(x) {
  if (x$message == "param can be negative only for dim = 2"){
    stop("The dependence parameter is negative for the dataset. For this copula can this be only the case if the dimension is 2. Therefore is this not an appropriate copula for the dataset. Please consider to use another one.")
  } else {
    stop(x$message)
  }
}