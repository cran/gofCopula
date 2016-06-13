gofADChisq = function(...) {
  warning("The test gofADChisq was renamed to gofRosenblattChisq. Please use the new phrase.")
  gofRosenblattChisq(...)
}

gofADGamma = function(...) {
  stop("The test gofADGamma was renamed to gofRosenblattGamma. Please use the new phrase.")
  gofRosenblattGamma(...)
}