.get.time <- function(x) {
  x.day <- floor(x / 86400)
  x.remainder <- x %% 86400
  x.hour <- floor(x.remainder / 3600)
  x.remainder <- x.remainder %% 3600
  x.min <- floor(x.remainder / 60)
  x.remainder <- x.remainder %% 60
  out <- list(x.day, x.hour, x.min, x.remainder)
  class(out) <- "goftime"
  if (any(!is.na(unlist(out)))) {
    return(out)
  }
}

print.goftime <- function(x, ...) {
  print(sprintf("The computation will take approximately %d d, %d h, %d min and %d sec.", x[[1]], x[[2]], x[[3]], x[[4]]))
}
