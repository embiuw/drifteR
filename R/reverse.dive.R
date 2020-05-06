#' @export

reverse.dive <- function(dat, centered=T) {
  t(apply(dat, 1, function(x) {
    c(rev(-1*x[1:6]), rev(x[7:12]))
  }))
}

