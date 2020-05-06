#' @export

standardise.dive <- function(dat, features=c(paste0('T', c(1:4)), paste0('D', c(1:4))),
                             center=T) {
  dat <- cbind(T0=rep(0, nrow(dat)), dat[,match(features[c(1:4)], names(dat))],
               T5=rep(100, nrow(dat)), D0=rep(0, nrow(dat)),
               dat[,match(features[c(5:8)], names(dat))], D5=rep(0, nrow(dat)))
  dat <- t(apply(dat, 1, function(x) {
    c((x[1:6]/100), x[7:12])
  }))
  dat <- t(apply(dat, 1, function(x) {
    max.d <- max(x[6:10])
    c(x[1:6], (x[7:12]/max.d))
  }))
  if(center) {
    dat-0.5
  } else {
    dat
  }
}


