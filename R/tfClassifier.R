#' @export

tfClassifier <- function(data, n.bs=NA) {
  require(keras)

  mod.nam <- system.file("extdata", "diveClassTF.h5", package = "drifteR", mustWork = TRUE)
  model <- load_model_hdf5(mod.nam)

  if('N.DEPTHS' %in% names(data)) {
    if(data$N.DEPTHS[1]==5) {
      dat <- strip.last(data)
    } else {
      dat <- data
    }
  } else {
    if(is.na(n.bs)) {
      stop('no input provided to n.bbs. Please specify number of inflection points in data')
    } else {
      if(n.bs==5) {
        dat <- strip.last(data)
      } else {
        dat <- data
      }
    }
  }

  dat <- standardise.dive(dat)
  attr(dat, 'dimnames') <- NULL

  predictions <- model %>% predict(dat)

  predictions <- data.frame(predictions)
  names(predictions) <- paste0('tf', c("SQ","U","V","W","DR","R"))

  cbind(data, predictions)
}
