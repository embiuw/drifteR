#' @export

tfClassifier <- function(data) {
  require(keras)

  mod.nam <- system.file("extdata", "diveClassTF.h5", package = "drifteR", mustWork = TRUE)
  model <- load_model_hdf5(mod.nam)

  if(data$N.DEPTHS[1]==5) {
    dat <- strip.last(data)
  } else {
    dat <- data
  }

  dat <- standardise.dive(dat)
  attr(dat, 'dimnames') <- NULL

  predictions <- model %>% predict(dat)

  predictions <- data.frame(predictions)
  names(predictions) <- paste0('tf', c("SQ","U","V","W","DR","R"))

  cbind(data, predictions)
}
