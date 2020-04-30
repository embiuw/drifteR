tfClassifier <- function(data) {
  require(keras)
  model <- load_model_hdf5('diveClassTF.h5')
  system.file("extdata", "diveClassTF.h5", package = "drifteR", mustWork = TRUE)
  if(dat$N.DEPTHS==5) dat <- strip.last(data)

  dat <- standardise.dive(dat)
  attr(dat, 'dimnames') <- NULL

  predictions <- model %>% predict(dat)

  predictions <- data.frame(predictions)
  names(predictions) <- paste0('tf', c("SQ","U","V","W","DR","R"))

  cbind(data, predictions)
}
