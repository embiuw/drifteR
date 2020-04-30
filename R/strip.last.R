strip.last <- function(dat=hp6) {
  ## Must have run 'rbs' first!
  new.dat <- matrix(NA, nrow=nrow(dat), ncol=8)
  pb <- txtProgressBar(title = "Reorganising remaining inflection point data", min = 0,
                       max = nrow(dat), style=3)
  
  for(i in 1:nrow(dat)) {
    setTxtProgressBar(pb, i, title=paste("Reorganising remaining inflection point data",
                                         round(i/nrow(dat)*100, 0),
                                         "% done"))
    if(!is.na(dat$order[i])) {
      last <- tail(unlist(strsplit(dat$order[i], '.', fixed=T)), 1)
      keep.names <- setdiff(orig.names, paste0(c('T', 'D'), last))
      new.dat[i,] <- as.numeric(dat[i,match(keep.names, names(dat))])
    }  
  }
  
  close(pb)
  new.dat <- as.data.frame(new.dat)
  names(new.dat) <- c(paste0('T', c(1:4)), paste0('D', c(1:4)))
  new.dat
}

