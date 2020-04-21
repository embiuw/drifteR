#' \code{rbs} Reverse Broken Stick algorithm
#' @param data Data frame with dive data output from Satellite Relay Data Loggerds (SRDLs)
#' @param num Which dive in data should be processed. If NA, the entire dataset is processed.
#' @param n.bs Number of valid broken stick time/depth points in data.
#' @details Takes a standard dive set of dive data from SRDLs (Sea Mammal Research Unit),
#' and returns three new variables:
#' \code{order}: Character, order in which inflction points were extracted
#' \code{max.res}: Numeric, maximum remaining residual (in meters) unaccounted for by summarized dive profile
#' \code{reordered}: Logical, if dive was reorganized from inconsistent ordering in original file
#' NOTE! Names in the data frame sould have underscores replaced by dots prior to running this function!
#' @family Drift dive functions
#' @seealso \code{\link{pDrift}} to calculate drift rate and probability weightings for drift dive detection
#' @seealso \code{\link{check.pDrift}} to interactively check drift dive weightings
#' @seealso \code{\link{fitDrate}} To fit state-space model to weighted data to estimate drift rate change over time
#' @seealso \code{\link{fitPlot}} To plot output from \code{fitDrate}
#' @author Martin Biuw
#' @examples
#' rbs(data=dive, num=100, n.bs=4)
#' dive <- rbs(data=dive, num=NA, n.bs=4)
#'
#' @export

rbs <- function(data=dive, num=100, n.bs=4) {
  rbfun <- function(plotit=T) {
    d <- as.numeric(data[num, match(paste0('D', c(1:n.bs)), names(data))])
    t <- as.numeric(data[num,match(paste0('T', c(1:n.bs)), names(data))])

    ## Check time ordering:
    if(any(order(t)!=c(1:n.bs))) {
      tOrd <- order(t)
      t <- t[tOrd]
      d <- d[tOrd]
    }
    d <- c(0, d, 0)

    if(data$DIVE.DUR[num]>0 & all(t[-1]>0) &
       all(t[-length(t)]<100) & all(d[-c(1, length(d))]>0)) {
      t <- c(0, t/100, 1)*data$DIVE.DUR[num]
      bs <- which(d==max(d))
      if(length(bs)>1) {
        if(length(bs)==n.bs) {
          bs <- which.min(abs(t-(0.5*tail(t, 1))))
        }
      }

      while(length(bs)<(n.bs-1)) {
        dApp <- c(0, d[sort(bs)], 0)
        tApp <- c(0, t[sort(bs)], tail(t, 1))
        app <- approx(tApp, dApp, t)
        bs <- c(bs, which.max(abs(d-app$y)))
        if(tail(bs, 1)==1) {
          bs <- bs[-length(bs)]
          bs <- c(bs, setdiff(c(1:n.bs)+1, bs))
        }

#        if(plotit) {
#          plot(t, -d, type='o', pch=21, bg='grey')
#          points(tApp, -dApp, type='o', pch=21, bg=2, lty=2)
#          segments(t[tail(bs,1)], -d[tail(bs,1)], t[tail(bs,1)], -app$y[tail(bs,1)], lty=3)
#        }
      }

      dApp <- c(0, d[sort(bs)], 0)
      tApp <- c(0, t[sort(bs)], tail(t, 1))
      app <- approx(tApp, dApp, t)
      bs <- c(bs, setdiff(c(2:5), bs))

      if(plotit) {
        plot(t, -d, type='o', pch=21, bg='grey')
        points(t[-tail(bs, 1)], -d[-tail(bs, 1)], type='o', pch=21, bg=2, lty=2)
        segments(t[tail(bs,1)], -d[tail(bs,1)], t[tail(bs,1)], -app$y[tail(bs,1)], lty=3)
        points(t, -d, pch=21, bg=2)
        for(i in 1:length(bs)) {
          text(t[bs[i]], -d[bs[i]], i, pos=3)
        }
      }

      reordered <- ifelse(exists('tOrd'), tOrd, NA)
      list(order=paste(bs-1, collapse='.'), max.res=max(abs(d-app$y)), reordered=reordered)
    } else {
      list(order=NA, max.res=NA, reordered=NA)
    }
  }

  if(!is.na(num)) {
    rbfun()
  } else {
    data$order <- rep(NA, nrow(data))
    data$max.res <- rep(NA, nrow(data))
    data$reordered <- rep(NA, nrow(data))

    cat('Calculating reverse broken stick\n')

    for(i in 1:nrow(data)) {
      pd <- round(50*(i/nrow(data)))
      pu <- 50-pd
      catstr <- paste0('##', paste0(rep('-', pd), collapse=''),
                       paste0(rep(' ', pu), collapse=''), '## ', 2*pd, '% done      \r')
      cat(catstr)
      flush.console()
      num <- i
      tmp <- try(rbfun(plotit=F), silent=T)
      if(class(tmp)!='try-error') {
        data$order[i] <- tmp$order
        data$max.res[i] <- tmp$max.res
        data$reordered[i] <- tmp$reordered[i]
      }
    }
    cat('\n\nDone!\n\n')
    flush.console()
    data
  }
}




