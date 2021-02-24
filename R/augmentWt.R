#' \code{augmentWt} Fit state-space model to weighted drift dive record
#' @param data Data frame with dive data output from the \code{pDrift} function.
#' @param binsize
#' @param windowsize
#' @param yLims
#' @param plotit
#' @details Takes a dive set of dive data from SRDLs (Sea Mammal Research Unit), with added variables from
#' \code{rbs} and \code{pDrift} and attempts to augment data to improve the state-space model fit by the \code{fitDrate}.
#' @return Same data frame as the input, but with added variables:
#' \item{\code{augWeight}}{Weights based on only the data augmentation}
#' \item{\code{combWeight}}{Weights based on the product of basic weights (from the \code{pDrift} function) and the augmented weights}
#' @family Drift dive functions
#' @seealso \code{\link{rbs}} Reverse Broken Stick algorithm
#' @seealso \code{\link{check.pDrift}} to interactively check drift dive weightings
#' @seealso \code{\link{pDrift}} to calculate drift rate and probability weightings for drift dive detection
#' @seealso \code{\link{fitPlot}} To plot output from \code{fitDrift}
#' @author Martin Biuw
#' @examples
#' rbs(data=dive, num=100, n.bs=4)
#' dive <- rbs(data=dive, num=NA, n.bs=4)
#' dive <- pDrift(dive)
#' dive <- augmentWt(dive)
#' fit <- fitDrate(dive, prefilter=0.01)
#' fit2 <- fitDrate(dive, weight='rRes*span*dWt', prefilter=0.01)
#' @export

augmentWt <- function(data=sdive %>% filter(ref==levels(sdive$ref)[i]),
                      binsize=0.02, windowsize=3, yLims=c(-0.5, 0.5), plotit=F) {

  if(plotit) {
    par(mfrow=c(2,1), mar=c(1,4,1,1))
    plot(drate~DE.DATE, data=data %>% filter(weight>0.1),
         ylim=yLims, cex=weight*2, pch=19, col=rgb(0,0,0,0.2))
  }

  mid.pt <- min(data$L.DAY, na.rm=T)
  new.wts <- data$weight
  while(mid.pt<max(data$L.DAY, na.rm=T)) {
    neigh <- mid.pt + (c(-1, 1)*(windowsize*86400))
    neigh[2] <- neigh[2]+86400
    neigh <- which(between(data$L.DAY, neigh[1], neigh[2]))
    if(length(neigh)>0) {
      ndr <- data$drate[neigh]
      nwt <- data$weight[neigh]
      #hdr <- hist(ndr, seq(-0.5, 0.5, length=30))
      bins <- seq(min(ndr, na.rm=T), max(ndr, na.rm=T), by=binsize)
      bins <- c(bins, tail(bins, 1)+binsize)
      hdr <- hist(ndr, bins, plot=F)
      haug <- hdr
      drbin <- cut(ndr, hdr$breaks)
      binwt <- aggregate(nwt, list(drbin), median, na.rm=T, drop=F)
      binwt$x[which(is.nan(binwt$x))] <- 0
      binwt$x[which(is.na(binwt$x))] <- 0

      haug$density <- haug$density*binwt$x
      haug$density <- haug$density/max(haug$density, na.rm=T)

      drbin <- as.numeric(cut(ndr, haug$breaks))
      new.wts[neigh] <- new.wts[neigh] * haug$density[drbin]
    }
    mid.pt <- mid.pt + (windowsize*86400)
  }

  new.wts <- new.wts/max(new.wts, na.rm=T)

  if(plotit) {
    plot(drate~DE.DATE, data=data,
       ylim=yLims, cex=2*((data$weight+new.wts)/2), pch=19, col=rgb(0,0,0,0.2))
  }

  data$augWeight <- new.wts
  data$combWeight <- (data$weight+new.wts)/2
  data
}



