#' \code{fitDrate} Fit state-space model to weighted drift dive record
#' @param data Data frame with dive data output from the \code{pDrift} function.
#' @param weight Select weightingbvariables to be used. Default is to use \code{weight} created by \code{pDrift}.
#' See example for how to specify other weighting
#' @param mixture Select which mixture distribution to use ('nt' for normal and t distribution, 'nn' for mixture of normals)
#' @param prefilter Should any prefiltering be done (based on selected weight) prior to running the model.
#' @details Takes a dive set of dive data from SRDLs (Sea Mammal Research Unit), with added variables from
#' \code{rbs} and \code{pDrift} and fits a state-space model to change over time (using TMB).
#' @return A list with several objects from fit:
#' \item{\code{data}}{Original input data}
#' \item{\code{w}}{Weighting used for probability of drift dive}
#' \item{\code{tSeq}}{Regular time series (e.g. whole days) where predictions are obtained}
#' \item{\code{obj}}{TMB object}
#' \item{\code{opt}}{TMB optimization}
#' \item{\code{par}}{Fitted parameters}
#' \item{\code{rep}}{TMB report}
#' \item{\code{allsd}}{Standard deviation of fitted values}
#' \item{\code{plsd}}{Fitted values}
#' \item{\code{pred}}{Data frame with regular time sequence, fitted falues, sd of the fit and lower and upper 95% confidence intervals}
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
#' fit <- fitDrate(dive, prefilter=0.01)
#' fit2 <- fitDrate(dive, weight='rRes*span*dWt', prefilter=0.01)
#' @export

fitDrate <- function(data, weight='default', mixture=1, prefilter=0) {
  require(TMB)
  ##compile("./src/drifteR.cpp")
  ##dyn.load(dynlib("drifteR"))

  if(weight[1]=='default') {
    data$wt <- data$weight
  } else {
    if(length(weight)==1) {
      data$wt <- data[,match(weight, names(data))]
    } else {
      nmatch <- unlist(lapply(weight, function(x) grep(x, names(data))))
      wmat <- data[,nmatch]
      w <- apply(wmat, 1, prod)
      w <- w/max(w, na.rm=T)
      data$wt <- w
    }
  }
  data <- data[which(!is.na(data$DE.DATE) & !is.na(data$wt)),]
  data <- data[order(data$DE.DATE),]
  data <- data[which(is.finite(data$drate)),]


  if(!is.na(prefilter)) data <- data[which(data$wt>prefilter),]

  ## Time regularization:
  ## Create vector of whole sequential days:

  Date <- as.POSIXct(strptime(format(data$DE.DATE, '%Y-%m-%d'), '%Y-%m-%d'), tz='GMT')
  tSeq <- seq(min(Date)-43200, max(Date)+(86400+43200), by='days')

  ## Find which whole day a particular observation falls within:
  tIDX <- unlist(lapply(data$DE.DATE, function(x) which(tSeq>=x)[1]))

  ## Time difference between time of observation and the next whole day:
  j <- as.numeric(difftime(tSeq[tIDX], data$DE.DATE, units='secs'))/86400

  ##w <- data$wt/max(data$wt)
  w <- data$wt
  dat <- list(slope=data$drate, idx=tIDX, j=j, w=w)

  if(dat$idx[1]>1) dat$idx <- dat$idx-(dat$idx[1]-1)

  newdat <- list(slope=dat$slope[1], idx=dat$idx[1], j=dat$j[1], w=dat$w[1])
  for(i in c(2:length(dat$slope))) {
    if(dat$idx[i]> newdat$idx[length(newdat$idx)]+1) {
      tdiff <- dat$idx[i]-newdat$idx[length(newdat$idx)]
      newdat$slope <- c(newdat$slope, rep(NA, tdiff))
      newdat$idx <- c(newdat$idx, newdat$idx[length(newdat$idx)]+c(1:tdiff))
      newdat$j <- c(newdat$j, rep(NA, tdiff))
      newdat$w <- c(newdat$w, rep(NA, tdiff))
    } else {
      newdat$slope <- c(newdat$slope, dat$slope[i])
      newdat$idx <- c(newdat$idx, dat$idx[i])
      newdat$j <- c(newdat$j, dat$j[i])
      newdat$w <- c(newdat$w, dat$w[i])
    }
  }

  dat <- newdat
  dat$mx <- ifelse(mixture=='nl', 1, 2)

  ##data$u0 <- data$slope[head(which(data$w>0.3), 1)]
  ##data$w <- rep(1, length(data$w))
  rm(newdat)

  parameters <- list(
    logSdProc=0,
    logSdObs=0,
    u0=0,
    u=rep(0,length(tSeq)),
    logitp=0
  )

  obj <- MakeADFun(dat,parameters,random="u",DLL="drifteR")

  opt<-nlminb(obj$par,obj$fn,obj$gr)

  pl <- try(obj$env$parList(), silent=T)
  rep<-try(sdreport(obj, getJointPrecision=TRUE), silent=T)
  allsd<-try(sqrt(diag(solve(rep$jointPrecision))), silent=T)
  if(class(allsd) != 'try-error') {
    plsd <- obj$env$parList(par=allsd)
  } else {
    plsd <- allsd
  }
  pred <- data.frame(DE.DATE=tSeq, fit=pl$u, sd.fit=plsd$u)
  pred$sd.fit[which(is.nan(pred$sd.fit))] <- NA
  pred$lwr <- pred$fit-(2*pred$sd.fit)
  pred$upr <- pred$fit+(2*pred$sd.fit)

  if(any(is.na(pred$sd.fit))) {
    pred$segment <- rep(NA, nrow(pred))
    seg <- 1
    if(!is.na(pred$sd.fit[1])) pred$segment[1] <- seg
    for(i in 2:nrow(pred)) {
      if(!is.na(pred$sd.fit[i])) {
        if(is.na(pred$sd.fit[i-1])) seg <- seg+1
        pred$segment[i] <- seg
      }
    }
    if(min(pred$segment, na.rm=T)>1) pred$segment <- pred$segment - min(pred$segment, na.rm=T) + 1
  } else {
    pred$segment <- rep(1, nrow(pred))
  }

  list(data=data, w=data$wt, tSeq=tSeq, obj=obj, opt=opt, par=pl,
       rep=rep, allsd=allsd, plsd=plsd, pred=pred)
}

#' \code{fitPlot} Plot output from \code{fitDrift}
#' @param fit Output list from the \code{fitDrift} function.
#' @param yLims Set range of y axis
#' @param fitCol Colour to be used for estimated drift trajectory
#' @param fitAlpa Transparency to e used for estimated drift trajectory
#' @param fitLwd Line width to be used for estimated mean drift trajectory
#' @param ptMag Scaling for point size
#' @param ptAlpha Scaling for point transparency (usually has to be 1 or less)
#' @return A plot with data and the fitted estimates with confidence intervals
#' @family Drift dive functions
#' @seealso \code{\link{rbs}} Reverse Broken Stick algorithm
#' @seealso \code{\link{check.pDrift}} to interactively check drift dive weightings
#' @seealso \code{\link{pDrift}} to calculate drift rate and probability weightings for drift dive detection
#' @seealso \code{\link{fitDrate}} To fit state-space model to weighted data to estimate drift rate change over time
#' @author Martin Biuw
#' @examples
#' rbs(data=dive, num=100, n.bs=4)
#' dive <- rbs(data=dive, num=NA, n.bs=4)
#' dive <- pDrift(dive)
#' fit <- fitDrate(dive, prefilter=0.01)
#' fit2 <- fitDrate(dive, weight='rRes*span*dWt', prefilter=0.01)
#' fitPlot(fit)
#' fitPlot(fit2, add=T, fitCol=3)
#' @export


fitPlot <- function(fit, yLims=c(-0.5, 0.2),
                    fitCol=2, fitAlpha=0.2, add=F, fitLwd=1,
                    ptMag=2, ptAlpha=1) {
  if(add==F) {
    plot(drate~DE.DATE, data=fit$data, pch=19, col=rgb(0,0,0,ptAlpha*fit$w),
         cex=ptMag*fit$w, ylim=yLims)
  }
  confCol <- as.vector(col2rgb(fitCol))/255
  confCol <- rgb(confCol[1], confCol[2], confCol[3], fitAlpha)

  for(i in unique(fit$pred$segment)) {
    tmp <- fit$pred[which(fit$pred$segment==i),]
    if(nrow(tmp)==1) {
      segments(tmp$DE.DATE, tmp$lwr, tmp$DE.DATE, tmp$upr, col=confCol)
      lines(fit~DE.DATE, data=tmp, col=fitCol, lwd=fitLwd)
    } else {
      polygon(c(tmp$DE.DATE, rev(tmp$DE.DATE)), c(tmp$lwr, rev(tmp$upr)), col=confCol, border=F)
      lines(fit~DE.DATE, data=tmp, col=fitCol, lwd=fitLwd)
    }
  }
  lines(fit~DE.DATE, data=fit$pred, lty=2, col=fitCol, lwd=1)
}


#' \code{plotGG} Create ggplot output from \code{fitDrift}
#' @param fit Either a dataframe of dives, including variables from \code{rbs} and \code{pDrift}, or
#' output list from the \code{fitDrift} function.
#' @param yLims Set range of y axis
#' @param fitCol Colour to be used for estimated drift trajectory
#' @param fitAlpa Transparency to e used for estimated drift trajectory
#' @param fitLwd Line width to be used for estimated mean drift trajectory
#' @param ptSize Vector of size range of points.
#' @param to.plotly Convert to interactive plot using \code{plotly::ggplotly}
#' @return A plot with data and (if used on fitted model list) the fitted estimates with confidence intervals
#' @family Drift dive functions
#' @seealso \code{\link{rbs}} Reverse Broken Stick algorithm
#' @seealso \code{\link{check.pDrift}} to interactively check drift dive weightings
#' @seealso \code{\link{pDrift}} to calculate drift rate and probability weightings for drift dive detection
#' @seealso \code{\link{fitDrate}} To fit state-space model to weighted data to estimate drift rate change over time
#' @author Martin Biuw
#' @examples
#' rbs(data=dive, num=100, n.bs=4)
#' dive <- rbs(data=dive, num=NA, n.bs=4)
#' dive <- pDrift(dive)
#' fit <- fitDrate(dive, prefilter=0.01)
#' fit2 <- fitDrate(dive, weight='rRes*span*dWt', prefilter=0.01)
#' plotGG(dive)
#' potGG(fit)
#' @export

plotGG <- function(fit, yLims=c(-0.5, 0.2), weight='weight',
                    fitCol=2, fitAlpha=0.2, fitLwd=0.5,
                    ptSize=c(0.5, 3), to.plotly=F) {
  require(ggplot2)

  if(class(fit)=='list') {
    dat <- fit$data
    dat$w <- fit$w
  } else {
    dat <- fit
    dat$w <- dat[,match(weight, names(dat))]
  }

  dp <- ggplot(dat, aes(x=DE.DATE, y=drate, size=w, alpha=w)) +
    geom_point() +
    scale_size(range = ptSize) +
    ylim(yLims) +
    theme_bw() +
    xlab('') +
    ylab(expression(paste('Drift rate (m', s^-1, ')')))

  if(class(fit)=='list') {
    if(to.plotly) fitAlpha <- min(c(1, 2*fitAlpha))

    confCol <- as.vector(col2rgb(fitCol))/255
    confCol <- rgb(confCol[1], confCol[2], confCol[3], fitAlpha)

    for(i in unique(fit$pred$segment)) {
      tmp <- fit$pred[which(fit$pred$segment==i),]
      dp <- dp + geom_line(data=tmp, aes(x=DE.DATE, y=fit),
                           colour=fitCol, size=fitLwd, show.legend=F, inherit.aes=F) +
                geom_ribbon(data=tmp, aes(x=DE.DATE, y=fit, ymin=lwr, ymax=upr, alpha=fitAlpha),
                           fill=confCol, show.legend=F, inherit.aes=F)
    }
  }

  if(to.plotly) {
    require(plotly)
    dp$labels$y <- 'Drift rate (m/s)'
    dpl <- ggplotly(dp)
    dpl
  } else {
    dp
  }
}

