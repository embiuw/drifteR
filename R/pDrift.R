#' \code{pDrift} Drift rate and weight calculations
#' @param data Data frame with dive data output from the \code{rbs} function.
#' @param depth.wpar Set mean ad standard deviation of half gaussian for depth weighting
#' @param plotit Logical, plot results or not
#' @param keep.all Logical, keep all data or only data containing valid drift rates and weights
#' @details Takes a dive set of dive data from SRDLs (Sea Mammal Research Unit), and calculates
#' drift rate based on segment between first two extracted inflection points, along with
#' probability weights for each dive being a drift dive.
#' @return Adds several new variables to teh input data frame:
#' \item{\code{direction}}{Direction of putative drift segment ('up' or 'down')}
#' \item{\code{drate}}{Rate of vertical drift (drift rate)}
#' \item{\code{skew}}{Position of deepest point relative to midpoint along dive duration}
#' \item{\code{span}}{Proportion of dive taken up by putative drift segment}
#' \item{\code{kink}}{Absolute difference between descent (for negative drates) or ascent (for positive drates) rate and drift rate}
#' \item{\code{rRes}}{Standardised 'reverse residual', (see details)}
#' \item{\code{dWt}}{Depth weighting based on input settings}
#' \item{\code{drWt}}{Weighting of drift rate probability given likely range for extreme range of lipid content (5% to 55%)}
#' \item{\code{weight}}{Combined standardized weight using all above weighting variables}
#' @family Drift dive functions
#' @seealso \code{\link{rbs}} Reverse Broken Stick algorithm
#' @seealso \code{\link{check.pDrift}} to interactively check drift dive weightings
#' @seealso \code{\link{fitDrift}} To fit state-space model to weighted data to estimate drift rate change over time
#' @seealso \code{\link{fitPlot}} To plot output from \code{fitDrift}
#' @author Martin Biuw
#' @examples
#' rbs(data=dive, num=100, n.bs=4)
#' dive <- rbs(data=dive, num=NA, n.bs=4)
#' dive <- pDrift(dive)
#' @export

pDrift <- function(data=dive, depth.wpar=c(100, 30), plotit=F, keep.all=T) {

  ## Which dives have valid data:
  which.valid <- which(!is.na(data$order))

  sub <- data[which.valid,]

  ## Extract first and second extracted inflection points:
  ifp1 <- unlist(lapply(sub$order, function(x) unlist(strsplit(x, '.', fixed=T))[1]))
  ifp2 <- unlist(lapply(sub$order, function(x) unlist(strsplit(x, '.', fixed=T))[2]))

  ## Second subsetting to omit dives where any inflection points are zero (i.e. start of dive'):
  sub <- sub[which(ifp1!='0' & ifp2!='0'),]

  ## Extract first and second extracted inflection points again on new subset:
  ifp1 <- unlist(lapply(sub$order, function(x) unlist(strsplit(x, '.', fixed=T))[1]))
  ifp2 <- unlist(lapply(sub$order, function(x) unlist(strsplit(x, '.', fixed=T))[2]))

  d2 <- d1 <- pr2 <- pr1 <- rep(NA, nrow(sub))

  ## Calculate time proportions for the two ifp's:

  pr1[which(ifp1=='1')] <- sub$T1[which(ifp1=='1')]/100
  pr1[which(ifp1=='2')] <- sub$T2[which(ifp1=='2')]/100
  pr1[which(ifp1=='3')] <- sub$T3[which(ifp1=='3')]/100
  pr1[which(ifp1=='4')] <- sub$T4[which(ifp1=='4')]/100

  pr2[which(ifp2=='1')] <- sub$T1[which(ifp2=='1')]/100
  pr2[which(ifp2=='2')] <- sub$T2[which(ifp2=='2')]/100
  pr2[which(ifp2=='3')] <- sub$T3[which(ifp2=='3')]/100
  pr2[which(ifp2=='4')] <- sub$T4[which(ifp2=='4')]/100

  ## Extract corresponding depths for the two ifp's:
  d1[which(ifp1=='1')] <- sub$D1[which(ifp1=='1')]
  d1[which(ifp1=='2')] <- sub$D2[which(ifp1=='2')]
  d1[which(ifp1=='3')] <- sub$D3[which(ifp1=='3')]
  d1[which(ifp1=='4')] <- sub$D4[which(ifp1=='4')]

  d2[which(ifp2=='1')] <- sub$D1[which(ifp2=='1')]
  d2[which(ifp2=='2')] <- sub$D2[which(ifp2=='2')]
  d2[which(ifp2=='3')] <- sub$D3[which(ifp2=='3')]
  d2[which(ifp2=='4')] <- sub$D4[which(ifp2=='4')]

  ## Determine direction of putative drifts:
  direction <- ifelse(pr1>pr2, 'down', 'up')
  direction[which(d1==d2)] <- 'flat'

  ## Calculate drift rates:
  drate <- rep(NA, nrow(sub))
  drate[which(direction=='down')] <- ((d2-d1)/(abs(pr1-pr2)*sub$DIVE.DUR))[which(direction=='down')]
  drate[which(direction=='up')] <- ((d1-d2)/(abs(pr1-pr2)*sub$DIVE.DUR))[which(direction=='up')]
  drate[which(direction=='flat')] <- 0


  ## Calculate skew index
  ## (i.e. how far from midway through the
  ## dive duration does max dive occur):

  skew <- 2*abs(pr1-0.5)

  ## Calculate time proportion between first two ifp's
  span <- abs(pr1-pr2)

  ## Calculate absolute difference between descent/ascent rate and drift rate
  kink <- rep(NA, nrow(sub))
  kink[which(direction=='down')] <- (abs(drate)/abs((sub$D1/(0.01*sub$T1*sub$DIVE.DUR))))[which(direction=='down')]
##  kink[which(direction=='down' & ifp2=='1')] <- (abs(drate)/abs((sub$D1/(0.01*sub$T1*sub$DIVE.DUR))))[which(direction=='down' & ifp2=='1')]
##  kink[which(direction=='down' & ifp2=='2')] <- (abs(drate)/abs((sub$D2/(0.01*sub$T2*sub$DIVE.DUR))))[which(direction=='down' & ifp2=='2')]
##  kink[which(direction=='down' & ifp2=='3')] <- (abs(drate)/abs((sub$D3/(0.01*sub$T3*sub$DIVE.DUR))))[which(direction=='down' & ifp2=='3')]

  kink[which(direction=='up')] <- (abs(drate)/abs((sub$D4/(0.01*(100-sub$T4)*sub$DIVE.DUR))))[which(direction=='up')]
##  kink[which(direction=='up' & ifp2=='2')] <- (abs(drate)/abs((sub$D2/(0.01*(100-sub$T2)*sub$DIVE.DUR))))[which(direction=='up' & ifp2=='2')]
##  kink[which(direction=='up' & ifp2=='3')] <- (abs(drate)/abs((sub$D3/(0.01*(100-sub$T3)*sub$DIVE.DUR))))[which(direction=='up' & ifp2=='3')]
##  kink[which(direction=='up' & ifp2=='4')] <- (abs(drate)/abs((sub$D4/(0.01*(100-sub$T4)*sub$DIVE.DUR))))[which(direction=='up' & ifp2=='4')]
  kink[which(kink>1)] <- 0

  ## Create weighting function based on skew, span and (standardised)
  ## maximum remaining residual:

##  rRes <- sub$max.res/sub$MAX.DEP

  internal <- which(abs(as.numeric(ifp1)-as.numeric(ifp2))>1)

  rRes <- rep(NA, nrow(sub))
  rRes[-internal] <- sub$max.res[-internal]/(abs(d1-d2)[-internal])

  for(i in internal) {
    n.bs <- length(unlist(strsplit(sub$order[i], '.', fixed=T)))
    if1 <- as.numeric(substr(sub$order[i], 1, 1))
    if2 <- as.numeric(substr(sub$order[i], 3, 3))
    d <- as.numeric(sub[i, match(paste0('D', c(1:n.bs)), names(sub))])
    t <- as.numeric(sub[i,match(paste0('T', c(1:n.bs)), names(sub))])
    t <- 0.01*t*sub$DIVE.DUR[i]
    app <- try(approx(t[c(if1, if2)], d[c(if1, if2)], t[c(if1:if2)]), silent=T)
    if(class(app)!='try-error') {
      res <- max(abs(d[c(if1:if2)] - app$y))
      rRes[i] <- res/diff(range(d[c(if1, if2)]))
    }
  }

  rRes[which(!is.finite(rRes))] <- NA
  rRes <- (1-plogis(rRes))/max(1-plogis(rRes), na.rm=T)

  ## Create second weighting function based on depth of shallowest of the two ifp's
  ## This uses a half gaussian down to depth threshold, specified in depth.par[1],
  ## with a standard deviation specified in depth.par[2]

  dWt <- rep(0, length(sub))
  dWt[which(direction=='down' | direction=='flat')] <- dnorm(d1[which(direction=='down' | direction=='flat')],
                                        depth.wpar[1], depth.wpar[2])/dnorm(depth.wpar[1], depth.wpar[1], depth.wpar[2])
  dWt[which(direction=='up')] <- dnorm(d2[which(direction=='up')],
                                      depth.wpar[1], depth.wpar[2])/dnorm(depth.wpar[1], depth.wpar[1], depth.wpar[2])
  dWt[which(direction=='down' & d1>depth.wpar[1])] <- 1
  dWt[which(direction=='flat' & d1>depth.wpar[1])] <- 1
  dWt[which(direction=='up' & d2>depth.wpar[1])] <- 1


  ## Make drift rate probability vector based on
  ## theoretical calculations of drift rate as a function of
  ## maximum likely range in lipid content (and morphometrics).
  ## Corresponding to lipid content from 5% to 55%

  drSeq <- seq(min(drate[which(is.finite(drate))], na.rm=T),
               max(drate[which(is.finite(drate))], na.rm=T), by=0.01)
  drSeq <- c(drSeq, tail(drSeq, 1)+0.01)
  drWt <- 0.5*dnorm(drSeq, -0.5, 0.1)+
          0.5*dnorm(drSeq, 0.4, 0.1)
  drWt <- drWt/max(drWt)
  drWt[which(drSeq>-0.5 & drSeq<0.4)] <- 1

  drWt <- approx(drSeq, drWt, drate)$y

  ## Add derived variables to subset data frame:
  sub$direction <- direction
  sub$drate <- drate
  sub$skew <- skew
  sub$span <- span
  sub$kink <- kink
  sub$rRes <- rRes
  sub$dWt <- dWt
  sub$drWt <- drWt
  ## Which way should kink go?
  ##sub$weight <- sub$skew*sub$span*(1-sub$kink)*sub$rRes*sub$dWt*sub$drWt
  sub$weight <- sub$skew*sub$span*sub$kink*sub$rRes*sub$dWt*sub$drWt

  sub$weight <- sub$weight/max(sub$weight[which(is.finite(sub$weight))], na.rm=T)
  sub$weight[which(!is.finite(sub$weight))] <- NA
  ##sub$weight[which(is.na(sub$weight))] <- 0

  if(keep.all) {
    ## And add new variables for corresponding dives in inout data frame_
    nNames <- setdiff(names(sub), names(data))

    nData <- data
    which.valid <- match(row.names(sub), row.names(data))
    for(i in nNames) eval(parse(text=paste0('nData$', i, ' <- rep(NA, nrow(nData))')))
    for(i in nNames) eval(parse(text=paste0('nData$', i, '[which.valid] <- sub$', i)))
    nData
  } else {
    sub
  }
}

#' \code{check.pDrift} Visual check of drift rate and weight calculations
#' @param data Data frame with dive data output from the \code{pDrift} function.
#' @param ptMag Scaling of point sizes to increase or decrease emphasis
#' @details Opens an interactive plot of entire dive record, where dives can be individually
#' checked by selecting them. Opens a second plot with the individual dive profile.
#' @return Interactive plot
#' @family Drift dive functions
#' @seealso \code{\link{rbs}} Reverse Broken Stick algorithm
#' @seealso \code{\link{pDrift}} to calculate drift rate and probability weightings for drift dive detection
#' @seealso \code{\link{fitDrift}} To fit state-space model to weighted data to estimate drift rate change over time
#' @seealso \code{\link{fitPlot}} To plot output from \code{fitDrift}
#' @author Martin Biuw
#' @examples
#' rbs(data=dive, num=100, n.bs=4)
#' dive <- rbs(data=dive, num=NA, n.bs=4)
#' dive <- pDrift(dive)
#' check.pDrift(dive)
#' @export

check.pDrift <- function(data, ptMag=2) {
  par(mar=c(2,3,1,1))
  layout(matrix(c(1,2), 1, 2), widths=c(0.6, 0.4))
  plot(data$DE.DATE, data$drate, cex=ptMag*data$weight, ylim=c(-0.5, 0.2), )

  id <- identify(data$DE.DATE, data$drate, n=1, labels='')
  points(data$DE.DATE[id], data$drate[id],
         cex=ptMag*data$weight[id], pch=21, bg=2)
  RBS <- rbs(data, id)
  text(mean(par('usr')[c(1,2)]), par('usr')[4], pos=1,
       paste('Dive', id,
             '\nmax.res =', round(RBS$max.res, 2),
             '\nrRes=', round(data$rRes[id], 2),
             '\nweight=',round(data$weight[id], 2)))
}

