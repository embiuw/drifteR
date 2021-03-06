% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitDrate.R
\name{fitPlot}
\alias{fitPlot}
\title{\code{fitPlot} Plot output from \code{fitDrift}}
\usage{
fitPlot(
  fit,
  yLims = c(-0.5, 0.2),
  fitCol = 2,
  fitAlpha = 0.2,
  add = F,
  fitLwd = 1,
  ptMag = 2,
  ptAlpha = 1
)
}
\arguments{
\item{fit}{Output list from the \code{fitDrift} function.}

\item{yLims}{Set range of y axis}

\item{fitCol}{Colour to be used for estimated drift trajectory}

\item{fitLwd}{Line width to be used for estimated mean drift trajectory}

\item{ptMag}{Scaling for point size}

\item{ptAlpha}{Scaling for point transparency (usually has to be 1 or less)}

\item{fitAlpa}{Transparency to e used for estimated drift trajectory}
}
\value{
A plot with data and the fitted estimates with confidence intervals
}
\description{
\code{fitPlot} Plot output from \code{fitDrift}
}
\examples{
rbs(data=dive, num=100, n.bs=4)
dive <- rbs(data=dive, num=NA, n.bs=4)
dive <- pDrift(dive)
fit <- fitDrate(dive, prefilter=0.01)
fit2 <- fitDrate(dive, weight='rRes*span*dWt', prefilter=0.01)
fitPlot(fit)
fitPlot(fit2, add=T, fitCol=3)
}
\seealso{
\code{\link{rbs}} Reverse Broken Stick algorithm

\code{\link{check.pDrift}} to interactively check drift dive weightings

\code{\link{pDrift}} to calculate drift rate and probability weightings for drift dive detection

\code{\link{fitDrate}} To fit state-space model to weighted data to estimate drift rate change over time

Other Drift dive functions: 
\code{\link{check.pDrift}()},
\code{\link{fitDrate}()},
\code{\link{pDrift}()},
\code{\link{rbs}()}
}
\author{
Martin Biuw
}
\concept{Drift dive functions}
