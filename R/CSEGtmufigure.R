#' Plot Forecast Uncertainty
#'
#' @description
#' Plot the uncertainty in the probability of hitting a percent threshold
#' (quasi-extinction) for a single random walk trajectory. This is the
#' quasi-extinction probability used in Population Viability Analysis. The
#' uncertainty is shown as a function of the forecast, where the forecast is
#' defined in terms of the forecast length (number of time steps) and
#' forecasted decline (percentage). This is a function used by one of the
#' vignettes in the [MARSS-package].
#'
#' @param N Time steps between the first and last population data point (positive integer)
#' @param u Per time-step decline (-0.1 means a 10\% decline per time step; 1 means a doubling per time step.)
#' @param s2p Process variance (Q). (a positive number)
#' @param make.legend Add a legend to the plot? (T/F)
#'
#' @details
#' This figure shows the region of high uncertainty in dark grey. In this
#' region, the minimum 95 percent confidence intervals on the probability of
#' quasi-extinction span 80 percent of the 0 to 1 probability. Green hashing
#' indicates where the 95 percent upper bound does not exceed 5\% probability
#' of quasi-extinction. The red hashing indicates, where the 95 percent lower
#' bound is above 95\% probability of quasi-extinction. The light grey lies
#' between these two certain/uncertain extremes. The extinction calculation is
#' based on Dennis et al. (1991). The minimum theoretical confidence interval
#' is based on Fieberg and Ellner (2000). This figure was developed in Ellner
#' and Holmes (2008).
#'
#' Examples using this figure are shown in the
#' [User Guide](https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf)
#' in the PVA chapter.
#'
#' @references
#' Dennis, B., P. L. Munholland, and J. M. Scott. (1991) Estimation of growth
#' and extinction parameters for endangered species. Ecological Monographs
#' 61:115-143.
#'
#' Fieberg, J. and Ellner, S.P. (2000) When is it meaningful to estimate an
#' extinction probability? Ecology, 81, 2040-2047.
#'
#' Ellner, S. P. and E. E. Holmes. (2008) Resolving the debate on when
#' extinction risk is predictable. Ecology Letters 11:E1-E5.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA, and Steve Ellner, Cornell Univ.
#'
#' @seealso [CSEGriskfigure()]
#'
#' @examples
#' CSEGtmufigure(N = 20, u = -0.1, s2p = 0.01)
#'
#' @keywords experimental hplot
#' @export
# This is a figure of the theoretical minimum uncertainty regions
# sensu Ellner and Holmes figure 1
# Code written by Steven Ellner and Eli Holmes
# The figure function
CSEGtmufigure <- function(N = 20, u = -0.1, s2p = 0.01, make.legend = TRUE) {
  mu <- as.vector(u) # function requires scalar not 1x1 matrix
  sigma2.b <- as.vector(s2p)
  if (s2p == 0) stop("Stopped in CSEGtmufigure(): function does not work with s2p=0.\n", call. = FALSE)
  # Set up some figure parameters
  ngrid <- 100
  Tvals <- seq(1, 110, length = ngrid)
  # win.graph(6,6); par(mfrow=c(1,1),cex.axis=1.35,cex.lab=1.35,yaxs="i",xaxs="i");
  par(yaxs = "i", xaxs = "i")

  # Define a number of functions
  Pext <- function(U, V) {
    pnorm(U - V) + exp(2 * U * V + pnorm(-(U + V), log.p = T))
  }

  TMU <- function(T, RelNe, side = "one") {
    qt <- switch(EXPR = side, one = qnorm(0.95), two = qnorm(0.975))
    a <- log(1 / RelNe)
    U <- -mu * sqrt(T / sigma2.b)
    V <- a / (sqrt(sigma2.b * T))
    upci <- Pext(U + qt * sqrt(T / N), V)
    lowci <- Pext(U - qt * sqrt(T / N), V)
    return(list(upper = upci, lower = lowci))
  }

  TMUpper <- function(T, RelNe, side = "one") TMU(T, RelNe, side)$upper
  TMLower <- function(T, RelNe, side = "one") TMU(T, RelNe, side)$lower

  CIWidth <- function(T, RelNe) {
    qt <- qnorm(0.975)
    a <- log(1 / RelNe)
    U <- -mu * sqrt(T / sigma2.b)
    V <- a / (sqrt(sigma2.b * T))
    upci <- Pext(U + qt * sqrt(T / N), V)
    lowci <- Pext(U - qt * sqrt(T / N), V)
    return(upci - lowci)
  }

  xlabs <- expression(paste("Projection interval ", italic(T), " time steps"))
  ylabs <- "xe = log10(N0/Ne)"
  safe.limits <- numeric(ngrid)
  dead.limits <- numeric(ngrid)
  minval <- 1e-16
  maxval <- 1 - minval
  for (j in 1:ngrid) {
    T <- Tvals[j]
    safe.limits[j] <- uniroot(function(x) {
      TMUpper(T, x) - 0.05
    }, lower = minval, upper = maxval)$root
    dead.limits[j] <- uniroot(function(x) {
      TMLower(T, x) - 0.95
    }, lower = minval, upper = maxval)$root
  }

  graphics::matplot(Tvals, log10(cbind(safe.limits, dead.limits)),
    ylim = c(-2.2, 0), type = "l", lty = 1, col = "white",
    xlim = c(1, 100), xlab = xlabs, ylab = ylabs
  )
  polygon(c(Tvals, rev(Tvals)), log10(c(safe.limits, rev(dead.limits))), col = "grey85")
  polygon(c(Tvals, rev(Tvals)), log10(c(dead.limits, rep(1, ngrid))), col = "black")
  polygon(c(Tvals, rev(Tvals)), c(log10(safe.limits), rep(-3, ngrid)), col = "white")

  logRelNe <- seq(-3.125, -0.001, length = ngrid)
  RelNe <- 10^(logRelNe)
  CIW <- outer(Tvals, RelNe, CIWidth)
  z <- contourLines(Tvals, log10(RelNe), CIW, levels = 0.8)
  if (length(z) != 0) {
    xvals <- z[[1]]$x
    yvals <- z[[1]]$y
    polygon(c(xvals, 100), c(yvals, yvals[1]), col = "grey45")
  }

  abline(h = -2.2)
  abline(h = 0)
  abline(v = 100)
  abline(v = 1)
  offset <- -0.05
  abline(h = -.3)
  text(5, -.3 + offset, "50%")
  abline(h = -1)
  text(5, -1 + offset, "90%")
  abline(h = -2)
  text(5, -2 + offset, "99%")
  # title(expression(paste("nyrs = ", N, ", hat(mu), = ", mu, sigma[b]," = ",sigma2.b)));
  title(paste("time steps = ", N, "\n mu = ", format(mu, digits = 2), "s2.p = ", format(sigma2.b, digits = 2)))
  if (make.legend) legend("topright", inset = 0.02, bg = "white", c("high certainty P<0.05", "high certainty P>0.95", "uncertain", "highly uncertain"), fill = c("white", "black", "grey85", "grey45"))
}
