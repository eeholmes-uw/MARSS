########################################################################################################################
#   MARSSinnovationsboot.r
#   This bootstraps multivariate time series data using Stoffer and Wall algorithm
#   It creates bootstrap data via sampling from the standardized innovations matrix
#   In the MARSS code, this is referred to as the nonparametric bootstrap.  Strictly speaking, it is not nonparametric.
########################################################################################################################
#' Bootstrapped Data using Stoffer and Wall's Algorithm
#'
#' @description
#' Creates bootstrap data via sampling from the standardized innovations matrix.
#' This is an internal function in the **MARSS** package and is not exported.
#' Users should access this with [MARSSboot()].
#'
#' @param MLEobj An object of class [marssMLE]. This object must have a `$par`
#'   element containing MLE parameter estimates from e.g. [MARSSkem()] or
#'   [MARSS()]. This algorithm may not be used if there are missing datapoints
#'   in the data.
#' @param nboot Number of bootstraps to perform.
#' @param minIndx Number of innovations to skip. Stoffer & Wall suggest not
#'   sampling from innovations 1-3.
#'
#' @details
#' Stoffer and Wall (1991) present an algorithm for generating CIs via a
#' non-parametric bootstrap for state-space models. The basic idea is that the
#' Kalman filter can be used to generate estimates of the residuals of the model
#' fit. These residuals are then standardized and resampled and used to generate
#' bootstrapped data using the MARSS model and its maximum-likelihood parameter
#' estimates. One of the limitations of the Stoffer and Wall algorithm is that
#' it cannot be used when there are missing data, unless all data at time
#' \eqn{t} are missing.
#'
#' @return
#' A list containing the following components:
#'
#' * `boot.states`: Array (dim is m x tSteps x nboot) of simulated state processes.
#' * `boot.data`: Array (dim is n x tSteps x nboot) of simulated data.
#' * `marss`: [marssMODEL] object element of the [marssMLE] object (`marssMLE$marss`) in "marss" form.
#' * `nboot`: Number of bootstraps performed.
#'
#' m is the number state processes (x in the MARSS model) and n is the number
#' of observation time series (y in the MARSS model).
#'
#' @references
#' Stoffer, D. S., and K. D. Wall. 1991. Bootstrapping state-space models:
#' Gaussian maximum likelihood estimation and the Kalman filter. Journal of
#' the American Statistical Association 86:1024-1033.
#'
#' @author
#' Eli Holmes and Eric Ward, NOAA, Seattle, USA.
#'
#' @seealso [stdInnov()], [MARSSparamCIs()], [MARSSboot()]
#'
#' @examples
#' dat <- t(kestrel)
#' dat <- dat[2:3, ]
#' fit <- MARSS(dat, model = list(U = "equal", Q = diag(.01, 2)))
#' boot.obj <- MARSSinnovationsboot(fit)
#' @export
MARSSinnovationsboot <- function(MLEobj, nboot = 1000, minIndx = 3) {
  if (any(is.na(MLEobj$marss$data))) {
    stop("Stopped in MARSSinnovationsboot() because this algorithm resamples from the innovations and doesn't allow missing values.\n", call. = FALSE)
  }

  # The following need to be present in model: parameter estimates + Sigma, Kt, Innov; the latter 3 are part of $kf
  if (any(is.null(MLEobj[["par"]]), is.null(MLEobj[["marss"]][["data"]]))) {
    stop("Stopped in MARSSinnovationsboot(). The passed in marssMLE object is missing par or model$data elements.\n", call. = FALSE)
  }
  if (is.null(MLEobj[["kf"]]) || !is.array(MLEobj[["kf"]][["Innov"]]) || !is.array(MLEobj[["kf"]][["Sigma"]])) {
    kf <- MARSSkfss(MLEobj)
  } else {
    kf <- MLEobj[["kf"]]
  }

  ## Rename things for code readability
  MODELobj <- MLEobj[["marss"]]
  par.dims <- attr(MODELobj, "model.dims")
  TT <- par.dims[["data"]][2]
  m <- par.dims[["x"]][1]
  n <- par.dims[["y"]][1]
  f <- MODELobj[["fixed"]]
  d <- MODELobj[["free"]]
  pari <- parmat(MLEobj, t = 1)

  #### make a list of time-varying parameters
  time.varying <- list()
  for (elem in attr(MODELobj, "par.names")) {
    if ((dim(d[[elem]])[3] == 1) & (dim(f[[elem]])[3] == 1)) {
      time.varying[[elem]] <- FALSE
    } else {
      time.varying[[elem]] <- TRUE
    } # not time-varying
  }

  ##### Set holders for output
  newData <- matrix(NA, n, TT) # store as written for state-space eqn with time across columns
  newStates <- matrix(NA, m, TT + 1) # store as written for state-space eqn with time across columns
  boot.data <- array(NA, dim = c(par.dims[["data"]], nboot))
  boot.states <- array(NA, dim = c(par.dims[["x"]], nboot))

  # Calculate the sqrt of sigma matrices, so they don't have to be computed 5000+ times
  sigma.Sqrt <- array(0, dim = c(n, n, TT))
  BKS <- array(0, dim = c(m, n, TT))
  for (i in 1:TT) {
    std.innovations <- stdInnov(kf$Sigma, kf$Innov) # standardized innovations; time across rows
    eig <- eigen(kf$Sigma[, , i]) # calculate inv sqrt(sigma[1])
    sigma.Sqrt[, , i] <- eig$vectors %*% tcrossprod(makediag(sqrt(eig$values)), eig$vectors)
    if (time.varying$B & i > 1) pari$B <- parmat(MLEobj, "B", t = i)$B
    BKS[, , i] <- pari$B %*% kf$Kt[, , i] %*% sigma.Sqrt[, , i] # this is m x n
  }

  for (i in 1:nboot) {
    # Then the bootstrap algorithm proceeds
    # Stoffer & Wall suggest not sampling from innovations 1-3 (not much data)
    minIndx <- ifelse(TT > 5, minIndx, 1)
    samp <- sample(seq(minIndx + 1, TT), size = (TT - minIndx), replace = TRUE)
    e <- as.matrix(std.innovations)
    e[, 1:minIndx] <- as.matrix(std.innovations[, 1:minIndx])
    e[, (minIndx + 1):TT] <- as.matrix(std.innovations[, samp])
    # newStates is a[] in the writeup by EH
    newStates[, 1] <- pari$x0
    # newStates[,2] = B %*% x0 + U + B %*% Kt[,,1] %*% sigma.Sqrt[,,1] %*% e[,1]
    for (t in 1:TT) {
      if (time.varying$B & i > 1) pari$B <- parmat(MLEobj, "B", t = i)$B
      if (time.varying$U & i > 1) pari$U <- parmat(MLEobj, "U", t = i)$U
      if (time.varying$Z & i > 1) pari$Z <- parmat(MLEobj, "Z", t = i)$Z
      if (time.varying$A & i > 1) pari$A <- parmat(MLEobj, "A", t = i)$A
      newStates[, t + 1] <- pari$B %*% newStates[, t, drop = FALSE] + pari$U + BKS[, , t] %*% e[, t, drop = FALSE]
      newData[, t] <- pari$Z %*% newStates[, t, drop = FALSE] + pari$A + sigma.Sqrt[, , t] %*% e[, t, drop = FALSE]
    }

    newStates <- newStates[, 2:(TT + 1)]
    boot.data[, , i] <- newData
    boot.states[, , i] <- newStates
    # reset newStates to its original dim
    newStates <- matrix(NA, m, TT + 1)
  }
  return(list(boot.states = boot.states, boot.data = boot.data, marss = MODELobj, nboot = nboot))
}

######################################################################################################################
#   stdInnov
######################################################################################################################
#' Standardized Innovations
#'
#' @description
#' Standardizes Kalman filter innovations. This is a helper function called by
#' [MARSSinnovationsboot()] in the **MARSS** package. Not exported.
#'
#' @param SIGMA n x n x T array of Kalman filter innovations variances. This is
#'   output from [MARSSkf()].
#' @param INNOV n x T matrix of Kalman filter innovations. This is output from
#'   [MARSSkf()].
#'
#' @details
#' n = number of observation (y) time series. T = number of time steps in the
#' time series.
#'
#' @return
#' n x T matrix of standardized innovations.
#'
#' @references
#' Stoffer, D. S., and K. D. Wall. 1991. Bootstrapping state-space models:
#' Gaussian maximum likelihood estimation and the Kalman filter. Journal of
#' the American Statistical Association 86:1024-1033.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSboot()], [MARSSkf()], [MARSSinnovationsboot()]
#'
#' @examples
#' \dontrun{
#' std.innovations <- stdInnov(kfList$Sigma, kfList$Innov)
#' }
#' @keywords internal
stdInnov <- function(SIGMA, INNOV) {
  # This function added by EW Nov 3, 2008
  # SIGMA is covariance matrix, E are original innovations
  TT <- dim(INNOV)[2]
  n <- dim(INNOV)[1]
  SI <- matrix(0, n, TT)
  for (i in 1:TT) {
    a <- SIGMA[, , i]
    a.eig <- eigen(a)
    a.sqrt <- a.eig$vectors %*% makediag((sqrt(a.eig$values))) %*% solve(a.eig$vectors)
    SI[, i] <- chol2inv(chol(a.sqrt)) %*% INNOV[, i] # eqn 1, p359 S&S
  }
  SI[which(INNOV == 0)] <- 0 # make sure things that should be 0 are 0
  return(SI)
}
