#######################################################################################################
#   Parameter estimation using R's optim function
#   Minimal error checking is done.  You should run is.marssMLE() before calling this.
#   Q and R are not allowed to be time-varying
#   Likelihood computation is via  Kalman filter
#######################################################################################################
#' @title Parameter estimation for MARSS models using optim
#'
#' @description
#' Parameter estimation for MARSS models using R's [optim()] function. This
#' allows access to R's quasi-Newton algorithms available in that function. The
#' `MARSSoptim()` function is called when [MARSS()] is called with
#' `method="BFGS"`. This is an internal function in the [MARSS-package].
#'
#' @param MLEobj An object of class [marssMLE].
#'
#' @details
#' Objects of class [marssMLE] may be built from scratch but are easier to
#' construct using [MARSS()] called with `MARSS(..., fit=FALSE, method="BFGS")`.
#'
#' Options for [optim()] are passed in using `MLEobj$control`. See [optim()]
#' for a list of that function's control options. If `lower` and `upper` for
#' [optim()] need to be passed in, they should be passed in as part of
#' `control` as `control$lower` and `control$upper`. Additional `control`
#' arguments affect printing and initial conditions.
#'
#' * `MLEobj$control$kf.x0`: The initial condition is at \eqn{t=0} if
#'   `kf.x0="x00"`. The initial condition is at \eqn{t=1} if `kf.x0="x10"`.
#' * `MLEobj$marss$diffuse`: If `diffuse=TRUE`, a diffuse initial condition is
#'   used. `MLEobj$par$V0` is then the scaling function for the diffuse part of
#'   the prior. Thus the prior is `V0*kappa` where `kappa-->Inf`. Note that
#'   setting a diffuse prior does not change the correlation structure within
#'   the prior. If `diffuse=FALSE`, a non-diffuse prior is used and
#'   `MLEobj$par$V0` is the non-diffuse prior variance on the initial states.
#'   The prior is `V0`.
#' * `MLEobj$control$silent`: Suppresses printing of progress bars, error
#'   messages, warnings and convergence information.
#'
#' @return
#' The [marssMLE] object which was passed in, with additional components:
#'
#' * `method`: String `"BFGS"`.
#' * `kf`: Kalman filter output.
#' * `iter.record`: If `MLEobj$control$trace = TRUE`, then this is the
#'   `$message` value from [optim].
#' * `numIter`: Number of iterations needed for convergence.
#' * `convergence`: Did estimation converge successfully?
#'   * `convergence=0`: Converged in less than `MLEobj$control$maxit`
#'     iterations and no evidence of degenerate solution.
#'   * `convergence=3`: No convergence diagnostics were computed because all
#'     parameters were fixed thus no fitting required.
#'   * `convergence=-1`: No convergence diagnostics were computed because the
#'     MLE object was not fit (called with `fit=FALSE`). This isn't a
#'     convergence error just information. There is no `par` element so no
#'     functions can be run with the object.
#'   * `convergence=1`: Maximum number of iterations `MLEobj$control$maxit`
#'     was reached before `MLEobj$control$abstol` condition was satisfied.
#'   * `convergence=10`: Some of the variance elements appear to be degenerate.
#'   * `convergence=52`: The algorithm was abandoned due to errors from the
#'     `"L-BFGS-B"` method.
#'   * `convergence=53`: The algorithm was abandoned due to numerical errors in
#'     the likelihood calculation from [MARSSkf]. If this happens with
#'     `"BFGS"`, it can sometimes be helped with a better initial condition.
#'     Try using the EM algorithm first (`method="kem"`), and then using the
#'     parameter estimates from that as initial conditions for `method="BFGS"`.
#'   * `convergence=54`: The algorithm successfully fit the model but the
#'     Kalman filter/smoother could not be run on the model. Consult
#'     `MARSSinfo('optimerror54')` for insight.
#' * `logLik`: Log-likelihood.
#' * `states`: State estimates from the Kalman smoother.
#' * `states.se`: Confidence intervals based on state standard errors, see
#'   caption of Fig 6.3 (p. 337) in Shumway & Stoffer (2006).
#' * `errors`: Any error messages.
#'
#' @section Discussion:
#' The function only returns parameter estimates. To compute CIs, use
#' [MARSSparamCIs] but if you use parametric or non-parametric bootstrapping
#' with this function, it will use the EM algorithm to compute the bootstrap
#' parameter estimates! The quasi-Newton estimates are too fragile for the
#' bootstrap routine since one often needs to search to find a set of initial
#' conditions that work (i.e. don't lead to numerical errors).
#'
#' Estimates from `MARSSoptim` (which come from [optim]) should be checked
#' against estimates from the EM algorithm. If the quasi-Newton algorithm
#' works, it will tend to find parameters with higher likelihood faster than
#' the EM algorithm. However, the MARSS likelihood surface can be multimodal
#' with sharp peaks at degenerate solutions where a \eqn{\mathbf{Q}}{Q} or
#' \eqn{\mathbf{R}}{R} diagonal element equals 0. The quasi-Newton algorithm
#' sometimes gets stuck on these peaks even when they are not the maximum.
#' Neither an initial conditions search nor starting near the known maximum (or
#' from the parameters estimates after the EM algorithm) will necessarily solve
#' this problem. Thus it is wise to check against EM estimates to ensure that
#' the BFGS estimates are close to the MLE estimates (and vis-a-versa, it's
#' wise to rerun with `method="BFGS"` after using `method="kem"`). Conversely,
#' if there is a strong flat ridge in your likelihood, the EM algorithm can
#' report early convergence while the BFGS may continue much further along the
#' ridge and find very different parameter values. Of course a likelihood
#' surface with strong flat ridges makes the MLEs less informative...
#'
#' Note this is mainly a problem if the time series are short or very gappy. If
#' the time series are long, then the likelihood surface should be nice with a
#' single interior peak. In this case, the quasi-Newton algorithm works well
#' but it can still be sensitive (and slow) if not started with a good initial
#' condition. Thus starting it with the estimates from the EM algorithm is
#' often desirable.
#'
#' One should be aware that the prior set on the variance of the initial states
#' at t=0 or t=1 can have catastrophic effects on one's estimates if the
#' presumed prior covariance structure conflicts with the structure implied by
#' the MARSS model. For example, if you use a diagonal variance-covariance
#' matrix for the prior but the model implies a variance-covariance matrix with
#' non-zero covariances, your MLE estimates can be strongly influenced by the
#' prior variance-covariance matrix. Setting a diffuse prior does not help
#' because the diffuse prior still has the correlation structure specified by
#' V0. One way to detect priors effects is to compare the BFGS estimates to the
#' EM estimates. Persistent differences typically signify a problem with the
#' correlation structure in the prior conflicting with the implied correlation
#' structure in the MARSS model.
#'
#' @author Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSS()], [MARSSkem()], [marssMLE()], [optim()]
#'
#' @examples
#' dat <- t(harborSealWA)
#' dat <- dat[2:4, ] # remove the year row
#'
#' # fit a model with EM and then use that fit as the start for BFGS
#' # fit a model with 1 hidden state where obs errors are iid
#' # R="diagonal and equal" is the default so not specified
#' # Q is fixed
#' kemfit <- MARSS(dat, model = list(Z = matrix(1, 3, 1), Q = matrix(.01)))
#' bfgsfit <- MARSS(dat,
#'   model = list(Z = matrix(1, 3, 1), Q = matrix(.01)),
#'   inits = coef(kemfit, form = "marss"), method = "BFGS"
#' )
#'
#' @export
MARSSoptim <- function(MLEobj) {
  # This function does not check if user specified a legal MLE object.
  ##
  # Define needed negLogLik function with chol transformed variances
  #
  neglogLik <- function(x, MLEobj = NULL) { # NULL assignment needed for optim call syntax
    # MLEobj is tmp.MLEobj so has altered free and fixed
    # x is the paramvector

    # update the MLEobj by putting the estimated pars from optim in
    MLEobj <- MARSSvectorizeparam(MLEobj, x)
    free <- MLEobj$marss$free
    pars <- MLEobj$par
    par.dims <- attr(MLEobj[["marss"]], "model.dims")
    for (elem in c("Q", "R", "V0")) {
      if (!is.fixed(free[[elem]])) # recompute par if needed since par in parlist is transformed
        {
          d <- sub3D(free[[elem]], t = 1) # this will be the one with the upper tri zero-ed out
          par.dim <- par.dims[[elem]][1:2]
          # t=1 since D not allowed to be time-varying; since code 4 lines down won't work otherwise
          L <- unvec(d %*% pars[[elem]], dim = par.dim) # this by def will have 0 row/col at the fixed values
          the.par <- tcrossprod(L) # L%*%t(L)
          # from f+Dm=M and if f!=0, D==0 so can leave off f
          MLEobj$par[[elem]] <- solve(crossprod(d)) %*% t(d) %*% vec(the.par)
          # solve(t(d)%*%d)%*%t(d)%*%vec(the.par)
        }
    } # end for over elem
    # This function is passed a special MLEobj with a marss.original element
    MLEobj$marss$fixed <- MLEobj$fixed.original
    MLEobj$marss$free <- MLEobj$free.original

    # kfsel selects the Kalman filter / smoother function based on MLEobj$fun.kf
    negLL <- MARSSkf(MLEobj, only.logLik = TRUE, return.lag.one = FALSE)$logLik

    -1 * negLL
  }

  if (!inherits(MLEobj, "marssMLE")) {
    stop("Stopped in MARSSoptim(). Object of class marssMLE is required.\n", call. = FALSE)
  }
  for (elem in c("Q", "R")) {
    if (dim(MLEobj$model$free[[elem]])[3] > 1) {
      stop(paste("Stopped in MARSSoptim() because this function does not allow estimated part of ", elem, " to be time-varying.\n", sep = ""), call. = FALSE)
    }
  }
  # the is.marssMODEL call is.validvarcov() which tests that the blocks are diagonal or unconstrained in the varcov matrices

  ## attach would be risky here since user might have one of these variables in their workspace
  MODELobj <- MLEobj[["marss"]]
  y <- MODELobj$data # must have time going across columns
  free <- MODELobj$free
  fixed <- MODELobj$fixed
  tmp.inits <- MLEobj$start
  control <- MLEobj$control
  par.dims <- attr(MODELobj, "model.dims")
  m <- par.dims[["x"]][1]
  n <- par.dims[["y"]][1]

  ## Set up the control list for optim; only pass in optim control elements
  control.names <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
  optim.control <- list()
  for (elem in control.names) {
    if (!is.null(control[[elem]])) optim.control[[elem]] <- control[[elem]]
  }
  if (is.null(control[["lower"]])) {
    lower <- -Inf
  } else {
    lower <- control[["lower"]]
  }
  if (is.null(control[["upper"]])) {
    upper <- Inf
  } else {
    upper <- control$upper
  }
  if (control$trace == -1) optim.control$trace <- 0

  # The code is used to set things up to use MARSSvectorizeparam to just select inits for the estimated parameters
  # Q=t(chol(Q)%*%chol(Q)); t(chol(Q)) has 0 in upper triangle
  tmp.MLEobj <- MLEobj
  # This is needed for the likelihood calculation
  tmp.MLEobj$fixed.original <- tmp.MLEobj$marss$fixed
  tmp.MLEobj$free.original <- tmp.MLEobj$marss$free

  tmp.MLEobj$par <- tmp.inits # set initial conditions for estimated parameters
  for (elem in c("Q", "R", "V0")) { # need the chol for these
    d <- sub3D(free[[elem]], t = 1) # free[[elem]] is required to be time constant
    f <- sub3D(fixed[[elem]], t = 1) # placeholder.  Need structure not actual values
    the.par <- unvec(f + d %*% tmp.inits[[elem]], dim = par.dims[[elem]][1:2])
    is.zero <- diag(the.par) == 0 # where the 0s on diagonal are
    if (any(is.zero)) diag(the.par)[is.zero] <- 1 # so the chol doesn't fail if there are zeros on the diagonal
    the.par <- t(chol(the.par)) # convert to transpose of chol
    if (any(is.zero)) diag(the.par)[is.zero] <- 0 # set back to 0
    if (!is.fixed(free[[elem]])) {
      # from f+Dm=M so m = solve(crossprod(d))%*%t(d)%*%(vec(the.par)-f)
      # but if d!=0,then f==0. if f!-0, then d==0.
      # Thus crossprod(d))%*%t(d) has 0 cols where fs appear in the.par and f is not needed
      tmp.MLEobj$par[[elem]] <- solve(crossprod(d)) %*% t(d) %*% vec(the.par)
    } else {
      tmp.MLEobj$par[[elem]] <- matrix(0, 0, 1)
    }
    # when being passed to optim, pars for var-cov mat is the chol, so need to reset free so we can get the L=t(chol) matrix
    # we don't need to reset fixed because it won't be used;
    # step 1, compute the D matrix corresponding to upper.tri=0 in t(chol)
    # note this only works because it is required that
    # a) if f!=0, d=0. so 1+a never appears in var-cov mat  b) a+b never appears in a var-cov mat, c) for BFGS, var-cov mat is time-invariant
    tmp.list.mat <- fixed.free.to.formula(f, d, par.dims[[elem]][1:2])
    tmp.list.mat[upper.tri(tmp.list.mat)] <- 0 # set upper tri to zero
    tmp.MLEobj$marss$free[[elem]] <- convert.model.mat(tmp.list.mat)$free
  }
  # will return the inits only for the estimated parameters
  pars <- MARSSvectorizeparam(tmp.MLEobj)

  if (substr(tmp.MLEobj$method, 1, 4) == "BFGS") {
    optim.method <- "BFGS"
  } else {
    optim.method <- "something wrong"
  }

  kf.function <- MLEobj$fun.kf # used for printing
  optim.output <- try(optim(pars, neglogLik, MLEobj = tmp.MLEobj, method = optim.method, lower = lower, upper = upper, control = optim.control, hessian = FALSE), silent = TRUE)

  if (inherits(optim.output, "try-error")) { # try MARSSkfss if the user did not use it
    if (MLEobj$fun.kf != "MARSSkfss") { # if user did not request MARSSkf
      cat("MARSSkfas returned error.  Trying MARSSkfss.\n")
      tmp.MLEobj$fun.kf <- "MARSSkfss"
      kf.function <- "MARSSkfss" # used for printing
      optim.output <- try(optim(pars, neglogLik, MLEobj = tmp.MLEobj, method = optim.method, lower = lower, upper = upper, control = optim.control, hessian = FALSE), silent = TRUE)
    }
  }

  # error returned
  if (inherits(optim.output, "try-error")) {
    optim.output <- list(convergence = 53, message = c("MARSSkfas and MARSSkfss tried to compute log likelihood and encountered numerical problems.\n", sep = ""))
  }


  MLEobj.return <- MLEobj
  MLEobj.return$iter.record <- optim.output$message
  #   MLEobj.return$control=MLEobj$control
  #   MLEobj.return$model=MLEobj$model
  MLEobj.return$start <- tmp.inits # set to what was used here
  MLEobj.return$convergence <- optim.output$convergence
  if (optim.output$convergence %in% c(1, 0)) {
    if ((!control$silent || control$silent == 2) && optim.output$convergence == 0) cat(paste("Success! Converged in ", optim.output$counts[2], " iterations.\n", "Function ", kf.function, " used for likelihood calculation.\n", sep = ""))
    if ((!control$silent || control$silent == 2) && optim.output$convergence == 1) cat(paste("Warning! Max iterations of ", control$maxit, " reached before convergence.\n", "Function ", kf.function, " used for likelihood calculation.\n", sep = ""))

    tmp.MLEobj <- MARSSvectorizeparam(tmp.MLEobj, optim.output$par)
    # par has the fixed and estimated values using t chol of Q and R

    # back transform Q, R and V0 if needed from chol form to usual form
    for (elem in c("Q", "R", "V0")) { # this works because by def fixed and free blocks of var-cov mats are independent
      if (!is.fixed(MODELobj$free[[elem]])) # get a new par if needed
        {
          d <- sub3D(tmp.MLEobj$marss$free[[elem]], t = 1) # this will be the one with the upper tri zero-ed out but ok since symmetric
          par.dim <- par.dims[[elem]][1:2]
          L <- unvec(tmp.MLEobj$marss$free[[elem]][, , 1] %*% tmp.MLEobj$par[[elem]], dim = par.dim) # this by def will have 0 row/col at the fixed values
          the.par <- tcrossprod(L) # L%*%t(L)
          tmp.MLEobj$par[[elem]] <- solve(crossprod(d)) %*% t(d) %*% vec(the.par)
        }
    } # end for

    pars <- MARSSvectorizeparam(tmp.MLEobj) # now the pars values have been adjusted back to normal scaling
    # now put the estimated values back into the original MLEobj; fixed and free matrices as in original
    MLEobj.return <- MARSSvectorizeparam(MLEobj.return, pars)
    kf.out <- try(MARSSkf(MLEobj.return), silent = TRUE)

    if (inherits(kf.out, "try-error")) {
      MLEobj.return$numIter <- optim.output$counts[2]
      MLEobj.return$logLik <- -1 * optim.output$value
      MLEobj.return$errors <- c(paste0("\nWARNING: optim() successfully fit the model but ", kf.function, " returned an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", sep = ""), "\nError: ", kf.out[1])
      MLEobj.return$convergence <- 54
      MLEobj.return <- MARSSaic(MLEobj.return)
      kf.out <- NULL
    }
  } else {
    if (optim.output$convergence == 10) optim.output$message <- c("degeneracy of the Nelder-Mead simplex\n", paste("Function ", kf.function, " used for likelihood calculation.\n", sep = ""), optim.output$message)
    optim.output$counts <- NULL
    if (!control$silent) cat("MARSSoptim() stopped with errors. No parameter estimates returned.\n")
    if (control$silent == 2) cat("MARSSoptim() stopped with errors. No parameter estimates returned. See $errors in output for details.\n")

    MLEobj.return$par <- NULL
    MLEobj.return$errors <- optim.output$message
    kf.out <- NULL
  }

  if (!is.null(kf.out)) {
    if (control$trace > 0) MLEobj.return$kf <- kf.out
    MLEobj.return$states <- kf.out$xtT
    MLEobj.return$numIter <- optim.output$counts[2]
    MLEobj.return$logLik <- kf.out$logLik
  }
  MLEobj.return$method <- MLEobj$method

  ## Add AIC and AICc to the object
  if (!is.null(kf.out)) MLEobj.return <- MARSSaic(MLEobj.return)

  return(MLEobj.return)
}
