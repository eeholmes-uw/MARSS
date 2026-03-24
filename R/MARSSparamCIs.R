####################################################################################
#   MARSSparamCIs function
#   This returns CIs for ML parameter estimates
#   If method='hessian', uses either Harvey1989 (analytical) or fdHess (numerical) or optim (numerical)
####################################################################################
#' @title Standard Errors, Confidence Intervals and Bias for MARSS Parameters
#'
#' @description
#' Computes standard errors, confidence intervals and bias for the
#' maximum-likelihood estimates of MARSS model parameters. If you want
#' confidence intervals on the estimated hidden states, see
#' [print.marssMLE()] and look for `states.cis`.
#'
#' @param MLEobj An object of class [marssMLE]. Must have a `$par` element
#'   containing the MLE parameter estimates.
#' @param method Method for calculating the standard errors: `"hessian"`,
#'   `"parametric"`, and `"innovations"` implemented currently.
#' @param alpha alpha level for the 1-alpha confidence intervals.
#' @param nboot Number of bootstraps to use for `"parametric"` and
#'   `"innovations"` methods.
#' @param hessian.fun The function to use for computing the Hessian. Options
#'   are `"Harvey1989"` (default analytical) or two numerical options:
#'   `"fdHess"` and `"optim"`. See [MARSShessian].
#' @param silent If false, a progress bar is shown for `"parametric"` and
#'   `"innovations"` methods.
#'
#' @details
#' Approximate confidence intervals (CIs) on the model parameters may be
#' calculated from the observed Fisher Information matrix ("Hessian CIs", see
#' [MARSSFisherI()]) or parametric or non-parametric (innovations) bootstrapping
#' using `nboot` bootstraps. The Hessian CIs are based on the asymptotic
#' normality of MLE parameters under a large-sample approximation. The Hessian
#' computation for variance-covariance matrices is a symmetric approximation
#' and the lower CIs for variances might be negative. Bootstrap estimates of
#' parameter bias are reported if method `"parametric"` or `"innovations"` is
#' specified.
#'
#' Note, these are added to the `par` elements of a [marssMLE] object but are
#' in `"marss"` form not `"marxss"` form. Thus the `MLEobj$par.upCI` and
#' related elements that are added to the [marssMLE] object may not look
#' familiar to the user. Instead the user should extract these elements using
#' `print(MLEobj)` and passing in the argument `what` set to `"par.se"`,
#' `"par.bias"`, `"par.lowCIs"`, or `"par.upCIs"`. See
#' [print.marssMLE][print](). Or use [tidy.marssMLE][tidy]().
#'
#' @return
#' `MARSSparamCIs` returns the [marssMLE] object passed in, with additional
#' components `par.se`, `par.upCI`, `par.lowCI`, `par.CI.alpha`,
#' `par.CI.method`, `par.CI.nboot` and `par.bias` (if method is `"parametric"`
#' or `"innovations"`).
#'
#' @references
#' Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of
#' multivariate time-series using the MARSS package. NOAA Fisheries, Northwest
#' Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112. Type
#' `RShowDoc("UserGuide", package = "MARSS")` to open a copy.
#'
#' @author Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSboot()], [MARSSinnovationsboot()], [MARSShessian()]
#'
#' @examples
#' dat <- t(harborSealWA)
#' dat <- dat[2:4, ]
#' kem <- MARSS(dat, model = list(
#'   Z = matrix(1, 3, 1),
#'   R = "diagonal and unequal"
#' ))
#' kem.with.CIs.from.hessian <- MARSSparamCIs(kem)
#' kem.with.CIs.from.hessian
#'
#' @export
MARSSparamCIs <- function(MLEobj, method = "hessian", alpha = 0.05, nboot = 1000, silent = TRUE, hessian.fun = "Harvey1989") {
  # this function expects a marssMLE object
  # it will add standard errors, biases, low/up CIs to MLEobj
  if (!inherits(MLEobj, "marssMLE")) {
    stop("MARSSparamCIs: This function needs a marssMLE object.\n", call. = FALSE)
  }
  if (!(method %in% c("hessian", "parametric", "innovations"))) stop("Stopped in MARSSparamCIs(). Current methods are hessian, parametric, innovations.\n", call. = FALSE)
  if (!(hessian.fun %in% c("Harvey1989", "fdHess", "optim"))) stop("Stopped in MARSSparamCIs(). Available functions for computing the Hessian are Harvey1989, fdHess, and optim.\n", call. = FALSE)
  if (is.null(MLEobj[["par"]])) {
    stop("Stopped in MARSSparamCIs(). The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (MLEobj[["convergence"]] == 54 && !(method=="parametric")) {
    stop("Stopped in MARSSparamCIs(). MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", call. = FALSE)
  }
  
  paramvec <- MARSSvectorizeparam(MLEobj)
  if (length(paramvec) == 0) stop("Stopped in MARSSparamCIs(). No estimated parameter elements.\n", call. = FALSE)

  paramnames <- names(paramvec)

  # default
  if (method == "hessian") {
    MLEobj <- MARSShessian(MLEobj, method = hessian.fun)

    vector.par.se <- rep(NA, length(paramvec))
    if (is.null(MLEobj$parSigma)) {
      warning("MARSSparamCIs: No parSigma element returned by Hessian function.  See marssMLE object errors (MLEobj$errors)")
    } else {
      is.neg <- diag(MLEobj$parSigma) < 0
      if (any(is.neg)) warning("MARSSparamCIs: parSigma element has negative values on the diagonal.  These are replaced with NA.")

      vector.par.se[!is.neg] <- try(sqrt(diag(MLEobj$parSigma)[!is.neg]), silent = TRUE) # wrap in a try() in case it fails
      if (inherits(vector.par.se, "try-error") || any(is.nan(vector.par.se))) {
        warning("MARSSparamCIs: Some of the Hessian diagonals cannot be square-rooted; NA is being returned")
      }
    }
    vector.par.upCI <- paramvec + qnorm(1 - alpha / 2) * vector.par.se
    vector.par.lowCI <- paramvec - qnorm(1 - alpha / 2) * vector.par.se
    vector.par.bias <- NULL
    par.CI.nboot <- NULL
  } # if method hessian

  if (method %in% c("parametric", "innovations")) {
    boot.params <- MARSSboot(MLEobj,
      nboot = nboot, output = "parameters", sim = method,
      param.gen = "MLE", silent = silent
    )$boot.params
    vector.par.lowCI <- apply(boot.params, 1, quantile, probs = alpha / 2)
    vector.par.upCI <- apply(boot.params, 1, quantile, probs = 1 - alpha / 2)
    vector.par.se <- sqrt(apply(boot.params, 1, var))
    paramvec <- MARSSvectorizeparam(MLEobj)
    vector.par.bias <- paramvec - apply(boot.params, 1, mean)
    par.CI.nboot <- nboot
  }

  # Finalize the output
  names(vector.par.se) <- paramnames
  if (!is.null(vector.par.bias)) names(vector.par.bias) <- paramnames
  names(vector.par.upCI) <- paramnames
  names(vector.par.lowCI) <- paramnames
  for (val in c("par.se", "par.bias", "par.upCI", "par.lowCI")) {
    tmp.vec <- get(paste("vector.", val, sep = ""))
    if (!is.null(tmp.vec)) {
      MLEobj[[val]] <- MARSSvectorizeparam(MLEobj, tmp.vec)[["par"]]
    } else {
      MLEobj[[val]] <- NULL
    }
  }

  # add on information about how the CIs were constructed
  if (method == "hessian") MLEobj$par.CI.info <- list(alpha = alpha, method = method, hessian.fun = hessian.fun)
  if (method %in% c("parametric", "innovations")) MLEobj$par.CI.info <- list(alpha = alpha, method = method, nboot = par.CI.nboot)
  return(MLEobj)
}
