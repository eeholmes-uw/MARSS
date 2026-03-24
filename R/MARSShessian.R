# Attaches Hessian, parSigma and parMean to MLEobj
# Computed at the values in MLEobj$par
# For confidence intervals, this should be the MLEs
#' Parameter Variance-Covariance Matrix from the Hessian Matrix
#'
#' @description
#' Calculates an approximate parameter variance-covariance matrix for the parameters using an inverse of the Hessian of the negative log-likelihood function at the MLEs (the observed Fisher Information matrix). It appends `$Hessian`, `$parMean`, `$parSigma` to the [marssMLE] object.
#'
#' @param MLEobj An object of class [marssMLE]. This object must have a `$par` element containing MLE parameter estimates from e.g. [MARSSkem].
#' @param method The method to use for computing the Hessian. Options are `Harvey1989` to use the Harvey (1989) recursion, which is an analytical solution, `fdHess` or `optim` which are two numerical methods. Although `optim` can be passed to this function, in the internal functions which call this function, `fdHess` will be used if a numerical estimate is requested.
#'
#' @details
#' See [MARSSFisherI] for a discussion of the observed Fisher Information matrix and references.
#'
#' Method `fdHess` uses [nlme::fdHess] from package nlme to numerically estimate the Hessian matrix (the matrix of partial 2nd derivatives of the negative log-likelihood function at the MLE). Method `optim` uses [optim] with `hessian=TRUE` and `list(maxit=0)` to ensure that the Hessian is computed at the values in the `par` element of the MLE object. Method `Harvey1989` (the default) uses the recursion in Harvey (1989) to compute the observed Fisher Information of a MARSS model analytically.
#'
#' Note that the parameter confidence intervals computed with the observed Fisher Information matrix are based on the asymptotic normality of maximum-likelihood estimates under a large-sample approximation.
#'
#' @return
#' `MARSShessian()` attaches `Hessian`, `parMean` and `parSigma` to the [marssMLE] object that is passed into the function.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSFisherI()], [MARSSharveyobsFI()], [MARSShessian.numerical()], [MARSSparamCIs()], [marssMLE]
#'
#' @examples
#' dat <- t(harborSeal)
#' dat <- dat[c(2, 11), ]
#' MLEobj <- MARSS(dat)
#' MLEobj.hessian <- MARSShessian(MLEobj)
#'
#' # show the approx Hessian
#' MLEobj.hessian$Hessian
#'
#' # generate a parameter sample using the Hessian
#' # this uses the rmvnorm function in the mvtnorm package
#' hess.params <- mvtnorm::rmvnorm(1,
#'   mean = MLEobj.hessian$parMean,
#'   sigma = MLEobj.hessian$parSigma
#' )
#'
#' @references
#' Harvey, A. C. (1989) Section 3.4.5 (Information matrix) in Forecasting, structural time series models and the Kalman filter. Cambridge University Press, Cambridge, UK.
#'
#' See also J. E. Cavanaugh and R. H. Shumway (1996) On computing the expected Fisher information matrix for state-space model parameters. Statistics & Probability Letters 26: 347-355. This paper discusses the Harvey (1989) recursion (and proposes an alternative).
#' @export
MARSShessian <- function(MLEobj, method = c("Harvey1989", "fdHess", "optim")) {
  method <- match.arg(method)
  paramvec <- MARSSvectorizeparam(MLEobj)
  if (length(paramvec) == 0) stop("Stopped in MARSShessian(). No estimated parameter elements thus no Hessian.\n", call. = FALSE)
  if (is.null(MLEobj[["par"]])) {
    stop("Stopped in MARSShessian(). The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (identical(MLEobj[["convergence"]], 54)) {
    stop("Stopped in MARSShessian(). MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", call. = FALSE)
  }
  
  MLEobj$parMean <- paramvec
  MLEobj$Hessian <- MARSSFisherI(MLEobj, method = method)

  # When inverting the Hessian, need to deal with NAs in the Hessian
  Hess.tmp <- MLEobj$Hessian
  na.diag <- is.na(diag(Hess.tmp))
  if (any(na.diag)) {
    msg <- "MARSShessian: Hessian has NAs due to numerical errors. parSigma returned but some elements will be NA. See MARSSinfo(\"HessianNA\").\n"
    warning(msg)
    MLEobj$errors <- c(MLEobj$errors, msg)
  }
  # set diags with NA to 1 and non-diag to 0 so that Hessian can be inverted
  diag(Hess.tmp)[na.diag] <- 1
  Hess.tmp[is.na(Hess.tmp)] <- 0 #

  hessInv <- try(solve(Hess.tmp), silent = TRUE)
  if (inherits(hessInv, "try-error")) {
    msg <- "MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix. parSigma set to NULL.  See MARSSinfo(\"HessianNA\").\n"
    warning(msg)
    MLEobj$parSigma <- NULL
    MLEobj$errors <- c(MLEobj$errors, msg)
  } else {
    # Set NAs values to 0
    diag(hessInv)[na.diag] <- NA
    MLEobj$parSigma <- hessInv
  }

  return(MLEobj)
}
