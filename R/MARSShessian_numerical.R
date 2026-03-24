#######################################################################################################
#   MARSShessian.numerical functions
#
#   Numerical estimation of the Hessian of the negative log-likelihood function
#     at the parameters values in MLEobj$par.  For the Hessian to be used
#     for confidence interval calculation, this should be the MLEs.
#     Two functions are available for the Hessian calculation, fdHess() and optim().
#
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#######################################################################################################
#' Hessian Matrix via Numerical Approximation
#'
#' @description
#' Calculates the Hessian of the log-likelihood function at the MLEs using either the [nlme::fdHess] function in the nlme package or the [optim] function. This is a utility function in the [MARSS-package] and is not exported. Use [MARSShessian] to access.
#'
#' @param MLEobj An object of class [marssMLE]. This object must have a `$par` element containing MLE parameter estimates from e.g. [MARSSkem].
#' @param fun The function to use for computing the Hessian. Options are 'fdHess' or 'optim'.
#'
#' @details
#' Method `fdHess` uses [nlme::fdHess] from package nlme to numerically estimate the Hessian matrix (the matrix of partial 2nd derivatives) of the negative log-likelihood function with respect to the parameters. Method `optim` uses [optim] with `hessian=TRUE` and `list(maxit=0)` to ensure that the Hessian is computed at the values in the `par` element of the MLE object.
#'
#' @return
#' The numerically estimated Hessian of the log-likelihood function at the maximum likelihood estimates.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSharveyobsFI()], [MARSShessian()], [MARSSparamCIs()]
#'
#' @examples
#' dat <- t(harborSeal)
#' dat <- dat[c(2, 11), ]
#' MLEobj <- MARSS(dat)
#' MARSS:::MARSShessian.numerical(MLEobj)
#' @keywords internal
MARSShessian.numerical <- function(MLEobj, fun = c("fdHess", "optim")) {
  fun <- match.arg(fun)
  kfNLL <- function(x, MLEobj = NULL) { # NULL assignment needed for optim call syntax
    MLEobj <- MARSSvectorizeparam(MLEobj, x)
    negLL <- MARSSkf(MLEobj, only.logLik = TRUE, return.lag.one = FALSE)$logLik
    -1 * negLL
  }

  paramvector <- MARSSvectorizeparam(MLEobj)

  # Hessian and gradient
  if (fun == "fdHess") {
    out <- fdHess(paramvector, function(paramvector, MLEobj) kfNLL(paramvector, MLEobj), MLEobj)
    Hessian <- out$Hessian
  }
  if (fun == "optim") {
    # maxit set to 0 so that the Hessian is computed at the values in the MLEobj
    out <- optim(paramvector, kfNLL, MLEobj = MLEobj, method = "BFGS", hessian = TRUE, control = list(maxit = 0))
    Hessian <- out$hessian
  }

  rownames(Hessian) <- names(paramvector)
  colnames(Hessian) <- names(paramvector)
  # MLEobj$gradient = out$gradient

  return(Hessian)
}
