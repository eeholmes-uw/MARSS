#######################################################################################################
#   MARSSaic function
#   Adds to marssMLE object:
#   elements in output arg
#   samp.size, num.params
#######################################################################################################

#' AIC for MARSS Models
#'
#' @description
#' Calculates AIC, AICc, a parametric bootstrap AIC (AICbp) and a
#' non-parametric bootstrap AIC (AICbb). If you simply want the AIC value for
#' a [marssMLE] object, you can use `AIC(fit)`.
#'
#' @param MLEobj An object of class [marssMLE]. This object must have a `$par`
#'   element containing MLE parameter estimates from e.g. `MARSSkem()`.
#' @param output A vector containing one or more of the following: "AIC",
#'   "AICc", "AICbp", "AICbb", "AICi", "boot.params". See Details.
#' @param Options A list containing:
#'   * `nboot` Number of bootstraps (positive integer)
#'   * `return.logL.star` Return the log-likelihoods for each bootstrap? (T/F)
#'   * `silent` Suppress printing of the progress bar during AIC bootstraps?
#'     (T/F)
#'
#' @details
#' When sample size is small, Akaike's Information Criterion (AIC)
#' under-penalizes more complex models. The most commonly used small sample
#' size corrector is AICc, which uses a penalty term of
#' \eqn{K n/(n-K-1)}, where \eqn{K} is the number of estimated parameters.
#' However, for time series models, AICc still under-penalizes complex models;
#' this is especially true for MARSS models.
#'
#' Two small-sample estimators specific for MARSS models have been developed.
#' Cavanaugh and Shumway (1997) developed a variant of bootstrapped AIC using
#' Stoffer and Wall's (1991) bootstrap algorithm ("AICbb"). Holmes and Ward
#' (2010) developed a variant on AICb ("AICbp") using a parametric bootstrap.
#' The parametric bootstrap permits AICb calculation when there are missing
#' values in the data, which Cavanaugh and Shumway's algorithm does not allow.
#' More recently, Bengtsson and Cavanaugh (2006) developed another
#' small-sample AIC estimator, AICi, based on fitting candidate models to
#' multivariate white noise.
#'
#' When the `output` argument passed in includes both `"AICbp"` and
#' `"boot.params"`, the bootstrapped parameters from `"AICbp"` will be added
#' to `MLEobj`.
#'
#' @return
#' Returns the [marssMLE] object that was passed in with additional AIC
#' components added on top as specified in the 'output' argument.
#'
#' @references
#' Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of
#' multivariate time-series using the MARSS package. NOAA Fisheries, Northwest
#' Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112. Type
#' `RShowDoc("UserGuide",package="MARSS")` to open a copy.
#'
#' Bengtsson, T., and J. E. Cavanaugh. 2006. An improved Akaike information
#' criterion for state-space model selection. Computational Statistics & Data
#' Analysis 50:2635-2654.
#'
#' Cavanaugh, J. E., and R. H. Shumway. 1997. A bootstrap variant of AIC for
#' state-space model selection. Statistica Sinica 7:473-496.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSboot()]
#'
#' @examples
#' dat <- t(harborSealWA)
#' dat <- dat[2:3, ]
#' kem <- MARSS(dat, model = list(
#'   Z = matrix(1, 2, 1),
#'   R = "diagonal and equal"
#' ))
#' kemAIC <- MARSSaic(kem, output = c("AIC", "AICc"))
#'
#' @export
MARSSaic <- function(MLEobj, output = c("AIC", "AICc"), Options = list(
                       nboot = 1000, return.logL.star = FALSE,
                       silent = FALSE
                     )) {
  # Options$nboot is the number of bootstrap replicates to do
  # mssm.model is a specified model (needs the structure and parameter elements of the list)
  # output tells what output to produce; this can be a vector  if multiple items should be returned
  #      AIC, AICc, AICbp, AICbb, AICi, boot.params
  # silent is a flag to indicate whether the progress bar should be printed

  if (!is.list(Options)) {
    msg <- " Options argument must be passed in as a list.\n"
    cat("\nErrors were caught in MARSSaic \n", msg, sep = "")
    stop("Stopped in MARSSaic() due to problem(s) with arguments.\n", call. = FALSE)
  }
  ## Set options if some not passed in
  if (is.null(Options[["nboot"]])) Options$nboot <- 1000
  if (is.null(Options[["return.logL.star"]])) Options$return.logL.star <- FALSE
  if (is.null(Options[["silent"]])) Options$silent <- FALSE
  if (!inherits(MLEobj, "marssMLE")) {
    stop("Stopped in MARSSaic(). An object of class marssMLE is required.\n", call. = FALSE)
  }

  if (is.null(MLEobj[["logLik"]])) {
    msg <- " No log likelihood.  This function expects a model fitted via maximum-likelihood.\n"
    cat("\n", "Errors were caught in MARSSaic \n", msg, sep = "")
    stop("Stopped in MARSSaic() due to problem(s) with the MLE object passed in.\n", call. = FALSE)
  }
  return.list <- list()

  ## Some renaming for readability
  model <- MLEobj$marss
  # kf = MLEobj$kf
  loglike <- MLEobj$logLik

  ##### AIC and AICc calculations
  if ("AIC" %in% output || "AICc" %in% output) {
    K <- 0
    for (elem in c("Z", "A", "B", "U", "x0", "R", "Q", "V0")) {
      pars <- dim(model$free[[elem]])[2]
      K <- K + pars
    }
    MLEobj$AIC <- -2 * loglike + 2 * K
    samp.size <- sum(!is.na(model$data))
    MLEobj$AICc <- ifelse(samp.size > (K + 1), -2 * loglike + 2 * K * (samp.size / (samp.size - K - 1)), "NA, number of data points less than K+1")
    MLEobj$samp.size <- samp.size
    MLEobj$num.params <- K
  }

  ##### AICbb & AICbp
  method <- c("AICbb", "AICbp")

  for (m in method[which(method %in% output)]) {
    bootstrap.method <- switch(m,
      AICbb = "innovations",
      AICbp = "parametric"
    )
    logL.star <- 0

    drawProgressBar <- FALSE # If the time library is not installed, no prog bar
    if (!Options$silent) { # then we can draw a progress bar
      cat(paste(m, "calculation in progress...\n"))
      prev <- progressBar()
      drawProgressBar <- TRUE
    }

    boot.params <- MARSSboot(MLEobj,
      nboot = Options$nboot, output = "parameters", sim = bootstrap.method,
      param.gen = "MLE", silent = TRUE
    )$boot.params

    if ((bootstrap.method == "parametric") && ("boot.params" %in% output)) MLEobj$boot.params <- boot.params

    for (i in 1:Options$nboot) {
      boot.model <- MARSSvectorizeparam(MLEobj, parvec = boot.params[, i]) # boot.model is a MLEobj
      logL.star[i] <- MARSSkf(boot.model, only.logLik = TRUE)$logLik
      if (drawProgressBar) prev <- progressBar(i / Options$nboot, prev)
    }
    MLEobj[[m]] <- -4 * (sum(logL.star)) / Options$nboot + 2 * loglike # -2*model$loglik + 2*(1/N)*(-2)*sum(boot$Yloglike-model$logLik)
    if (Options$return.logL.star == TRUE) {
      tmp <- paste(m, ".logL.star", sep = "")
      MLEobj[[tmp]] <- logL.star
    }
  }

  return(MLEobj)
}
