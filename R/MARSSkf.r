#######################################################################################################
#   MARSSkf function
#   Utility function to choose the Kalman filter and smoother
#######################################################################################################
#' Kalman Filtering and Smoothing
#'
#' Provides Kalman filter and smoother output for MARSS models with (or without) time-varying parameters. `MARSSkf()` is a small helper function to select which Kalman filter/smoother function to use based on the value in `MLEobj$fun.kf`.  The choices are `MARSSkfas()` which uses the filtering and smoothing algorithms in the [KFAS](https://CRAN.R-project.org/package=KFAS) package based on algorithms in Koopman and Durbin (2001-2003), and `MARSSkfss()` which uses the algorithms in Shumway and Stoffer. The default function is `MARSSkfas()` which is faster and generally more stable (fewer matrix inversions), but there are some cases where `MARSSkfss()` might be more stable and it returns a variety of diagnostics that `MARSSkfas()` does not.
#'
#' @param MLEobj A [marssMLE()] object with the `par` element of estimated parameters, `marss` element with the model description (in marss form) and data, and `control` element for the fitting algorithm specifications.  `control$debugkf` specifies that detailed error reporting will be returned (only used by `MARSSkf()`).  `model$diffuse=TRUE` specifies that a diffuse prior be used (only used by `MARSSkfas()`). See [KFAS::KFS()] documentation. When the diffuse prior is set, `V0` should be non-zero since the diffuse prior variance is `V0*kappa`, where kappa goes to infinity.
#' @param smoother Used by `MARSSkfss()`.  If set to FALSE, only the Kalman filter is run. The output `xtT`, `VtT`, `x0T`, `Vtt1T`, `V0T`, and `J0` will be NULL.
#' @param only.logLik Used by `MARSSkfas()`.  If set, only the log-likelihood is returned using the [KFAS::KFAS()] package function [KFAS::logLik.SSModel()].  This is much faster if only the log-likelihood is needed.
#' @param return.lag.one Used by `MARSSkfas()`.  If set to FALSE, the smoothed lag-one covariance values are not returned (output `Vtt1T` is set to NULL).  This speeds up `MARSSkfas()` because to return the smoothed lag-one covariance a stacked MARSS model is used with twice the number of state vectors---thus the state matrices are larger and take more time to work with.
#' @param return.kfas.model Used by `MARSSkfas()`.  If set to TRUE, it returns the MARSS model in [KFAS::KFAS()] model form (class [KFAS::SSModel()]).  This is useful if you want to use other KFAS functions or write your own functions to work with [optim()] to do optimization.  This can speed things up since there is a bit of code overhead in [MARSSoptim()] associated with the [marssMODEL()] model specification needed for the constrained EM algorithm (but not strictly needed for [optim()]; useful but not required.).
#' @param newdata A new matrix of data to use in place of the data used to fit the model (in the `model$data` and `marss$data` elements of a [marssMLE()] object). If the initial \eqn{x} was estimated (in `x0`) then this estimate will be used for `newdata` and this may not be appropriate.
#'
#' @details
#' For state-space models, the Kalman filter and smoother provide optimal (minimum mean square error) estimates of the hidden states. The Kalman filter is a forward recursive algorithm which computes estimates of the states \eqn{\mathbf{x}_t}{x(t)} conditioned on the data up to time \eqn{t} (`xtt`). The Kalman smoother is a backward recursive algorithm which starts at time \eqn{T} and works backwards to \eqn{t = 1} to provide estimates of the states conditioned on all data (`xtT`).    The data may contain missing values (NAs).  All parameters may be time varying.
#'
#' The initial state is either an estimated parameter or treated as a prior (with mean and variance). The initial state can be specified at \eqn{t=0} or \eqn{t=1}.  The EM algorithm in the MARSS package ([MARSSkem()]) provides both Shumway and Stoffer's derivation that uses \eqn{t=0} and Ghahramani et al algorithm which uses \eqn{t=1}.  The `MLEobj$model$tinitx` argument specifies whether the initial states (specified with `x0` and `V0` in the `model` list) is at \eqn{t=0} (`tinitx=0`) or \eqn{t=1} (`tinitx=1`). If `MLEobj$model$tinitx=0`, `x0` is defined as \eqn{\textrm{E}[\mathbf{X}_0|\mathbf{y}_0]}{E[X(0)|y(0)]} and `V0` is defined as \eqn{\textrm{E}[\mathbf{X}_0\mathbf{X}_0|\mathbf{y}_0]}{E[X(0)X(0)|y(0)]} which appear in the Kalman filter at \eqn{t=1} (first set of equations). If `MLEobj$model$tinitx=1`, `x0` is defined as \eqn{\textrm{E}[\mathbf{X}_1|\mathbf{y}_0]}{E[X(1)|y(0)]} and `V0` is defined as \eqn{\textrm{E}[\mathbf{X}_1\mathbf{X}_1|\mathbf{y}_0]}{E[X(1)X(1)|y(0)]} which appear in the Kalman filter at \eqn{t=1} (and the filter starts at t=1 at the 3rd and 4th equations in the Kalman filter recursion). Thus if `MLEobj$model$tinitx=1`, `x0=xtt1[,1]` and `V0=Vtt1[,,1]` in the Kalman filter output while if `MLEobj$model$tinitx=0`, the initial condition will not be in the filter output since time starts at 1 not 0 in the output.
#'
#' `MARSSkfss()` is a native R implementation based on the Kalman filter and smoother equation as shown in Shumway and Stoffer (sec 6.2, 2006).  The equations have been altered to allow the initial state distribution to be to be specified at \eqn{t=0} or \eqn{t=1} (data starts at \eqn{t=1}) per per Ghahramani and Hinton (1996).  In addition, the filter and smoother equations have been altered to allow partially deterministic models (some or all elements of the \eqn{\mathbf{Q}}{Q} diagonals equal to 0), partially perfect observation models (some or all elements of the \eqn{\mathbf{R}}{R} diagonal equal to 0) and fixed (albeit unknown) initial states (some or all elements of the \eqn{\mathbf{V0}}{V0} diagonal equal to 0) (per Holmes 2012).  The code includes numerous checks to alert the user if matrices are becoming ill-conditioned and the algorithm unstable.
#'
#' `MARSSkfas()` uses the (Fortran-based) Kalman filter and smoother function ([KFAS::KFS()]) in the [KFAS](https://cran.r-project.org/package=KFAS) package (Helske 2012) based on the algorithms of Koopman and Durbin (2000, 2001, 2003).  The Koopman and Durbin algorithm is faster and more stable since it avoids matrix inverses.  Exact diffuse priors are also allowed in the KFAS Kalman filter function.  The standard output from the KFAS functions do not include the lag-one covariance smoother needed for the EM algorithm.  `MARSSkfas` computes the smoothed lag-one covariance using the Kalman filter applied to a stacked MARSS model as described on page 321 in Shumway and Stoffer (2000). Also the KFAS model specification only has the initial state at \eqn{t=1} (as \eqn{\mathbf{X}_1}{X(1)} conditioned on \eqn{\mathbf{y}_0}{y(0)}, which is missing).  When the initial state is specified at \eqn{t=0} (as \eqn{\mathbf{X}_0}{X(0)} conditioned on \eqn{\mathbf{y}_0}{y(0)}), `MARSSkfas()` computes the required \eqn{\textrm{E}[\mathbf{X}_1|\mathbf{y}_0}{E[X(1)|y(0)} and \eqn{\textrm{var}[\mathbf{X}_1|\mathbf{y}_0}{var[X(1)|y(0)} using the Kalman filter equations per Ghahramani and Hinton (1996).
#'
#' The likelihood returned for both functions is the exact likelihood when there are missing values rather than the approximate likelihood sometimes presented in texts for the missing values case.  The functions return the same filter, smoother and log-likelihood values.  The differences are that `MARSSkfas()` is faster (and more stable) but `MARSSkfss()` has many internal checks and error messages which can help debug numerical problems (but slow things down).  Also `MARSSkfss()` returns some output specific to the traditional filter algorithm (`J` and `Kt`).
#'
#' @return
#' A list with the following components. \eqn{m} is the number of state processes and \eqn{n} is the number of observation time series. "V" elements are called "P" in Shumway and Stoffer (2006, eqn 6.17 with s=T).  The output is referenced against equations in Shumway and Stoffer (2006) denoted S&S; the Kalman filter and smoother implemented in MARSS is for a more general MARSS model than that shown in S&S but the output has the same meaning.  In the expectations below, the parameters are left off; \eqn{\textrm{E}[\mathbf{X} | \mathbf{y}_1^t]}{E[X | y(1:t)]} is really \eqn{\textrm{E}[\mathbf{X} | \Theta, \mathbf{Y}_1^t=\mathbf{y}_1^t]}{E[X | Theta, Y(1:t)=y(1:t)]} where \eqn{\Theta}{Theta} is the parameter list. \eqn{\mathbf{y}_1^t}{y(1:t)} denotes the data from \eqn{t=1} to \eqn{t=t}.
#'
#' The notation for the conditional expectations is \eqn{\mathbf{x}_t^t}{xtt(t)} = \eqn{\textrm{E}[\mathbf{X} | \mathbf{y}_1^t]}{E[X | y(1:t)]}, \eqn{\mathbf{x}_t^{t-1}}{xtt1(t)} = \eqn{\textrm{E}[\mathbf{X} | \mathbf{y}_1^{t-1}]}{E[X | y(1:t-1)]} and \eqn{\mathbf{x}_t^T}{xtT(t)} = \eqn{\textrm{E}[\mathbf{X} | \mathbf{y}_1^T]}{E[X | y(1:T)]}. The conditional variances and covariances use similar notation. Note that in the Holmes (2012), the EM Derivation, \eqn{\mathbf{x}_t^T}{xtT(t)} and \eqn{\mathbf{V}_t^T}{VtT(t)} are given special symbols because they appear repeatedly: \eqn{\tilde{\mathbf{x}}_t}{tildex(t)} and \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} but here the more general notation is used.
#'
#' * xtT: \eqn{\mathbf{x}_t^T}{xtT(t)} State first moment conditioned on \eqn{\mathbf{y}_1^T}{y(1:T)}: \eqn{\textrm{E}[\mathbf{X}_t|\mathbf{y}_1^T]}{E[X(t) | y(1:T)]} (m x T matrix). Kalman smoother output.
#' * VtT: \eqn{\mathbf{V}_t^T}{VtT(t)} State variance matrix conditioned on \eqn{\mathbf{y}_1^T}{y(1:T)}: \eqn{\textrm{E}[(\mathbf{X}_t-\mathbf{x}_t^T)(\mathbf{X}_t-\mathbf{x}_t^T)^\top|\mathbf{y}_1^T]}{E[(X(t)-xtT(t))(x(t)-xtT(t))'| | y(1:T)]} (m x m x T array). Kalman smoother output. Denoted \eqn{P_t^T}{P_t^T} in S&S (S&S eqn 6.18 with \eqn{s=T}, \eqn{t1=t2=t}). The state second moment \eqn{\textrm{E}[\mathbf{X}_t\mathbf{X}_t^\top|\mathbf{y}_1^T]}{E[X(t)X(t)'| y(1:T)]} is equal to \eqn{\mathbf{V}_t^T + \mathbf{x}_t^T(\mathbf{x}_t^T)^\top}{VtT(t)+xtT(t)xtT(t)'}.
#' * Vtt1T: \eqn{\mathbf{V}_{t,t-1}^T}{Vtt1T(t)} State lag-one cross-covariance matrix \eqn{\textrm{E}[(\mathbf{X}_t-\mathbf{x}_t^T)(\mathbf{X}_{t-1}-\mathbf{x}_{t-1}^T)^\top|\mathbf{y}_1^T]}{E[(X(t)-xtT(t))(X(t-1)-xtT(t-1))' | y(1:T)]} (m x m x T). Kalman smoother output. \eqn{P_{t,t-1}^T} in S&S (S&S eqn 6.18 with \eqn{s=T}, \eqn{t1=t}, \eqn{t2=t-1}). State lag-one second moments \eqn{\textrm{E}[\mathbf{X}_t\mathbf{X}_{t-1}^\top|\mathbf{y}_1^T]}{E[X(t)X(t-1)'| y(1:T)]} is equal to \eqn{\mathbf{V}_{t, t-1}^T + \mathbf{x}_t^T(\mathbf{x}_{t-1}^T)^\top}{Vtt1T(t)+xtT(t)xtT(t-1)'}.
#' * x0T: Initial smoothed state estimate \eqn{\textrm{E}[\mathbf{X}_{t0}|\mathbf{y}_1^T}{E[X(t0) | y(1:T)]} (m x 1). If `model$tinitx=0`, \eqn{t0=0}; if `model$tinitx=1`, \eqn{t0=1}. Kalman smoother output.
#' * x01T: Smoothed state estimate \eqn{\textrm{E}[\mathbf{X}_1|\mathbf{y}_1^T}{E[X(1) | y(1:T)]} (m x 1).
#' * x00T: Smoothed state estimate \eqn{\textrm{E}[\mathbf{X}_0 |\mathbf{y}_1^T}{E[X(0) | y(1:T)]} (m x 1). If `model$tinitx=1`, this will be NA.
#' * V0T: Initial smoothed state covariance matrix \eqn{\textrm{E}[\mathbf{X}_{t0}\mathbf{X}_0^\top | \mathbf{y}_1^T}{E[X(t0)X(0)' | y(1:T)]} (m x m). If `model$tinitx=0`, \eqn{t0=0} and `V0T=V00T`; if `model$tinitx=1`, \eqn{t0=1} and `V0T=V10T`.  In the case of `tinitx=0`, this is \eqn{P_0^T} in S&S.
#' * V10T: Smoothed state covariance matrix \eqn{\textrm{E}[\mathbf{X}_1\mathbf{X}_0^\top | \mathbf{y}_1^T}{E[X(1)X(0)' | y(1:T)]} (m x m).
#' * V00T: Smoothed state covariance matrix \eqn{\textrm{E}[\mathbf{X}_0\mathbf{X}_0^\top | \mathbf{y}_1^T}{E[X(0)X(0)' | y(1:T)]} (m x m). If `model$tinitx=1`, this will be NA.
#' * J: (m x m x T) Kalman smoother output.  Only for `MARSSkfss()`. (ref S&S eqn 6.49)
#' * J0: J at the initial time (t=0 or t=1) (m x m x T). Kalman smoother output. Only for `MARSSkfss()`.
#' * xtt: State first moment conditioned on \eqn{\mathbf{y}_1^t}{y(1:t)}: \eqn{\textrm{E}[\mathbf{X}_t | \mathbf{y}_1^t}{E[X(t) | y(1:t)]} (m x T).  Kalman filter output. (S&S eqn 6.17 with \eqn{s=t})
#' * xtt1: State first moment conditioned on \eqn{\mathbf{y}_1^{t-1}}{y(1:t-1)}: \eqn{\textrm{E}[\mathbf{X}_t | \mathbf{y}_1^{t-1}}{E[X(t) | y(1:t-1)]} (m x T).  Kalman filter output. (S&S eqn 6.17 with \eqn{s=t-1})
#' * Vtt: State variance conditioned on \eqn{\mathbf{y}_1^t}{y(1:t)}: \eqn{\textrm{E}[(\mathbf{X}_t-\mathbf{x}_t^t)(\mathbf{X}_t-\mathbf{x}_t^t)^\top|\mathbf{y}_1^t]}{E[(X(t)-xtt(t))(X(t)-xtt(t))'| | y(1:t)]} (m x m x T array). Kalman filter output. \eqn{P_t^t} in S&S (S&S eqn 6.18 with s=t, t1=t2=t). The state second moment \eqn{\textrm{E}[\mathbf{X}_t\mathbf{X}_t^\top|\mathbf{y}_1^t]}{E[X(t)X(t)'| y(1:t)]} is equal to \eqn{\mathbf{V}_t^t + \mathbf{x}_t^t(\mathbf{x}_t^t)^\top}{Vtt(t)+xtt(t)xtt(t)'}.
#' * Vtt1: State variance conditioned on \eqn{\mathbf{y}_1^{t-1}}{y(1:t-1)}: \eqn{\textrm{E}[(\mathbf{X}_t-\mathbf{x}_t^{t-1})(\mathbf{X}_t-\mathbf{x}_t^{t-1})^\top|\mathbf{y}_1^{t-1}]}{E[(X(t)-xtt1(t))(X(t)-xtt1(t))'| | y(1:t-1)]} (m x m x T array). Kalman filter output. The state second moment \eqn{\textrm{E}[\mathbf{X}_t\mathbf{X}_t^\top|\mathbf{y}_1^{t-1}]}{E[X(t)X(t)'| y(1:t-1)]} is equal to \eqn{\mathbf{V}_t^{t-1} + \mathbf{x}_t^{t-1}(\mathbf{x}_t^{t-1})^\top}{Vtt1(t)+xtt1(t)xtt1(t)'}.
#' * Kt: Kalman gain (m x m x T). Kalman filter output (ref S&S eqn 6.23). Only for `MARSSkfss()`.
#' * Innov: Innovations \eqn{\mathbf{y}_t-\textrm{E}[\mathbf{Y}_t|\mathbf{y}_1^{t-1}]}{y(t) - E[Y(t) | y(1:t-1)]} (n x T). Kalman filter output. Only returned with `MARSSkfss()`. (ref page S&S 339).
#' * Sigma: Innovations covariance matrix. Kalman filter output. Only returned with `MARSSkfss()`. (ref S&S eqn 6.61)
#' * logLik: Log-likelihood logL(y(1:T) | Theta) (ref S&S eqn 6.62)
#' * kfas.model: The model in [KFAS::KFAS()] model form (class [KFAS::SSModel()]). Only for `MARSSkfas`.
#' * errors: Any error messages.
#'
#' @references
#' A. C. Harvey (1989).  Chapter 5, Forecasting, structural time series models and the Kalman filter.  Cambridge University Press.
#'
#' R. H. Shumway and D. S. Stoffer (2006).  Time series analysis and its applications: with R examples.  Second Edition. Springer-Verlag, New York.
#'
#' Ghahramani, Z. and Hinton, G.E. (1996) Parameter estimation for linear dynamical systems.  University of Toronto Technical Report CRG-TR-96-2.
#'
#' Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME] `RShowDoc("EMDerivation",package="MARSS")` to open a copy.
#'
#' Jouni Helske (2012). KFAS: Kalman filter and smoother for exponential family state space models. <https://CRAN.R-project.org/package=KFAS>
#'
#' Koopman, S.J. and Durbin J. (2000). Fast filtering and smoothing for non-stationary time series models, Journal of American Statistical Association, 92, 1630-38.
#'
#' Koopman, S.J. and Durbin J. (2001). Time series analysis by state space methods. Oxford: Oxford University Press.
#'
#' Koopman, S.J. and Durbin J. (2003). Filtering and smoothing of state vector for diffuse state space models, Journal of Time Series Analysis, Vol. 24, No. 1.
#'
#' The MARSS User Guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type `RShowDoc("UserGuide",package="MARSS")` to open a copy.
#'
#' @author Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSS()], [marssMODEL()], [MARSSkem()], [KFAS::KFAS()]
#'
#' @aliases MARSSkfas MARSSkfss
#'
#' @examples
#' dat <- t(harborSeal)
#' dat <- dat[2:nrow(dat), ]
#' # you can use MARSS to construct a marssMLE object
#' # MARSS calls MARSSinits to construct default initial values
#' # with fit = FALSE, the $par element of the marssMLE object will be NULL
#' fit <- MARSS(dat, fit = FALSE)
#' # MARSSkf needs a marssMLE object with the par element set
#' fit$par <- fit$start
#' # Compute the kf output at the params used for the inits
#' kfList <- MARSSkf(fit)
#'
#' @export
MARSSkf <- function(MLEobj, only.logLik = FALSE, return.lag.one = TRUE, return.kfas.model = FALSE, newdata = NULL, smoother = TRUE) {
  if (is.null(MLEobj$par)) {
    stop("Stopped in MARSSkf(): par element of marssMLE object is required.\n")
  }
  if (!missing(newdata)) {
    data.dim <- attr(MLEobj[["model"]], "model.dims")[["data"]]
    if (!identical(dim(newdata), data.dim)) {
      stop("Stopped in MARSSkf(): newdata must be the same dimensions as the data used to fit the model.\n")
    }
    if (!is.numeric(newdata)) {
      stop("Stopped in MARSSkf(): newdata must be numeric.\n")
    }
    if (dim(MLEobj[["model"]][["free"]][["x0"]])[2] != 0) {
      message("MARSSkf(): x0 was estimated and this x0 will be used with newdata.\n")
    }
    MLEobj[["model"]][["data"]] <- newdata
    MLEobj[["marss"]][["data"]] <- newdata
  }
  if (MLEobj$fun.kf == "MARSSkfss") {
    return(MARSSkfss(MLEobj, smoother = smoother))
  }
  if (MLEobj$fun.kf == "MARSSkfas") {
    return(MARSSkfas(MLEobj, only.logLik = only.logLik, return.lag.one = return.lag.one, return.kfas.model = return.kfas.model))
  }
  return(list(ok = FALSE, errors = "kf.function does not specify a valid Kalman filter and smoother function."))
}
