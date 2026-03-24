#' Fit a MARSS Model via Maximum-Likelihood Estimation
#'
#' @description
#' This is the main function for fitting multivariate autoregressive
#' state-space (MARSS) models with linear constraints. Scroll down to the
#' bottom to see some short examples. To open a guide to show you how to get
#' started quickly, type `RShowDoc("Quick_Start",package="MARSS")`. To open
#' the MARSS User Guide from the command line, type
#' `RShowDoc("UserGuide",package="MARSS")`. To get an overview of the package
#' and all its main functions and how to get output (parameter estimates,
#' fitted values, residuals, Kalman filter or smoother output, or plots), go
#' to [MARSS-package]. If `MARSS()` is throwing errors or warnings that you
#' don't understand, try the Troubleshooting section of the user guide or type
#' [MARSSinfo()] at the command line.
#'
#' The default MARSS model form is "marxss", which is Multivariate
#' Auto-Regressive(1) eXogenous inputs State-Space model:
#' \deqn{\mathbf{x}_{t} = \mathbf{B}_t \mathbf{x}_{t-1} + \mathbf{u}_t + \mathbf{C}_t \mathbf{c}_t + \mathbf{G}_t \mathbf{w}_t, \textrm{ where } \mathbf{W}_t \sim \textrm{MVN}(0,\mathbf{Q}_t)}{x(t) = B(t) x(t-1) + u(t) + C(t) c(t) + G(t) w(t), where W(t) ~ MVN(0,Q(t))}
#' \deqn{\mathbf{y}_t = \mathbf{Z}_t \mathbf{x}_t + \mathbf{a}_t + \mathbf{D}_t \mathbf{d}_t + \mathbf{H}_t \mathbf{v}_t, \textrm{ where } \mathbf{V}_t \sim \textrm{MVN}(0,\mathbf{R}_t)}{y(t) = Z(t) x(t) + a(t) + D(t) d(t) + H(t) v(t), where V(t) ~ MVN(0,R(t))}
#' \deqn{\mathbf{X}_1 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) \textrm{ or } \mathbf{X}_0 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) }{X(1) ~ MVN(x0, V0) or X(0) ~  MVN(x0, V0) }
#' The parameters are everything except \eqn{\mathbf{x}}{x}, \eqn{\mathbf{y}}{y},
#' \eqn{\mathbf{v}}{v}, \eqn{\mathbf{w}}{w}, \eqn{\mathbf{c}}{c} and
#' \eqn{\mathbf{d}}{d}. \eqn{\mathbf{y}}{y} are data (missing values allowed).
#' \eqn{\mathbf{c}}{c} and \eqn{\mathbf{d}}{d} are inputs (no missing values
#' allowed). All parameters (except \eqn{\mathbf{x0}}{x0} and
#' \eqn{\mathbf{V0}}{V0}) can be time-varying but by default, all are
#' time-constant (and the MARSS equation is generally written without the
#' \eqn{t} subscripts on the parameter matrices). All parameters can be zero,
#' including the variance matrices.
#'
#' The parameter matrices can have fixed values and linear constraints. This
#' is an example of a 3x3 matrix with linear constraints. All matrix elements
#' can be written as a linear function of \eqn{a}, \eqn{b}, and \eqn{c}:
#' \deqn{\left[\begin{array}{c c c} a+2b & 1 & a\\ 1+3a+b & 0 & b \\ 0 & -2 & c\end{array}\right]}{[a+2b  1   a \n 1+3a+b   0   b \n 0   -2    c ]}
#'
#' Values such as \eqn{a b} or \eqn{a^2} or \eqn{log(a)} are not linear constraints.
#'
#' @param y A n x T matrix of n time series over T time steps. Only y is
#'   required for the function. A ts object (univariate or multivariate) can
#'   be used and this will be converted to a matrix with time in the columns.
#' @param inits A list with the same form as the list outputted by `coef(fit)`
#'   that specifies initial values for the parameters. See also [MARSS.marxss()].
#' @param model Model specification using a list of parameter matrix text
#'   shortcuts or matrices. See Details and [MARSS.marxss()] for the default
#'   form. Or better yet open the Quick Start Guide
#'   `RShowDoc("Quick_Start",package="MARSS")`.
#' @param miss.value Deprecated. Denote missing values by NAs in your data.
#' @param method Estimation method. MARSS provides an EM algorithm
#'   (`method="kem"`) (see [MARSSkem()]) and the BFGS algorithm
#'   (`method="BFGS"`) (see [MARSSoptim()]).
#' @param form The equation form used in the `MARSS()` call. The default is
#'   "marxss". See [MARSS.marxss()] or [MARSS.dfa()].
#' @param fit TRUE/FALSE Whether to fit the model to the data. If FALSE, a
#'   [marssMLE] object with only the model is returned.
#' @param silent Setting to TRUE(1) suppresses printing of full error messages,
#'   warnings, progress bars and convergence information. Setting to FALSE(0)
#'   produces error output. Setting silent=2 will produce more verbose error
#'   messages and progress information.
#' @param fun.kf What Kalman filter function to use. MARSS has two:
#'   [MARSSkfas()] which is based on the Kalman filter in the
#'   [KFAS](https://cran.r-project.org/package=KFAS) package based on Koopman
#'   and Durbin and [MARSSkfss()] which is a native R implementation of the
#'   Kalman filter and smoother in Shumway and Stoffer. The KFAS filter is
#'   much faster. [MARSSkfas()] modifies the input and output in order to
#'   output the lag-one covariance smoother needed for the EM algorithm (per
#'   page 321 in Shumway and Stoffer (2000).
#' @param control Estimation options for the maximization algorithm. The
#'   typically used control options for method="kem" are below but see
#'   [marssMLE] for the full list of control options. Note many of these are
#'   not allowed if method="BFGS"; see [MARSSoptim()] for the allowed control
#'   options for this method.
#'
#'   * `minit` The minimum number of iterations to do in the maximization
#'     routine (if needed by method). If `method="kem"`, this is an easy way
#'     to up the iterations and see how your estimates are converging.
#'     (positive integer)
#'   * `maxit` Maximum number of iterations to be used in the maximization
#'     routine (if needed by method) (positive integer).
#'   * `min.iter.conv.test` Minimum iterations to run before testing
#'     convergence via the slope of the log parameter versus log iterations.
#'   * `conv.test.deltaT=9` Number of iterations to use for the testing
#'     convergence via the slope of the log parameter versus log iterations.
#'   * `conv.test.slope.tol` The slope of the log parameter versus log
#'     iteration to use as the cut-off for convergence. The default is 0.5
#'     which is a bit high. For final analyses, this should be set lower. If
#'     you want to only use abstol as your convergence test, then set to
#'     something very large, for example `conv.test.slope.tol=1000`. Type
#'     `MARSSinfo(11)` to see some comments on when you might want to do this.
#'   * `abstol` The logLik.(iter-1)-logLik.(iter) convergence tolerance for
#'     the maximization routine. To meet convergence both the abstol and slope
#'     tests must be passed.
#'   * `allow.degen` Whether to try setting \eqn{\mathbf{Q}}{Q} or
#'     \eqn{\mathbf{R}}{R} elements to zero if they appear to be going to zero.
#'   * `trace` An integer specifying the level of information recorded and
#'     error-checking run during the algorithms. `trace=0` specifies basic
#'     error-checking and brief error-messages; `trace>0` will print full
#'     error messages. In addition if trace>0, the Kalman filter output will
#'     be added to the outputted [marssMLE] object. Additional information
#'     recorded depends on the method of maximization. For the EM algorithm, a
#'     record of each parameter estimate for each EM iteration will be added.
#'     See [optim()] for trace output details for the BFGS method.
#'     `trace=-1` will turn off most internal error-checking and most error
#'     messages. The internal error checks are time expensive so this can speed
#'     up model fitting. This is particularly useful for bootstrapping and
#'     simulation studies. It is also useful if you get an error saying that
#'     `MARSS()` stops in [MARSSkfss()] due to a `chol()` call. `MARSSkfss()`
#'     uses matrix inversions and for some models these are unstable (high
#'     condition value). `MARSSkfss()` is used for error-checks and does not
#'     need to be called normally.
#'   * `safe` Setting `safe=TRUE` runs the Kalman smoother after each
#'     parameter update rather than running the smoother only once after
#'     updating all parameters. The latter is faster but is not a strictly
#'     correct EM algorithm. In most cases, `safe=FALSE` (default) will not
#'     change the fits. If this setting does cause problems, you will know
#'     because you will see an error regarding the log-likelihood dropping and
#'     it will direct you to set `safe=TRUE`.
#' @param ... Optional arguments passed to function specified by form.
#'
#' @details
#' The `model` argument specifies the structure of your model. There is a
#' one-to-one correspondence between how you would write your model in matrix
#' form on the whiteboard and how you specify the model for `MARSS()`. Many
#' different types of multivariate time-series models can be converted to the
#' MARSS form. See the
#' [User Guide](https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf)
#' and [Quick Start Guide](https://cran.r-project.org/package=MARSS/vignettes/Quick_Start.html)
#' for examples.
#'
#' The MARSS package has two forms for standard users: marxss and dfa.
#'
#' * [MARSS.marxss()] This is the default form. This is a MARSS model with
#'   (optional) inputs \eqn{\mathbf{c}_t}{c(t)} or \eqn{\mathbf{d}_t}{d(t)}.
#'   Most users will want this help page.
#' * [MARSS.dfa()] This is a model form to allow easier specification of
#'   models for Dynamic Factor Analysis. The \eqn{\mathbf{Z}}{Z} parameters
#'   has a specific form and the \eqn{\mathbf{Q}}{Q} is set at i.i.d
#'   (diagonal) with variance of 1.
#'
#' Those looking to modify or understand the base code, should look at
#' [MARSS.marss()] and `MARSS.vectorized()`. These describe the forms used by
#' the base functions. The EM algorithm uses the MARSS model written in
#' vectorized form. This form is what allows linear constraints.
#'
#' The likelihood surface for MARSS models can be multimodal or with strong
#' ridges. It is recommended that for final analyses the estimates are checked
#' by using a Monte Carlo initial conditions search; see the chapter on initial
#' conditions searches in the User Guide. This requires more computation time,
#' but reduces the chance of the algorithm terminating at a local maximum and
#' not reaching the true MLEs. Also it is wise to check the EM results against
#' the BFGS results (if possible) if there are strong ridges in the likelihood.
#' Such ridges seems to slow down the EM algorithm considerably and can cause
#' the algorithm to report convergence far from the maximum-likelihood values.
#' EM steps up the likelihood and the convergence test is based on the rate of
#' change of the log-likelihood in each step. Once on a strong ridge, the
#' steps can slow dramatically. You can force the algorithm to keep working by
#' setting `minit`. BFGS seems less hindered by the ridges but can be
#' prodigiously slow for some multivariate problems. BFGS tends to work better
#' if you give it good initial conditions (see Examples below for how to do
#' this).
#'
#' If you are working with models with time-varying parameters, it is important
#' to notice the time-index for the parameters in the process equation (the
#' \eqn{\mathbf{x}}{x} equation). In some formulations (e.g. in
#' [KFAS::KFAS]), the process equation is
#' \eqn{\mathbf{x}_t=\mathbf{B}_{t-1}\mathbf{x}_{t-1}+\mathbf{w}_{t-1}}{x(t)=B(t-1)x(t-1)+w(t-1)}
#' so \eqn{\mathbf{B}_{t-1}}{B(t-1)} goes with \eqn{\mathbf{x}_t}{x(t)} not
#' \eqn{\mathbf{B}_t}{B(t)}. Thus one needs to be careful to line up the time
#' indices when passing in time-varying parameters to `MARSS()`. See the User
#' Guide for examples.
#'
#' @return
#' An object of class [marssMLE]. The structure of this object is discussed
#' below, but if you want to know how to get specific output (like residuals,
#' coefficients, smoothed states, confidence intervals, etc), see
#' [print.marssMLE()], [tidy.marssMLE()], [MARSSresiduals()] and
#' [plot.marssMLE()].
#'
#' The outputted [marssMLE] object has the following components:
#'
#' * `model` MARSS model specification. It is a [marssMODEL] object in the
#'   form specified by the user in the `MARSS()` call. This is used by print
#'   functions so that the user sees the expected form.
#' * `marss` The [marssMODEL] object in marss form. This form is needed for
#'   all the internal algorithms, thus is a required part of a [marssMLE]
#'   object.
#' * `call` All the information passed in in the `MARSS()` call.
#' * `start` List with specifying initial values that were used for each
#'   parameter matrix.
#' * `control` A list of estimation options, as specified by arguments
#'   `control`.
#' * `method` Estimation method.
#'
#' If `fit=TRUE`, the following are also added to the [marssMLE] object.
#' If `fit=FALSE`, a [marssMLE] object ready for fitting via the specified
#' `method` is returned.
#'
#' * `par` A list of estimated parameter values in marss form. Use
#'   [print.marssMLE()], [tidy.marssMLE()] or [coef.marssMLE()] for outputting
#'   the model estimates in the `MARSS()` call (e.g. the default "marxss" form).
#' * `states` The expected value of \eqn{\mathbf{X}}{X} conditioned on all
#'   the data, i.e. smoothed states.
#' * `states.se` The standard errors of the expected value of
#'   \eqn{\mathbf{X}}{X}.
#' * `ytT` The expected value of \eqn{\mathbf{Y}}{Y} conditioned on all the
#'   data. Note this is just \eqn{y} for those \eqn{y} that are not missing.
#' * `ytT.se` The standard errors of the expected value of
#'   \eqn{\mathbf{Y}}{Y}. Note this is 0 for any non-missing \eqn{y}.
#' * `numIter` Number of iterations required for convergence.
#' * `convergence` Convergence status. 0 means converged successfully, 3
#'   means all parameters were fixed (so model did not need to be fit) and -1
#'   means call was made with `fit=FALSE` and parameters were not fixed (thus
#'   no `$par` element and Kalman filter/smoother cannot be run). Anything
#'   else is a warning or error. 2 means the [marssMLE] object has an error;
#'   the object is returned so you can debug it. The other numbers are errors
#'   during fitting. The error code depends on the fitting method. See
#'   [MARSSkem()] and [MARSSoptim()].
#' * `logLik` Log-likelihood.
#' * `AIC` Akaike's Information Criterion.
#' * `AICc` Sample size corrected AIC.
#'
#' If `control$trace` is set to 1 or greater, the following are also added to
#' the [marssMLE] object.
#'
#' * `kf` A list containing Kalman filter/smoother output from [MARSSkf()].
#'   This is not normally added to a [marssMLE] object since it is verbose,
#'   but can be added using [MARSSkf()].
#' * `Ey` A list containing output from [MARSShatyt()]. This isn't normally
#'   added to a [marssMLE] object since it is verbose, but can be computed
#'   using [MARSShatyt()].
#'
#' @references
#' The MARSS User Guide: Holmes, E. E., E. J. Ward, and M. D. Scheuerell
#' (2012) Analysis of multivariate time-series using the MARSS package. NOAA
#' Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E.,
#' Seattle, WA 98112. Type `RShowDoc("UserGuide",package="MARSS")` to open a
#' copy.
#'
#' Holmes, E. E. (2012). Derivation of the EM algorithm for constrained and
#' unconstrained multivariate autoregressive state-space (MARSS) models.
#' Technical Report. arXiv:1302.3919 [stat.ME]
#'
#' Holmes, E. E., E. J. Ward and K. Wills. (2012) MARSS: Multivariate
#' autoregressive state-space models for analyzing time-series data. R Journal
#' 4: 11-19.
#'
#' @author
#' Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.
#'
#' @seealso
#' [marssMLE], [MARSSkem()], [MARSSoptim()], [MARSSkf()], [MARSS-package],
#' [print.marssMLE()], [plot.marssMLE()], [print.marssMODEL()],
#' [MARSS.marxss()], [MARSS.dfa()], [fitted.marssMLE()],
#' [residuals.marssMLE()], [MARSSresiduals()], [predict.marssMLE()],
#' [tsSmooth.marssMLE()], [tidy.marssMLE()], [coef.marssMLE()]
#'
#' @examples
#' dat <- t(harborSealWA)
#' dat <- dat[2:4, ] # remove the year row
#' # fit a model with 1 hidden state and 3 observation time series
#' kemfit <- MARSS(dat, model = list(
#'   Z = matrix(1, 3, 1),
#'   R = "diagonal and equal"
#' ))
#' kemfit$model # This gives a description of the model
#' print(kemfit$model) # same as kemfit$model
#' summary(kemfit$model) # This shows the model structure
#'
#' # add CIs to a marssMLE object
#' # default uses an estimated Hessian matrix
#' kem.with.hess.CIs <- MARSSparamCIs(kemfit)
#' kem.with.hess.CIs
#'
#' # fit a model with 3 hidden states (default)
#' kemfit <- MARSS(dat, silent = TRUE) # suppress printing
#' kemfit
#'
#' # Fit the above model with BFGS using a short EM fit as initial conditions
#' kemfit <- MARSS(dat, control=list(minit=5, maxit=5))
#' bffit <- MARSS(dat, method="BFGS", inits=kemfit)
#'
#' # fit a model with 3 correlated hidden states
#' # with one variance and one  covariance
#' # maxit set low to speed up example, but more iters are needed for convergence
#' kemfit <- MARSS(dat, model = list(Q = "equalvarcov"), control = list(maxit = 50))
#' # use Q="unconstrained" to allow different variances and covariances
#'
#' # fit a model with 3 independent hidden states
#' # where each observation time series is independent
#' # the hidden trajectories 2-3 share their U parameter
#' kemfit <- MARSS(dat, model = list(U = matrix(c("N", "S", "S"), 3, 1)))
#'
#' # same model, but with fixed independent observation errors
#' # and the 3rd x processes are forced to have a U=0
#' # Notice how a list matrix is used to combine fixed and estimated elements
#' # all parameters can be specified in this way using list matrices
#' kemfit <- MARSS(dat, model = list(U = matrix(list("N", "N", 0), 3, 1), R = diag(0.01, 3)))
#'
#' # fit a model with 2 hidden states (north and south)
#' # where observation time series 1-2 are north and 3 is south
#' # Make the hidden state process independent with same process var
#' # Make the observation errors different but independent
#' # Make the growth parameters (U) the same
#' # Create a Z matrix as a design matrix that assigns the "N" state to the first 2 rows of dat
#' # and the "S" state to the 3rd row of data
#' Z <- matrix(c(1, 1, 0, 0, 0, 1), 3, 2)
#' # You can use factor is a shortcut making the above design matrix for Z
#' # Z <- factor(c("N","N","S"))
#' # name the state vectors
#' colnames(Z) <- c("N", "S")
#' kemfit <- MARSS(dat, model = list(
#'   Z = Z,
#'   Q = "diagonal and equal", R = "diagonal and unequal", U = "equal"
#' ))
#'
#' # print the model followed by the marssMLE object
#' kemfit$model
#'
#' \dontrun{
#' # simulate some new data from our fitted model
#' sim.data <- MARSSsimulate(kemfit, nsim = 10, tSteps = 10)
#'
#' # Compute bootstrap AIC for the model; this takes a long, long time
#' kemfit.with.AICb <- MARSSaic(kemfit, output = "AICbp")
#' kemfit.with.AICb
#' }
#'
#' \dontrun{
#' # Many more short examples can be found in the
#' # Quick Examples chapter in the User Guide
#' RShowDoc("UserGuide", package = "MARSS")
#'
#' # You can find the R scripts from the chapters by
#' # going to the index page
#' RShowDoc("index", package = "MARSS")
#' }
#'
#' @export
MARSS <- function(y,
                  model = NULL,
                  inits = NULL,
                  miss.value = as.numeric(NA),
                  method = c("kem", "BFGS", "TMB", "BFGS_TMB", "nlminb_TMB"),
                  form = c("marxss", "dfa", "marss"),
                  fit = TRUE,
                  silent = FALSE,
                  control = NULL,
                  fun.kf = c("MARSSkfas", "MARSSkfss"),
                  ...) {
  # If user did not pass in fun.kf, then MARSSkfas will be used by default. MARSSkfss will be tried if that fails
  pkg <- "MARSS"
  if (missing(fun.kf)) missing.fun.kf <- FALSE else missing.fun.kf <- TRUE
  fun.kf <- match.arg(fun.kf)
  method <- match.arg(method)
  form <- match.arg(form)
  allowed.methods <- get("allowed.methods", envir = pkg_globals)
  # Some error checks depend on an allowable method
  if (length(grep("TMB", method)) > 0) {
    if(!requireNamespace("marssTMB")){
      message("Fitting with TMB requires the  'TMB' package. Please install marssTMB from https://atsa-es.github.io/marssTMB/")
      return(invisible())
    }
  }
  ## Start by checking the data, since if the data have major problems then the rest of the code
  ## will have problems
  if (!missing(miss.value)) {
    stop("miss.value is deprecated in MARSS.  Replace missing values in y with NA.\n")
  }
  if (is.null(y)) {
    stop("MARSS: No data (y) passed in.", call. = FALSE)
  }
  if (!(is.vector(y) | is.matrix(y) | inherits(y, "ts"))) stop("MARSS: Data (y) must be a vector, matrix (time going across columns) or ts/mts object.", call. = FALSE)
  if (length(y) == 0) stop("MARSS: Data (y) is length 0.", call. = FALSE)
  if (is.vector(y)) y <- matrix(y, nrow = 1)
  if (inherits(y, "ts")) {
    model.tsp <- stats::tsp(y)
    y <- t(y)
  } else {
    model.tsp <- c(1, ncol(y), 1)
  }
  attr(y, "model.tsp") <- model.tsp
  if (any(is.nan(y))) cat("MARSS: NaNs in data are being replaced with NAs.  There might be a problem if NaNs shouldn't be in the data.\nNA is the normal missing value designation.\n")
  y[is.na(y)] <- as.numeric(NA)

  if (inherits(model, "marssMLE")) model <- c(coef(model, type = "matrix"), tinitx = model$model$tinitx, diffuse = model$model$diffuse)
  if (inherits(model, "marssMODEL")) model <- c(marssMODEL.to.list(model), tinitx = model$tinitx, diffuse = model$diffuse)

  MARSS.call <- list(data = y, inits = inits, model = model, control = control, method = method, form = form, silent = silent, fit = fit, fun.kf = fun.kf, ...)

  # First make sure specified equation form has a corresponding function to do the conversion to marssMODEL (form=marss) object
  as.marss.fun <- paste("MARSS.", form[1], sep = "")
  tmp <- try(exists(as.marss.fun, mode = "function"), silent = TRUE)
  if (!isTRUE(tmp)) {
    msg <- paste(" MARSS.", form[1], "() function to construct a marssMODEL (form=marss) object does not exist.\n", sep = "")
    cat("\n", "Errors were caught in MARSS \n", msg, sep = "")
    stop("Stopped in MARSS() due to problem(s) with required arguments.\n", call. = FALSE)
  }

  # Build the marssMODEL object from the model argument to MARSS()
  ## The as.marss.fun() call adds the following to MARSS.inputs list:
  ## $marss Translate model strucuture names (shortcuts) into a marssMODEL (form=marss) object put in $marss
  ## $alt.forms with any alternate marssMODEL objects in other forms that might be needed later
  ## a marssMODEL object is a list(data, fixed, free, tinitx, diffuse)
  ## with attributes model.dims, X.names, form, equation
  ## error checking within the function is a good idea though not required
  ## if changes to the control values are wanted these can be set by changing MARSS.inputs$control
  if (silent == 2) cat("Building the marssMODEL object from the model argument to MARSS().\n")
  MARSS.inputs <- eval(call(as.marss.fun, MARSS.call))
  marss.object <- MARSS.inputs$marss
  if (silent == 2) cat(paste("Resulting model has ", attr(marss.object, "model.dims")$x[1], " state processes and ", attr(marss.object, "model.dims")$y[1], " observation processes.\n", sep = ""))

  ## Check that the marssMODEL object output by MARSS.form() is ok
  ## More checking on the control list is done by is.marssMLE() to make sure the MLEobj is ready for fitting
  if (!identical(MARSS.call$control$trace, -1)) { # turn off all error checking
    if (silent == 2) cat(paste("Checking that the marssMODEL object output by MARSS.", form, "() is ok.\n", sep = ""))
    tmp <- is.marssMODEL(marss.object, method = method)
    if (!isTRUE(tmp)) {
      if (!silent || silent == 2) {
        cat(tmp)
        return(marss.object)
      }
      stop("Stopped in MARSS() due to problem(s) with model specification. If silent=FALSE, marssMODEL object (in marss form) will be returned.\n", call. = FALSE)
    } else {
      if (silent == 2) cat("  marssMODEL is ok.\n")
    }
  }

  ## checkMARSSInputs() call does the following:
  ## Check that the user didn't pass in any illegal arguments
  ## and fill in defaults if some params left off
  ## This does not check model since the marssMODEL object is constructed
  ## via the MARSS.form() function above
  if (silent == 2) cat("Running checkMARSSInputs().\n")
  MARSS.inputs <- checkMARSSInputs(MARSS.inputs, silent = FALSE)


  ###########################################################################################################
  ##  MODEL FITTING
  ###########################################################################################################

  ## MLE estimation
  kem.methods <- get("kem.methods", envir = pkg_globals)
  optim.methods <- get("optim.methods", envir = pkg_globals)
  if (method %in% allowed.methods) {
    ## Create the marssMLE object

    MLEobj <- list(marss = marss.object, model = MARSS.inputs$model, control = c(MARSS.inputs$control, silent = silent), method = method, fun.kf = fun.kf)
    # Set the call form since that info needed for MARSSinits
    if (MLEobj$control$trace != -1) {
      MLEobj$call <- MARSS.call
    }

    # This is a helper function to set simple inits for a marss MLE model object
    if (silent == 2) cat("Running MARSSinits() to set start conditions.\n")
    MLEobj$start <- MARSSinits(MLEobj, MARSS.inputs$inits)

    class(MLEobj) <- c("marssMLE", method)

    ## Check the marssMLE object
    ## is.marssMLE() calls is.marssMODEL() to check the model,
    ## then checks dimensions of initial value matrices.
    ## it also checks the control list and add defaults if some values are NULL
    if (MLEobj$control$trace != -1) {
      if (silent == 2) cat("Checking the marssMLE object before fitting.\n")
      tmp <- is.marssMLE(MLEobj)
    }

    # if errors, tmp will not be true, it will be error messages
    if (!isTRUE(tmp)) {
      if (!silent || silent == 2) {
        cat(tmp)
        cat(" The incomplete/inconsistent MLE object is being returned.\n")
      }
      cat("Error: Stopped in MARSS() due to marssMLE object incomplete or inconsistent. \nPass in silent=FALSE to see the errors.\n\n")
      MLEobj$convergence <- 2
      return(MLEobj)
    }

    # MLEobj is ok.
    if (!fit) {
      # will be set to 3 if all fixed
      MLEobj$convergence <- -1
    }

    # If all parameters fixed. Set convergence=3 whether or not fit=TRUE
    model.is.fixed <- FALSE
    if (all(unlist(lapply(MLEobj[["marss"]][["free"]], is.fixed)))) {
      model.is.fixed <- TRUE
      if (silent == 2) cat("All parameters fixed. No estimation done.\n")
      MLEobj$convergence <- 3
      MLEobj$par <- list()
      for (el in attr(MLEobj[["marss"]], "par.names")) MLEobj[["par"]][[el]] <- matrix(0, 0, 1)
      kf.out <- try(MARSSkf(MLEobj, only.logLik = TRUE), silent = TRUE)
      if (inherits(kf.out, "try-error")) {
        MLEobj$convergence <- 53
        MLEobj$logLik <- NA
      } else {
        # rest of the output will be added below
        MLEobj$logLik <- kf.out$logLik
        MLEobj <- MARSSaic(MLEobj)
        MLEobj$coef <- coef(MLEobj, type = "vector")
        kf.out <- try(MARSSkf(MLEobj), silent = TRUE)
        if (inherits(kf.out, "try-error")) {
          MLEobj$convergence <- 54
        }
      }
    }

    # If fitting is needed, check that Kalman filter/smoother will run
    if (fit && MLEobj$control$trace != -1 && !model.is.fixed) {
      MLEobj.test <- MLEobj
      MLEobj.test$par <- MLEobj$start
      kftest <- try(MARSSkf(MLEobj.test), silent = TRUE)
      if (inherits(kftest, "try-error")) {
        cat("Error: Stopped in MARSS() before fitting because", fun.kf, "stopped.  Something in the model structure prevents the Kalman filter/smoother (KF) from running.\n Try setting fun.kf to use a different KF function (MARSSkfss or MARSSkfas) or use fit=FALSE and check the model you are trying to fit. You can also try trace=1 to get more progress output. You could try trace=-1 to bypass the initial KF check if you are using method='BFGS' and know the logLik function will run.\n\n", kftest$condition)
        MLEobj.test$convergence <- 2
        return(MLEobj.test)
      }
      if (!kftest$ok) {
        cat(kftest$msg)
        cat(paste("Error: Stopped in MARSS() before fitting because", fun.kf, "stopped.  Something in the model structure prevents the Kalman filter or smoother running.\n Try setting fun.kf to use a different KF function (MARSSkfss or MARSSkfas) or use fit=FALSE and check the model you are trying to fit. You can also try trace=1 to get more progress output.\n", kftest$errors, "\n", sep = ""))
        MLEobj.test$convergence <- 2
        return(MLEobj.test)
      }
      MLEobj.test$kf <- kftest
    }
    # Ey is needed for method=kem
    if (fit && !model.is.fixed && MLEobj$control$trace != -1 && MLEobj$method %in% kem.methods) {
      Eytest <- try(MARSShatyt(MLEobj.test), silent = TRUE)
      if (inherits(Eytest, "try-error")) {
        cat("Error: Stopped in MARSS() before fitting because MARSShatyt() stopped.  Something is wrong with the model structure that prevents MARSShatyt() running.\n\n", Eytest$condition)
        MLEobj.test$convergence <- 2
        MLEobj.test$Ey <- Eytest
        return(MLEobj.test)
      }
      if (!Eytest$ok) {
        cat(Eytest$msg)
        cat("Error: Stopped in MARSS() before fitting because MARSShatyt returned errors.  Something is wrong with the model structure that prevents function running.\n\n")
        MLEobj.test$convergence <- 2
        MLEobj.test$Ey <- Eytest
        return(MLEobj.test)
      }
    }

    # fit and not all parameters estimated
    if (fit && !model.is.fixed) {
      if (silent == 2) cat(paste("Fitting model with ", method, ".\n", sep = ""))
      ## Fit and add param estimates to the object
      MLEobj <- try(MARSSfit(MLEobj), silent = TRUE)
      if (inherits(MLEobj, "try-error")) {
        cat(paste("Error: Stopped in MARSS() when trying to fit.  Try control$trace=1 for more information. You can also try another fitting method by passing in the method argument. This is especially helpful if you are seeing a chol error.\n\n", MLEobj$condition))
        MLEobj.test$convergence <- 2
        return(MLEobj.test)
      }
    }

    ## Add AIC and AICc and coef to the object
    ## Add states.se and ytT.se if no errors.  Return kf and Ey if trace>0
    if (MLEobj$convergence %in% c(0, 1, 3) || (MLEobj$convergence %in% c(10, 11) && MLEobj$method %in% kem.methods)) {
      kf <- MARSSkf(MLEobj) # use kf function requested by user; default smoother=TRUE
      if (silent == 2) cat("Adding logLik, AIC and coefficients.\n")
      MLEobj$logLik <- kf$logLik
      MLEobj <- MARSSaic(MLEobj)
      MLEobj$coef <- coef(MLEobj, type = "vector")
      if (silent == 2) cat("Adding states and states.se.\n")
      MLEobj$states <- kf$xtT
      if (!is.null(kf[["VtT"]])) {
        m <- attr(MLEobj$marss, "model.dims")[["x"]][1]
        TT <- attr(MLEobj$marss, "model.dims")[["data"]][2]
        states.se <- apply(kf[["VtT"]], 3, function(x) takediag(x))
        # KFS
        if (any(states.se < 0 & states.se > -1 * sqrt(.Machine$double.eps))) {
          states.se[states.se < 0 & states.se > -1 * sqrt(.Machine$double.eps)] <- 0
          MLEobj$errors <- c(MLEobj$errors, paste("\nAlert:", MLEobj$fun.kf, "returned negative values close to machine tolerance on diagonal of VtT and the states.se values for these are set to 0. You may want to use", ifelse(MLEobj$fun.kf == "MARSSkfas", "MARSSkfss()", "MARSSkfas()"), "to compute the states standard errors (states.se). See MARSSinfo('negVt') for insight.\n"))
        }
        if (any(states.se < 0)) {
          states.se[states.se < 0] <- NA
          MLEobj$errors <- c(MLEobj$errors, paste("\nAlert:", MLEobj$fun.kf, "returned negative values on diagonal of VtT and the states.se values for these are set to NA. You may want to use", ifelse(MLEobj$fun.kf == "MARSSkfas", "MARSSkfss()", "MARSSkfas()"), "to compute the states standard errors (states.se). See MARSSinfo('negVt') for insight.\n"))
        }
        states.se <- sqrt(states.se)
        if (m == 1) states.se <- matrix(states.se, 1, TT)
        rownames(states.se) <- attr(MLEobj$marss, "X.names")
      } else {
        states.se <- NULL
      }
      MLEobj[["states.se"]] <- states.se
      Ey <- MARSShatyt(MLEobj)
      MLEobj$ytT <- Ey[["ytT"]]
      if (!is.null(Ey[["OtT"]])) {
        n <- attr(MLEobj$marss, "model.dims")[["y"]][1]
        TT <- attr(MLEobj$marss, "model.dims")[["data"]][2]
        if (n == 1) yy <- matrix(Ey[["OtT"]][, , 1:TT], nrow = 1)
        if (n > 1) {
          yy <- apply(Ey[["OtT"]], 3, function(x) {
            takediag(x)
          })
        }
        ytT.se <- sqrt(yy - MLEobj$ytT^2)
        rownames(ytT.se) <- attr(MLEobj$marss, "Y.names")
      } else {
        ytT.se <- NULL
      }
      MLEobj$ytT.se <- ytT.se

      # Return kf and Ey, if trace = 1
      if (MLEobj$control$trace > 0) {
        if (fun.kf == "MARSSkfas") {
          kfss <- try(MARSSkfss(MLEobj, smoother = FALSE), silent = TRUE)
          if (inherits(kfss, "try-error") || !kfss$ok) {
            msg <- c("Not available. MARSSkfss() returned error.", kfss$condition)
            kfss <- list(Innov = msg, Sigma = msg, J = msg, Kt = msg)
          }
        } else {
          kfss <- kf
        }
        MLEobj$kf <- kf # from above will use function requested by user
        MLEobj$Ey <- Ey # from above
        # these are only returned by MARSSkfss
        MLEobj$Innov <- kfss$Innov
        MLEobj$Sigma <- kfss$Sigma
        MLEobj$J <- kfss$J
        MLEobj$Kt <- kfss$Kt
        if (fun.kf == "MARSSkfss") {
          MLEobj$J0 <- kfss$J0
        } else {
          # From kfss smoother so won't be available if fun.kf=MARSSkfas
          # Line above used smoother = FALSE; here smoother = TRUE
          J0 <- try(MARSSkfss(MLEobj), silent = TRUE)
          if (!inherits(J0, "try-error") && J0$ok) MLEobj$J0 <- J0$J0 else MLEobj$J0 <- "Not available. MARSSkfss() smoother returned error."
        }
      }
      # apply X and Y names various X and Y related elements
      MLEobj <- MARSSapplynames(MLEobj)
    }
    # END Adding info to output #################

    if ((!silent || silent == 2) && MLEobj[["convergence"]] %in% c(0, 1, 3, 10, 11, 12, 54)) {
      print(MLEobj)
    }
    if ((!silent || silent == 2) && !(MLEobj[["convergence"]] %in% c(0, 1, 3, 10, 11, 12, 54))) {
      cat(MLEobj$errors)
      if (!is.null(MLEobj$iter.record$message)) cat("Optimizer message: ", MLEobj$iter.record$message)
    } # 3 added since doesn't print if fit=FALSE
    if ((!silent || silent == 2) && !fit) {
      print(MLEobj$model)
    }

    return(MLEobj)
  } # end MLE methods

  return("method allowed but it's not in kem.methods, optim.methods or nlminb.methods so marssMLE object was not created")
}
