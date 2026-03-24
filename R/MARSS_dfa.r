###################################################################################
# Helper functions to create and work with a DFA model sensu Zuur
# x(t)=x(t-1) + w(t), W~MVN(0,1)
# y(t)=Z x(t) + A(t) + D(t) d(t) + v(t), V~MVN(0,R)
# x(t0) = x0 + l, L ~ MVN(0,5)

# MARSS.dfa: The conversion functions have 2 parts
# Part 1 Set up the DFA model in MARSS.marxss form
# Part 2 Call MARSS.marxss to finish the set-up and checking

# Functions that are called by the generic functions print.marssMLE, MARSSresiduals,
# predict.marssMLE, coef.marssMLE, MARSSinits.marssMLE
# print_dfa, residuals_dfa, predict_dfa, coef_dfa, MARSSinits_dfa
###################################################################################

#' Multivariate Dynamic Factor Analysis
#'
#' @description
#' The Dynamic Factor Analysis model in MARSS is specified via
#' `form="dfa"` in a [MARSS()] function call. This is a MARSS(1) model of
#' the form:
#' \deqn{\mathbf{x}_{t} = \mathbf{x}_{t-1} + \mathbf{w}_t, \textrm{ where } \mathbf{W}_t \sim \textrm{MVN}(0,\mathbf{I})}{x(t) = x(t-1) + w(t), where W(t) ~ MVN(0,I)}
#' \deqn{\mathbf{y}_t = \mathbf{Z}_t \mathbf{x}_t + \mathbf{D}_t \mathbf{d}_t + \mathbf{v}_t, \textrm{ where } \mathbf{V}_t \sim \textrm{MVN}(0,\mathbf{R}_t)}{y(t) = Z(t) x(t) + D(t) d(t) + v(t), where V(t) ~ MVN(0,R(t))}
#' \deqn{\mathbf{X}_1 \sim \textrm{MVN}(\mathbf{x0}, 5\mathbf{I})}{X(1) ~ MVN(x0, 5I) }
#' Note, by default \eqn{\mathbf{x}_1}{x(1)} is treated as a diffuse prior.
#'
#' Passing in `form="dfa"` to [MARSS()] invokes a helper function to create
#' that model and creates the \eqn{\mathbf{Z}}{Z} matrix for the user.
#' \eqn{\mathbf{Q}}{Q} is by definition identity, \eqn{\mathbf{x}_0}{x0} is
#' zero and \eqn{\mathbf{V_0}}{V0} is diagonal with large variance (5).
#' \eqn{\mathbf{u}}{U} is zero, \eqn{\mathbf{a}}{A} is zero, and covariates
#' only enter the \eqn{\mathbf{y}}{Y} equation. Because \eqn{\mathbf{u}}{U}
#' and \eqn{\mathbf{a}}{A} are 0, the data should have mean 0 (demeaned)
#' otherwise one is likely to be creating a structurally inadequate model
#' (i.e. the model implies that the data have mean = 0, yet data do not have
#' mean = 0).
#'
#' @section Usage:
#' ```
#' MARSS(y,
#'     inits = NULL,
#'     model = NULL,
#'     miss.value = as.numeric(NA),
#'     method = "kem",
#'     form = "dfa",
#'     fit = TRUE,
#'     silent = FALSE,
#'     control = NULL,
#'     fun.kf = "MARSSkfas",
#'     demean = TRUE,
#'     z.score = TRUE)
#' ```
#'
#' @param MARSS.call A list of arguments from a [MARSS()] call with
#'   `form="dfa"`.
#'
#' @details
#' Some arguments are common to all forms: "y" (data), "inits", "control",
#' "method", "form", "fit", "silent", "fun.kf". See [MARSS()] for information
#' on these arguments.
#'
#' In addition, `form="dfa"` has some special arguments that can be passed in:
#'
#' * `demean` Logical. Default is TRUE, which means the data will be demeaned.
#' * `z.score` Logical. Default is TRUE, which means the data will be
#'   z-scored (demeaned and variance standardized to 1).
#' * `covariates` Covariates (\eqn{d}) for the \eqn{y} equation. No missing
#'   values allowed and must be a matrix with the same number of time steps as
#'   the data. An unconstrained \eqn{D} matrix will be estimated.
#'
#' The `model` argument of the [MARSS()] call is constrained in terms of what
#' parameters can be changed and how they can be changed. An additional
#' element, `m`, can be passed into the `model` argument that specifies the
#' number of hidden state variables. It is not necessary for the user to
#' specify `Z` as the helper function will create a `Z` appropriate for a DFA
#' model.
#'
#' The `model` argument is a list. The following details what list elements
#' can be passed in:
#'
#' * `B` "Identity". The standard (and default) DFA model has B="identity".
#'   However it can be "identity", "diagonal and equal", "diagonal and
#'   unequal" or a time-varying fixed or estimated diagonal matrix.
#' * `U` "Zero". Cannot be changed or passed in via model argument.
#' * `Q` "Identity". The standard (and default) DFA model has Q="identity".
#'   However, it can be "identity", "diagonal and equal", "diagonal and
#'   unequal" or a time-varying fixed or estimated diagonal matrix.
#' * `Z` Can be passed in as a (list) matrix if the user does not want a
#'   default DFA `Z` matrix. There are many equivalent ways to construct a
#'   DFA `Z` matrix. The default is Zuur et al.'s form (see User Guide).
#' * `A` Default="zero". Can be "unequal", "zero" or a matrix.
#' * `R` Default="diagonal and equal". Can be set to "identity", "zero",
#'   "unconstrained", "diagonal and unequal", "diagonal and equal",
#'   "equalvarcov", or a (list) matrix to specify general forms.
#' * `x0` Default="zero". Can be "unconstrained", "unequal", "zero", or a
#'   (list) matrix.
#' * `V0` Default=diagonal matrix with 5 on the diagonal. Can be "identity",
#'   "zero", or a matrix.
#' * `tinitx` Default=0. Can be 0 or 1. Tells MARSS whether x0 is at t=0 or
#'   t=1.
#' * `m` Default=1. Can be 1 to n (the number of y time-series). Must be
#'   integer.
#'
#' See the
#' [User Guide](https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf)
#' chapter on Dynamic Factor Analysis for examples of using `form="dfa"`.
#'
#' @return
#' An object of class [marssMLE]. See [print.marssMLE()] for a discussion of
#' the various output available for [marssMLE] objects (coefficients,
#' residuals, Kalman filter and smoother output, imputed values for missing
#' data, etc.). See [MARSSsimulate()] for simulating from [marssMLE] objects.
#' [MARSSboot()] for bootstrapping, [MARSSaic()] for calculation of various
#' AIC related model selection metrics, and [MARSSparamCIs()] for calculation
#' of confidence intervals and bias.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSS()], [MARSS.marxss()]
#'
#' @references
#' The MARSS User Guide: Holmes, E. E., E. J. Ward, and M. D. Scheuerell
#' (2012) Analysis of multivariate time-series using the MARSS package. NOAA
#' Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E.,
#' Seattle, WA 98112. Type `RShowDoc("UserGuide",package="MARSS")` to open a
#' copy.
#'
#' @examples
#' \dontrun{
#' dat <- t(harborSealWA[,-1])
#' # DFA with 3 states; used BFGS because it fits much faster for this model
#' fit <- MARSS(dat, model = list(m=3), form="dfa", method="BFGS")
#'
#' # See the Dynamic Factor Analysis chapter in the User Guide
#' RShowDoc("UserGuide", package = "MARSS")
#' }
#'
#' @export
MARSS.dfa <- function(MARSS.call) {
  # MARSS(data, model=list(), covariates=NULL, z.score=TRUE, demean=TRUE, control=list())
  # model.defaults =list(A="zero", R="diagonal and equal", D="zero", x0="zero", V0=diag(5,1), tinitx=0, diffuse=FALSE, m=1)

  # load needed package globals
  common.allowed.in.MARSS.call <- get("common.allowed.in.MARSS.call", envir = pkg_globals)


  # Part 1 Set up defaults and check that what the user passed in is allowed
  # 1 Check for form dependent user inputs for method and reset defaults for inits, and control if desired
  # 2 Specify the text shortcuts and whether factors or matrices can be passed in
  #   The names in the allowed list do not need to be A, B, Q .... as used in the marss form
  #   Other names can be used if you want the user to use those names; then in the MARSS.form function,
  #   you convert the user passed in names into the marss form with the A, B, Q, R, ... names
  #   checkModelList() will check what the user passes in against these allowed values, so
  #   so you need to make sure each name in model.defaults has a model.allowed value here

  # Check that no args were passed into MARSS that are not allowed
  dfa.allowed.in.MARSS.call <- c("model", "z.score", "demean", "covariates")
  allowed.in.call <- c(dfa.allowed.in.MARSS.call, common.allowed.in.MARSS.call)
  if (any(!(names(MARSS.call) %in% allowed.in.call))) {
    bad.names <- names(MARSS.call)[!(names(MARSS.call) %in% allowed.in.call)]
    msg <- paste("Argument ", paste(bad.names, collapse = ", "), "  not allowed MARSS call for form ", MARSS.call$form, ". See ?MARSS.dfa\n", sep = "")
    cat("\n", "Errors were caught in MARSS.dfa \n", msg, sep = "")
    stop("Stopped in MARSS.dfa() due to problem(s) with model specification.\n", call. = FALSE)
  }

  # Set up some defaults
  if (is.null(MARSS.call[["z.score"]])) MARSS.call$z.score <- TRUE
  if (is.null(MARSS.call[["demean"]])) MARSS.call$demean <- TRUE

  # Start error checking
  problem <- FALSE
  msg <- c()
  # check that data and covariates elements are matrix or vector, no dataframes, and is numeric
  for (el in c("data", "covariates"[!is.null(MARSS.call[["covariates"]])])) {
    if (!(is.matrix(MARSS.call[[el]]) || is.vector(MARSS.call[[el]]) || inherits(MARSS.call[[el]], "ts"))) {
      problem <- TRUE
      msg <- c(msg, paste(el, " must be a matrix, vector or ts/mts ojbect (not a data frame or 3D array).\n"))
    } else {
      if (is.vector(MARSS.call[[el]])) MARSS.call[[el]] <- matrix(MARSS.call[[el]], 1)
      if (inherits(MARSS.call[[el]], "ts")) MARSS.call[[el]] <- t(MARSS.call[[el]])
    }
  }
  if (!is.null(MARSS.call[["covariates"]])) {
    if (dim(MARSS.call[["data"]])[2] != dim(MARSS.call[["covariates"]])[2]) {
      problem <- TRUE
      msg <- c(msg, "data and covariates must have the same number of time steps.\n")
    }
  }
  if (!is.null(MARSS.call[["covariates"]])) {
    if (any(is.na(MARSS.call[["covariates"]]))) {
      problem <- TRUE
      msg <- c(msg, "covariates cannot have any missing values in a standard DFA.\n See User Guide section on DFA for alternate approaches when covariates have missing values.\n")
    }
  }
  if (!is.null(MARSS.call[["model"]][["m"]])) {
    if (length(MARSS.call[["model"]][["m"]]) != 1) {
      problem <- TRUE
      msg <- c(msg, "model$m must be an integer between 1 and n.\n")
    } else {
      if (!is.numeric(MARSS.call[["model"]][["m"]]) || !is.wholenumber(MARSS.call$model$m)) {
        problem <- TRUE
        msg <- c(msg, "model$m must be an integer between 1 and n.\n")
      } else {
        if (MARSS.call[["model"]][["m"]] > dim(MARSS.call[["data"]])[1]) {
          problem <- TRUE
          msg <- c(msg, "model$m must be an integer between 1 and n.\n")
        }
      }
    }
  }
  if (!is.null(MARSS.call[["model"]][["Z"]])) {
      problem <- TRUE
      msg <- c(msg, "If using form='dfa', specify the Z matrix using m (number of trends). If you need a custom Z, then use the default MARSS model and specify all your matrices. See ?MARSS.dfa for defaults for the matrices.\n")
  }
  

  if (problem) {
    cat("\n", "Errors were caught in MARSS.dfa \n", msg, sep = "")
    stop("Stopped in MARSS.dfa() due to specification problem(s).\n", call. = FALSE)
  }

  n <- dim(MARSS.call[["data"]])[1]
  model.allowed <- list(
    # if it is a length 1 vector then the value must be one of these.  All elements in your model list must be here
    A = c("unequal", "zero"),
    R = c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    D = c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    x0 = c("unconstrained", "unequal", "zero"),
    B = c("identity", "diagonal and equal", "diagonal and unequal"),
    Q = c("identity", "diagonal and equal", "diagonal and unequal"),
    V0 = c("identity", "zero"),
    tinitx = c(0, 1),
    diffuse = c(TRUE, FALSE),
    m = 1:n,
    # This line says what is allowed to be a matrix
    # Z matrix needs to be allowed here because a Z matrix is computed from
    # the m passed in by user
    matrices = c("Z", "A", "R", "D", "x0", "V0", "Q", "B")
  )
 
  ## Set-up model defaults
  if (!is.null(MARSS.call[["covariates"]])) D <- "unconstrained" else D <- "zero"
  # Set up m and Z
  if (is.null(MARSS.call[["model"]][["m"]])) m <- 1 else m <- MARSS.call[["model"]][["m"]]

  # Set up default Z
    Z <- matrix(list(), nrow = n, ncol = m)
    # insert row (i) & col (j) indices
    for (i in seq(n)) {
      Z[i, ] <- paste(i, seq(m), sep = "")
    }
    # set correct i,j values in Z to numeric 0
    if (m > 1) {
      for (i in 1:(m - 1)) {
        Z[i, (i + 1):m] <- 0
      }
    }
  
  # defaults for any missing model list elements
  model.defaults <- list(
    A = "zero",
    R = "diagonal and equal",
    D = D,
    B = "identity",
    Q = "identity",
    x0 = "zero",
    V0 = diag(5, m),
    tinitx = 0,
    diffuse = FALSE,
    m = 1,
    Z = Z
  )

  # This checks that what user passed in model list can be interpreted and converted to form marss
  # if no errors, it updates the model list by filling in missing elements with the defaults
  MARSS.call$model <- checkModelList(MARSS.call[["model"]], model.defaults, model.allowed)
  model <- MARSS.call[["model"]]

  if (!(MARSS.call$z.score %in% c(TRUE, FALSE))) {
    stop("Stopped in MARSS.dfa: z.score must be TRUE/FALSE.\n", call. = FALSE)
  }
  if (!(MARSS.call$demean %in% c(TRUE, FALSE))) {
    stop("Stopped in MARSS.dfa: demean must be TRUE/FALSE.\n", call. = FALSE)
  }

  # Set up U; always 0 for dfa
  U <- matrix(0, m, 1)

  # Set up D and d
  if (is.null(MARSS.call[["covariates"]])) d <- matrix(0, 1, 1) else d <- MARSS.call$covariates

  # Set up Q  & B; always fixed for dfa
  # Allow user to use diagonal matrices
  # Q=diag(1,m); B=diag(1,m)
  Q <- model$Q
  if (is.array(Q)) { # 2D or 3D
    if (length(dim(Q)) == 3) {
      Q.is.diagonal <- all(apply(Q, 3, is.diagonal))
    }
    if (length(dim(Q)) == 2) {
      Q.is.diagonal <- is.diagonal(Q)
    }
    if (!Q.is.diagonal) stop("Stopped in MARSS.dfa: Q must be diagonal.\n", call. = FALSE)
  }
  B <- model$B
  if (is.array(B)) { # 2D or 3D
    if (length(dim(B)) == 3) {
      B.is.diagonal <- all(apply(B, 3, is.diagonal))
    }
    if (length(dim(B)) == 2) {
      B.is.diagonal <- is.diagonal(B)
    }
    if (!B.is.diagonal) stop("Stopped in MARSS.dfa: B must be diagonal.\n", call. = FALSE)
  }

  # set up list of model components for a marxss model
  dfa.model <- list(Z = Z, A = model$A, D = model$D, d = d, R = model$R, B = model$B, U = U, Q = model$Q, x0 = model$x0, V0 = model$V0, tinitx = model$tinitx)

  dat <- MARSS.call[["data"]]
  if (MARSS.call$demean) {
    y.bar <- apply(dat, 1, mean, na.rm = TRUE)
    dat <- (dat - y.bar)
  }
  if (MARSS.call[["z.score"]]) {
    Sigma <- sqrt(apply(dat, 1, var, na.rm = TRUE))
    dat <- dat * (1 / Sigma)
  }

  MARSS.call <- list(data = dat, inits = MARSS.call$inits, control = MARSS.call$control, method = MARSS.call$method, form = "dfa", silent = MARSS.call$silent, fit = MARSS.call$fit, fun.kf = MARSS.call$fun.kf)

  # dfa is a type of marxss model, so use MARSS.marxss to test it and set up the marss object
  tmp <- MARSS.marxss(list(data = dat, model = dfa.model, method = MARSS.call$method, silent = MARSS.call$silent))
  # marss is the name for the form=marss model object that MARSS.form functions return
  # need to add "dfa" to attribute form
  marxss_object <- tmp$model
  attr(marxss_object, "form") <- c("dfa", "marxss")
  MARSS.call$model <- marxss_object
  MARSS.call$marss <- tmp$marss

  ## Return MARSS inputs as list
  MARSS.call
}
# This works since dfa just creates a marxss object so x$model is marssMODEL form=c("dfa", "marxss")
# probably want to customize for dfa later
print_dfa <- function(x) {
  return(print_marxss(x))
}
coef_dfa <- function(x) {
  return(coef_marxss(x))
}
MARSSinits_dfa <- function(MLEobj, inits) {
  return(MARSSinits_marxss(MLEobj, inits))
}
predict_dfa <- function(x, newdata, n.ahead, t.start) {
  predict_marxss(x, newdata, n.ahead, t.start)
}
describe_dfa <- function(MODELobj) {
  describe_marss(MODELobj)
}
is.marssMODEL_dfa <- function(MODELobj, method = "kem") {
  is.marssMODEL_marxss(MODELobj, method = method)
}
