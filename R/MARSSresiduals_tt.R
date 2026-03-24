#' MARSS Contemporaneous Residuals
#'
#' @description
#' Calculates the standardized (or auxiliary) contemporaneous residuals, aka
#' the residuals and their variance conditioned on the data up to time \eqn{t}.
#' Contemporaneous residuals are only for the observations. Not exported.
#' Access this function with `MARSSresiduals(object, type="tt")`.
#'
#' @param object An object of class [marssMLE].
#' @param method Algorithm to use. Currently only "SS".
#' @param normalize TRUE/FALSE See details.
#' @param silent If TRUE, don't print inversion warnings.
#' @param fun.kf Can be ignored. This will change the Kalman filter/smoother
#'   function from the value in `object$fun.kf` if desired.
#'
#' @return
#' A list with the following components:
#'
#' * `model.residuals`: The observed contemporaneous model residuals: data minus
#'   the model predictions conditioned on the data 1 to t. A n x T matrix. NAs
#'   will appear where the data are missing.
#' * `state.residuals`: All NA. There are no contemporaneous residuals for the
#'   states.
#' * `residuals`: The residuals. `model.residuals` are in rows 1:n and
#'   `state.residuals` are in rows n+1:n+m.
#' * `var.residuals`: The joint variance of the residuals conditioned on
#'   observed data from 1 to t. This only has values in the 1:n,1:n upper block
#'   for the model residuals.
#' * `std.residuals`: The Cholesky standardized residuals as a n+m x T matrix.
#'   Rows n+1:n+m are all NA.
#' * `mar.residuals`: The marginal standardized residuals as a n+m x T matrix.
#' * `bchol.residuals`: Because state residuals do not exist, this will be
#'   equivalent to the Cholesky standardized residuals, `std.residuals`.
#' * `E.obs.residuals`: The expected value of the model residuals conditioned on
#'   the observed data 1 to t. Returned as a n x T matrix.
#' * `var.obs.residuals`: The variance of the model residuals conditioned on the
#'   observed data. Returned as a n x n x T matrix. For observed data, this will
#'   be 0. See [MARSSresiduals.tT()] for a discussion.
#' * `msg`: Any warning messages.
#'
#' @details
#' This function returns the conditional expected value (mean) and variance of
#' the model contemporaneous residuals. 'conditional' means conditioned on the
#' observed data up to time \eqn{t} and a set of parameters.
#'
#' **Model residuals**
#'
#' \eqn{\mathbf{v}_t}{v(t)} is the difference between the data and the
#' predicted data at time \eqn{t} given \eqn{\mathbf{x}_t}{x(t)}:
#' \deqn{ \mathbf{v}_t = \mathbf{y}_t - \mathbf{Z} \mathbf{x}_t - \mathbf{a} - \mathbf{d}\mathbf{d}_{t}}{ v(t) = y(t) - Z x(t) - a - D d(t)}
#' The observed model residuals use the data up to time \eqn{t}:
#' \deqn{ \hat{\mathbf{v}}_t = \mathbf{y}_t - \mathbf{Z}\mathbf{x}_t^{t} - \mathbf{a} - \mathbf{D}\mathbf{d}_{t}}{ hatv(t) = y(t) - Z xtt - a - D d(t)}
#'
#' The conditional variance is:
#' \deqn{ \hat{\Sigma}_t = \mathbf{R}+\mathbf{Z} \mathbf{V}_t^{t} \mathbf{Z}^\top }{hatSigma(t) = R + Z Vtt t(Z)}
#'
#' **Normalized residuals**
#'
#' If `normalize=FALSE`, the model is:
#' \deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{v}_t}{ y(t) = Z x(t) + a + v(t)}
#' \deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{w}_t}{ x(t) = B x(t-1) + u + w(t)}
#' If `normalize=TRUE`:
#' \deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{H}\mathbf{v}_t}{ y(t) = Z x(t) + a + Hv(t)}
#' \deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{G}\mathbf{w}_t}{ x(t) = B x(t-1) + u + Gw(t)}
#' with the variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)}
#' equal to \eqn{\mathbf{I}}{I} (identity).
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSresiduals.tT()], [MARSSresiduals.tt1()], [fitted.marssMLE()], [plot.marssMLE()]
#'
#' @examples
#' dat <- t(harborSeal)
#' dat <- dat[c(2, 11), ]
#' fit <- MARSS(dat)
#'
#' # Returns a matrix
#' MARSSresiduals(fit, type = "tt")$std.residuals
#' # Returns a data frame in long form
#' residuals(fit, type = "tt")
#'
#' @references
#' Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models.
#' Technical Report. arXiv:1411.0045.
#' @export
MARSSresiduals.tt <- function(object, method = c("SS"), normalize = FALSE, silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  # These are the residuals and their variance conditioned on the data up to time t
  # state residuals do not exist for this case

  ######################################
  # Set up variables
  MLEobj <- object
  if (missing(fun.kf)) {
    fun.kf <- MLEobj$fun.kf
  } else {
    MLEobj$fun.kf <- fun.kf
    # to ensure that MARSSkf, MARSShatyt and fitted() use the spec'd fun
  }
  method <- match.arg(method)
  model.dims <- attr(MLEobj$marss, "model.dims")
  TT <- model.dims[["x"]][2]
  m <- model.dims[["x"]][1]
  n <- model.dims[["y"]][1]
  y <- MLEobj$marss$data
  # set up holders
  et <- st.et <- mar.st.et <- matrix(NA, n + m, TT)
  model.et <- model.st.et <- model.mar.st.et <- matrix(NA, n, TT)
  model.var.et <- array(0, dim = c(n, n, TT))
  var.et <- array(0, dim = c(n + m, n + m, TT))
  msg <- NULL

  #### list of time-varying parameters
  time.varying <- is.timevarying(MLEobj)

  kf <- MARSSkf(MLEobj)
  Ey <- MARSShatyt(MLEobj, only.kem = FALSE)
  Rt <- parmat(MLEobj, "R", t = 1)$R # returns matrix
  Ht <- parmat(MLEobj, "H", t = 1)$H
  Rt <- Ht %*% tcrossprod(Rt, Ht)
  Zt <- parmat(MLEobj, "Z", t = 1)$Z
  Qtp <- parmat(MLEobj, "Q", t = 2)$Q

  if (method == "SS") {
    # We could set model.et to 0 where no data, but Kt will have a 0 column
    # for any missing y.
    model.et <- Ey$ytt - fitted(MLEobj, type = "ytt", output = "matrix") # model residuals
    et[1:n, ] <- model.et

    cov.et <- matrix(NA, n, m)
    for (t in 1:TT) {
      # model residuals
      if (time.varying$R) Rt <- parmat(MLEobj, "R", t = t)$R # returns matrix
      if (time.varying$H) Ht <- parmat(MLEobj, "H", t = t)$H
      if (time.varying$R || time.varying$H) Rt <- Ht %*% tcrossprod(Rt, Ht)
      if (time.varying$Z) Zt <- parmat(MLEobj, "Z", t = t)$Z
      # compute the variance of the residuals and state.et
      St <- Ey$yxtt[, , t] - tcrossprod(Ey$ytt[, t, drop = FALSE], kf$xtt[, t, drop = FALSE])
      model.var.et[, , t] <- Rt - tcrossprod(Zt %*% kf$Vtt[, , t], Zt) + tcrossprod(St, Zt) + tcrossprod(Zt, St)

      tmpvar.state.et <- matrix(NA, m, m)
      var.et[1:n, , t] <- cbind(model.var.et[, , t], cov.et)
      var.et[(n + 1):(n + m), , t] <- cbind(t(cov.et), tmpvar.state.et)

      if (normalize) {
        Rinv <- psolve(t(pchol(Rt)))
        model.et[, t] <- Rinv %*% model.et[, t]
        et[1:n, t] <- model.et[, t]
        model.var.et[, , t] <- tcrossprod(Rinv %*% model.var.et[, , t], Rinv)
        var.et[1:n, , t] <- cbind(model.var.et[, , t], cov.et)
      }
    }
  }

  # prepare standardized residuals; only model residuals since state residuals don't exist for this case
  for (t in 1:TT) {
    tmpvar <- sub3D(model.var.et, t = t)
    resids <- model.et[, t, drop = FALSE]
    # don't include values for resids if there is no residual (no data)
    # replace NAs with 0s
    is.miss <- is.na(y[, t])
    resids[is.miss] <- 0

    tmpvar[abs(tmpvar) < sqrt(.Machine$double.eps)] <- 0

    # inverse of diagonal of variance matrix for marginal standardization
    tmpvarinv <- try(psolve(makediag(takediag(tmpvar))), silent = TRUE)
    if (inherits(tmpvarinv, "try-error")) {
      model.mar.st.et[, t] <- NA
      msg <- c(msg, paste('MARSSresiduals.tt warning: the diagonal matrix of the variance of the residuals at t =", t, "is not invertible.  NAs returned for mar.residuals at t =", t, "\n'))
    } else { # inverse of diagonal ok
      model.mar.st.et[, t] <- sqrt(tmpvarinv) %*% resids
      model.mar.st.et[is.miss, t] <- NA
    }

    # psolve and pchol deal with 0s on diagonal
    # wrapped in try to prevent crashing if inversion not possible
    tmpchol <- try(pchol(tmpvar), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      model.st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
      next
    }
    # chol is ok
    # chol() returns the upper triangle. We need to lower triangle
    tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
    if (inherits(tmpcholinv, "try-error")) {
      model.st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      next
    }
    # both chol and inverse are ok
    model.st.et[, t] <- tmpcholinv %*% resids
    model.st.et[is.miss, t] <- NA
  }
  st.et[1:n, ] <- model.st.et
  mar.st.et[1:n, ] <- model.mar.st.et

  # Do not need to NA out the state residuals since they are all NA in this case

  # et is the expected value of the residuals conditioned on y(1,t)-the observed data up to time t
  E.obs.v <- et[1:n, , drop = FALSE]
  var.obs.v <- array(0, dim = c(n, n, TT))
  # this will be 0 for observed data
  for (t in 1:TT) var.obs.v[, , t] <- Ey$Ott[, , t] - tcrossprod(Ey$ytt[, t])

  # the observed model residuals are data - E(data), so NA for missing data.
  model.et <- et[1:n, , drop = FALSE]
  model.et[is.na(y)] <- NA
  et[1:n, ] <- model.et

  # add rownames
  Y.names <- attr(MLEobj$model, "Y.names")
  X.names <- attr(MLEobj$model, "X.names")
  rownames(et) <- rownames(st.et) <- rownames(mar.st.et) <- rownames(var.et) <- colnames(var.et) <- c(Y.names, X.names)
  rownames(E.obs.v) <- Y.names
  rownames(var.obs.v) <- colnames(var.obs.v) <- Y.names

  # output any warnings
  if (!is.null(msg) && object[["control"]][["trace"]] >= 0 & !silent) cat("MARSSresiduals.tt reported warnings. See msg element of returned residuals object.\n")

  return(list(
    model.residuals = et[1:n, , drop = FALSE],
    state.residuals = et[(n + 1):(n + m), , drop = FALSE],
    residuals = et,
    var.residuals = var.et,
    std.residuals = st.et,
    mar.residuals = mar.st.et,
    bchol.residuals = st.et,
    E.obs.residuals = E.obs.v,
    var.obs.residuals = var.obs.v,
    msg = msg
  ))
}
