#' MARSS One-Step-Ahead Residuals
#'
#' @description
#' Calculates the standardized (or auxiliary) one-step-ahead residuals, aka the
#' innovations residuals and their variance. Not exported. Access this function
#' with `MARSSresiduals(object, type="tt1")`. To get the residuals as a data
#' frame in long-form, use [residuals()][residuals.marssMLE] with
#' `type="tt1"`.
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
#' * `model.residuals`: The observed one-step-ahead model residuals: data minus
#'   the model predictions conditioned on the data \eqn{t=1} to \eqn{t-1}.
#'   These are termed innovations. A n x T matrix. NAs will appear where the
#'   data are missing.
#' * `state.residuals`: The one-step-ahead state residuals
#'   \eqn{ \mathbf{x}_{t+1}^{t+1} - \mathbf{B}\mathbf{x}_{t}^t - \mathbf{u} }{ xtt(t+1) - B xtt(t) - u}.
#'   Note, state residual at time \eqn{t} is the transition from time \eqn{t=t}
#'   to \eqn{t+1}.
#' * `residuals`: The residuals conditioned on the observed data up to time
#'   \eqn{t-1}. Returned as a (n+m) x T matrix with `model.residuals` in rows
#'   1 to n and `state.residuals` in rows n+1 to n+m. NAs will appear in rows
#'   1 to n in the places where data are missing.
#' * `var.residuals`: The joint variance of the one-step-ahead residuals.
#'   Returned as a n+m x n+m x T matrix.
#' * `std.residuals`: The Cholesky standardized residuals as a n+m x T matrix.
#'   This is `residuals` multiplied by the inverse of the lower triangle of the
#'   Cholesky decomposition of `var.residuals`. The model standardized residuals
#'   associated with the missing data are replaced with NA.
#' * `mar.residuals`: The marginal standardized residuals as a n+m x T matrix.
#'   This is `residuals` multiplied by the inverse of the diagonal matrix formed
#'   by the square-root of the diagonal of `var.residuals`. The model marginal
#'   residuals associated with the missing data are replaced with NA.
#' * `bchol.residuals`: The Block Cholesky standardized residuals as a (n+m) x
#'   T matrix.
#' * `E.obs.residuals`: The expected value of the model residuals conditioned on
#'   the observed data \eqn{t=1} to \eqn{t-1}. Returned as a n x T matrix.
#'   This will be all 0s. Included for completeness.
#' * `var.obs.residuals`: For one-step-ahead residuals, this will be the same
#'   as the 1:n, 1:n upper diagonal block in `var.residuals`. Included for
#'   completeness and as a code check.
#' * `msg`: Any warning messages.
#'
#' @details
#' This function returns the conditional expected value (mean) and variance of
#' the one-step-ahead residuals. 'conditional' means conditioned on the observed
#' data up to time \eqn{t-1} and a set of parameters.
#'
#' **Model residuals**
#'
#' \eqn{\mathbf{v}_t}{v_t} is the difference between the data and the predicted
#' data at time \eqn{t} given \eqn{\mathbf{x}_t}{x(t)}:
#' \deqn{ \mathbf{v}_t = \mathbf{y}_t - \mathbf{Z} \mathbf{x}_t - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ v(t) = y(t) - Z x(t) - a - D d(t)}
#' The observed model residuals \eqn{\hat{\mathbf{v}}_t}{hatv(t)} use the data
#' up to time \eqn{t-1}:
#' \deqn{ \hat{\mathbf{v}}_t = \mathbf{y}_t - \mathbf{Z}\mathbf{x}_t^{t-1} - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ hatv(t) = y(t) - Z xtt1(t) - a - D d(t)}
#'
#' **State residuals**
#'
#' The estimated state residuals:
#' \deqn{ \hat{\mathbf{w}}_{t+1} = \mathbf{x}_{t+1}^{t+1} - \mathbf{B}\mathbf{x}_{t}^t - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ hatw(t+1) = xtt(t+1) - B xtt(t) - u - C c(t+1)}
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
#' @seealso [MARSSresiduals.tT()], [MARSSresiduals.tt()], [fitted.marssMLE()], [plot.marssMLE()]
#'
#' @examples
#' dat <- t(harborSeal)
#' dat <- dat[c(2, 11), ]
#' fit <- MARSS(dat)
#'
#' MARSSresiduals(fit, type = "tt1")$std.residuals
#' residuals(fit, type = "tt1")
#'
#' @references
#' R. H. Shumway and D. S. Stoffer (2006). Section on the calculation of the
#' likelihood of state-space models in Time series analysis and its
#' applications. Springer-Verlag, New York.
#'
#' Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models.
#' Technical Report. arXiv:1411.0045.
#' @export
MARSSresiduals.tt1 <- function(object, method = c("SS"), normalize = FALSE, silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  # These are the residuals and their variance conditioned on the data up to time t-1

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
  et <- st.et <- mar.st.et <- bchol.st.et <- matrix(NA, n + m, TT)
  model.et <- matrix(NA, n, TT)
  model.var.et <- array(0, dim = c(n, n, TT))
  var.et <- array(0, dim = c(n + m, n + m, TT))
  msg <- NULL

  #### list of time-varying parameters
  time.varying <- is.timevarying(MLEobj)

  kf <- MARSSkf(MLEobj)
  Kt.ok <- TRUE
  if (MLEobj$fun.kf == "MARSSkfas") {
    Kt <- try(MARSSkfss(MLEobj), silent = TRUE)
    if (inherits(Kt, "try-error") || !Kt$ok) {
      Kt <- array(0, dim = c(m, n, TT))
      Kt.ok <- FALSE
    } else {
      Kt <- Kt$Kt
    }
    kf$Kt <- Kt
  }
  Ey <- MARSShatyt(MLEobj, only.kem = FALSE)
  Rt <- parmat(MLEobj, "R", t = 1)$R # returns matrix
  Ht <- parmat(MLEobj, "H", t = 1)$H
  Rt <- Ht %*% tcrossprod(Rt, Ht)
  Zt <- parmat(MLEobj, "Z", t = 1)$Z
  Qtp <- parmat(MLEobj, "Q", t = 2)$Q
  ######################################

  ######################################
  # Compute residuals
  if (method == "SS") {
    # We could set model.et to 0 where no data, but Kt will have a 0 column
    # for any missing y.
    model.et <- Ey$ytt - fitted(MLEobj, type = "ytt1", output = "matrix") # model residuals
    et[1:n, ] <- model.et

    cov.et <- matrix(0, n, m)
    # get the model residual variance
    for (t in 1:TT) {
      # model residuals
      if (time.varying$R) Rt <- parmat(MLEobj, "R", t = t)$R # returns matrix
      if (time.varying$H) Ht <- parmat(MLEobj, "H", t = t)$H
      if (time.varying$R || time.varying$H) Rt <- Ht %*% tcrossprod(Rt, Ht)
      if (time.varying$Z) Zt <- parmat(MLEobj, "Z", t = t)$Z
      model.var.et[, , t] <- Rt + tcrossprod(Zt %*% kf$Vtt1[, , t], Zt)
    }
    # Then the states since that needs model var at t+1
    for (t in 1:TT) {
      if (t < TT) {
        Ktp <- sub3D(kf$Kt, t = t + 1)
        et[(n + 1):(n + m), t] <- Ktp %*% model.et[, t + 1]
        tmpvar.state.et <- tcrossprod(Ktp %*% model.var.et[, , t + 1], Ktp)
      } else {
        tmpvar.state.et <- matrix(0, m, m)
      }
      var.et[1:n, , t] <- cbind(model.var.et[, , t], cov.et)
      var.et[(n + 1):(n + m), , t] <- cbind(t(cov.et), tmpvar.state.et)

      if (normalize) {
        Qpinv <- matrix(0, m, n + m)
        Rinv <- matrix(0, n, n + m)
        if (t < TT) {
          if (time.varying$Q) Qtp <- parmat(MLEobj, "Q", t = t + 1)$Q
          Qpinv[, (n + 1):(n + m)] <- psolve(t(pchol(Qtp)))
        }
        Rinv[, 1:n] <- psolve(t(pchol(Rt)))
        RQinv <- rbind(Rinv, Qpinv) # block diag matrix
        et[, t] <- RQinv %*% et[, t]
        var.et[, , t] <- tcrossprod(RQinv %*% var.et[, , t], RQinv)
      }
    }
  }
  ######################################

  ######################################
  # prepare standardized residuals
  for (t in 1:TT) {
    tmpvar <- sub3D(var.et, t = t)
    resids <- et[, t, drop = FALSE]
    # don't include values for resids if there is no residual (no data)
    # replace NAs with 0s
    is.miss <- c(is.na(y[, t]), rep(FALSE, m))
    resids[is.miss] <- 0

    tmpvar[abs(tmpvar) < sqrt(.Machine$double.eps)] <- 0

    # Marginal
    # inverse of diagonal of variance matrix for marginal standardization
    tmpvarinv <- try(psolve(makediag(takediag(tmpvar))), silent = TRUE)
    if (inherits(tmpvarinv, "try-error")) {
      mar.st.et[, t] <- NA
      msg <- c(msg, paste('MARSSresiduals.tt1 warning: the diagonal matrix of the variance of the residuals at t =", t, "is not invertible.  NAs returned for mar.residuals at t =", t, "\n'))
    } else { # inv of diagonal ok, can compute marginal std residuals
      mar.st.et[, t] <- sqrt(tmpvarinv) %*% resids
      mar.st.et[is.miss, t] <- NA
    }

    # Block Cholesky
    # psolve and pchol deal with 0s on diagonal
    tmpchol <- try(pchol(tmpvar[(n + 1):(n + m), (n + 1):(n + m), drop = FALSE]), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      bchol.st.et[(n + 1):(n + m), t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tT warning: the variance of the state residuals at t =", t, "is not invertible.  NAs returned for bchol.std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
    } else {
      # chol() returns the upper triangle. We need to lower triangle to t()
      tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
      if (inherits(tmpcholinv, "try-error")) {
        bchol.st.et[(n + 1):(n + m), t] <- NA
        msg <- c(msg, paste("MARSSresiduals.tT warning: the variance of the state residuals at t =", t, "is not invertible.  NAs returned for bchol.std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      } else {
        bchol.st.et[(n + 1):(n + m), t] <- tmpcholinv %*% resids[(n + 1):(n + m), 1, drop = FALSE]
      }
    }

    # Cholesky
    # psolve and pchol deal with 0s on diagonal
    tmpchol <- try(pchol(tmpvar), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt1 warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
      next # got to next t since next inversion will fail
    }
    # chol() returns the upper triangle. We need to lower triangle
    tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
    if (inherits(tmpcholinv, "try-error")) {
      st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tt1 warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      next # go to next t
    }
    # all ok, can compute Cholesky std.residuals
    st.et[, t] <- tmpcholinv %*% resids
    st.et[is.miss, t] <- NA
  }
  bchol.st.et[1:n, ] <- st.et[1:n, ] # because the upper right block of the lower tri is 0
  ######################################

  ######################################
  # the state.residual at the last time step is NA because it is x(T+1) - f(x(T)) and T+1 does not exist.  For the same reason, the var.residuals at TT will have NAs
  et[(n + 1):(n + m), TT] <- NA
  var.et[, (n + 1):(n + m), TT] <- NA
  var.et[(n + 1):(n + m), , TT] <- NA
  st.et[, TT] <- NA
  mar.st.et[(n + 1):(n + m), TT] <- NA
  ######################################

  ######################################
  # et is the expected value of the residuals conditioned on y(1)-the observed data
  E.obs.v <- et[1:n, , drop = FALSE]
  var.obs.v <- array(0, dim = c(n, n, TT))
  # this will be 0 for observed data
  for (t in 1:TT) var.obs.v[, , t] <- Ey$Ott1[, , t] - tcrossprod(Ey$ytt1[, t])
  ######################################

  ######################################
  # the observed model residuals are data - E(data), so NA for missing data.
  model.et <- et[1:n, , drop = FALSE]
  model.et[is.na(y)] <- NA
  et[1:n, ] <- model.et
  ######################################

  ######################################
  # add rownames
  Y.names <- attr(MLEobj$model, "Y.names")
  X.names <- attr(MLEobj$model, "X.names")
  rownames(et) <- rownames(st.et) <- rownames(mar.st.et) <- rownames(bchol.st.et) <- rownames(var.et) <- colnames(var.et) <- c(Y.names, X.names)
  rownames(E.obs.v) <- Y.names
  rownames(var.obs.v) <- colnames(var.obs.v) <- Y.names
  ######################################

  ######################################
  # output any warnings
  if (!is.null(msg) && object[["control"]][["trace"]] >= 0 & !silent) cat("MARSSresiduals.tt1 reported warnings. See msg element of returned residuals object.\n")
  ######################################

  if (!Kt.ok) {
    et[(n + 1):(n + m), ] <- NA
    var.et[(n + 1):(n + m), (n + 1):(n + m), ] <- NA
    st.et[(n + 1):(n + m), ] <- NA
    mar.st.et[(n + 1):(n + m), ] <- NA
    bchol.st.et[(n + 1):(n + m), ] <- NA
  }
  ret <- list(
    model.residuals = et[1:n, , drop = FALSE],
    state.residuals = et[(n + 1):(n + m), , drop = FALSE],
    residuals = et,
    var.residuals = var.et,
    std.residuals = st.et,
    mar.residuals = mar.st.et,
    bchol.residuals = bchol.st.et,
    E.obs.residuals = E.obs.v,
    var.obs.residuals = var.obs.v,
    msg = msg
  )
  return(ret)
}
