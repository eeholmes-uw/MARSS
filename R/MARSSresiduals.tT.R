#' MARSS Smoothed Residuals
#'
#' Calculates the standardized (or auxiliary) smoothed residuals sensu Harvey,
#' Koopman and Penzer (1998). The expected values and variance for missing (or
#' left-out) data are also returned (Holmes 2014). Not exported. Access this
#' function with `MARSSresiduals(object, type="tT")`. At time \eqn{t} (in the
#' returned matrices), the model residuals are for time \eqn{t}, while the state
#' residuals are for the transition from \eqn{t} to \eqn{t+1} following the
#' convention in Harvey, Koopman and Penzer (1998).
#'
#' @param object An object of class [marssMLE()].
#' @param Harvey TRUE/FALSE. Use the Harvey et al. (1998) algorithm or use the
#'   Holmes (2014) algorithm. The values are the same except for missing values.
#' @param normalize TRUE/FALSE See details.
#' @param silent If TRUE, don't print inversion warnings.
#' @param fun.kf Kalman filter function to use. Can be ignored.
#'
#' @return A list with the following components
#' * `model.residuals`: The the observed smoothed model residuals: data minus the
#'   model predictions conditioned on all observed data. This is different than
#'   the Kalman filter innovations which use on the data up to time \eqn{t-1}
#'   for the predictions. See details.
#' * `state.residuals`: The smoothed state residuals
#'   \eqn{\mathbf{x}_{t+1}^T - \mathbf{Z} \mathbf{x}_{t}^T - \mathbf{u}}{xtT(t+1) - Z xtT(t) - u}.
#'   The last time step will be NA because the last step would be for T to T+1
#'   (past the end of the data).
#' * `residuals`: The residuals conditioned on the observed data. Returned as a
#'   (n+m) x T matrix with `model.residuals` in rows 1 to n and
#'   `state.residuals` in rows n+1 to n+m. NAs will appear in rows 1 to n in
#'   the places where data are missing.
#' * `var.residuals`: The joint variance of the model and state residuals
#'   conditioned on observed data. Returned as a (n+m) x (n+m) x T matrix. For
#'   Harvey=FALSE, this is Holmes (2014) equation 57. For Harvey=TRUE, this is
#'   the residual variance in eqn. 24, page 113, in Harvey et al. (1998). They
#'   are identical except for missing values, for those Harvey=TRUE returns 0s.
#'   For the state residual variance, the last time step will be all NA because
#'   the last step would be for T to T+1 (past the end of the data).
#' * `std.residuals`: The Cholesky standardized residuals as a (n+m) x T matrix.
#'   This is `residuals` multiplied by the inverse of the lower triangle of the
#'   Cholesky decomposition of `var.residuals`. The model standardized residuals
#'   associated with the missing data are replaced with NA.
#' * `mar.residuals`: The marginal standardized residuals as a (n+m) x T matrix.
#'   This is `residuals` multiplied by the inverse of the diagonal matrix formed
#'   by the square-root of the diagonal of `var.residuals`. The model marginal
#'   residuals associated with the missing data are replaced with NA.
#' * `bchol.residuals`: The Block Cholesky standardized residuals as a (n+m) x T
#'   matrix. This is `model.residuals` multiplied by the inverse of the lower
#'   triangle of the Cholesky decomposition of `var.residuals[1:n,1:n,]` and
#'   `state.residuals` multiplied by the inverse of the lower triangle of the
#'   Cholesky decomposition of
#'   `var.residuals[(n+1):(n+m),(n+1):(n+m),]`.
#' * `E.obs.residuals`: The expected value of the model residuals conditioned on
#'   the observed data. Returned as a n x T matrix. For observed data, this will
#'   be the observed residuals (values in `model.residuals`). For unobserved
#'   data, this will be 0 if \eqn{\mathbf{R}}{R} is diagonal but non-zero if
#'   \eqn{\mathbf{R}}{R} is non-diagonal. See details.
#' * `var.obs.residuals`: The variance of the model residuals conditioned on the
#'   observed data. Returned as a n x n x T matrix. For observed data, this will
#'   be 0. See details.
#' * `msg`: Any warning messages. This will be printed unless
#'   `Object$control$trace = -1` (suppress all error messages).
#'
#' @details
#'
#' This function returns the raw, the Cholesky standardized and the marginal
#' standardized smoothed model and state residuals. 'smoothed' means conditioned
#' on all the observed data and a set of parameters. These are the residuals
#' presented in Harvey, Koopman and Penzer (1998) pages 112-113, with the
#' addition of the values for unobserved data (Holmes 2014). If Harvey=TRUE, the
#' function uses the algorithm on page 112 of Harvey, Koopman and Penzer (1998)
#' to compute the conditional residuals and variance of the residuals. If
#' Harvey=FALSE, the function uses the equations in the technical report (Holmes
#' 2014). Unlike the innovations residuals, the smoothed residuals are
#' autocorrelated (section 4.1 in Harvey and Koopman 1992) and thus an ACF test
#' on these residuals would not reveal model inadequacy.
#'
#' The residuals matrix has a value for each time step. The residuals in column
#' \eqn{t} rows 1 to n are the model residuals associated with the data at time
#' \eqn{t}. The residuals in rows n+1 to n+m are the state residuals associated
#' with the transition from \eqn{\mathbf{x}_{t}}{x(t)} to
#' \eqn{\mathbf{x}_{t+1}}{x(t+1)}, not the transition from
#' \eqn{\mathbf{x}_{t-1}}{x(t-1)} to \eqn{\mathbf{x}_{t}}{x(t)}. Because
#' \eqn{\mathbf{x}_{t+1}}{x(t+1)} does not exist at time \eqn{T}, the state
#' residuals and associated variances at time \eqn{T} are NA.
#'
#' Below the conditional residuals and their variance are discussed. The random
#' variables are capitalized and the realizations from the random variables are
#' lower case. The random variables are \eqn{\mathbf{X}}{X},
#' \eqn{\mathbf{Y}}{Y}, \eqn{\mathbf{V}}{V} and \eqn{\mathbf{W}}{W}. There are
#' two types of \eqn{\mathbf{Y}}{Y}. The observed \eqn{\mathbf{Y}}{Y} that are
#' used to estimate the states \eqn{\mathbf{x}}{x}. These are termed
#' \eqn{\mathbf{Y}^{(1)}}{Y(1)}. The unobserved \eqn{\mathbf{Y}}{Y} are termed
#' \eqn{\mathbf{Y}^{(2)}}{Y(2)}. These are not used to estimate the states
#' \eqn{\mathbf{x}}{x} and we may or may not know the values of
#' \eqn{\mathbf{y}^{(2)}}{y(2)}. Typically we treat
#' \eqn{\mathbf{y}^{(2)}}{y(2)} as unknown but it may be known but we did not
#' include it in our model fitting. Note that the model parameters
#' \eqn{\Theta}{Theta} are treated as fixed or known. The 'fitting' does not
#' involve estimating \eqn{\Theta}{Theta}; it involves estimating
#' \eqn{\mathbf{x}}{x}. All MARSS parameters can be time varying but the
#' \eqn{t} subscripts are left off parameters to reduce clutter.
#'
#' **Model residuals**
#'
#' \eqn{\mathbf{v}_{t}}{v(t)} is the difference between the data and the
#' predicted data at time \eqn{t} given \eqn{\mathbf{x}_{t}}{x(t)}:
#' \deqn{ \mathbf{v}_{t} = \mathbf{y}_{t} - \mathbf{Z} \mathbf{x}_{t} - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ v(t) = y(t) - Z x(t) - a - D d(t)}
#' \eqn{\mathbf{x}_{t}}{x(t)} is unknown (hidden) and our data are one
#' realization of \eqn{\mathbf{y}_{t}}{y(t)}. The observed model residuals
#' \eqn{\hat{\mathbf{v}}_{t}}{hatv(t)} are the difference between the observed
#' data and the predicted data at time \eqn{t} using the fitted model.
#' `MARSSresiduals.tT` fits the model using all the data, thus
#' \deqn{ \hat{\mathbf{v}}_{t} = \mathbf{y}_{t} - \mathbf{Z}\mathbf{x}_{t}^T - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ hatv(t) = y(t) - Z xtT(t) - a - D d(t)}
#' where \eqn{\mathbf{x}_{t}^T}{xtT(t)} is the expected value of
#' \eqn{\mathbf{X}_{t}}{X(t)} conditioned on the data from 1 to \eqn{T} (all
#' the data), i.e. the Kalman smoother estimate of the states at time \eqn{t}.
#' \eqn{\mathbf{y}_{t}}{y(t)} are your data and missing values will appear as
#' NA in the observed model residuals. These are returned as `model.residuals`
#' and rows 1 to \eqn{n} of `residuals`.
#'
#' `res1` and `res2` in the code below will be the same.
#' ```
#' dat = t(harborSeal)[2:3,]
#' fit = MARSS(dat)
#' Z = coef(fit, type="matrix")$Z
#' A = coef(fit, type="matrix")$A
#' res1 = dat - Z %*% fit$states - A %*% matrix(1,1,ncol(dat))
#' res2 = MARSSresiduals(fit, type="tT")$model.residuals
#' ```
#'
#' **State residuals**
#'
#' \eqn{\mathbf{w}_{t+1}}{w(t+1)} are the difference between the state at time
#' \eqn{t+1} and the expected value of the state at time \eqn{t+1} given the
#' state at time \eqn{t}:
#' \deqn{ \mathbf{w}_{t+1} = \mathbf{x}_{t+1} - \mathbf{B} \mathbf{x}_{t} - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ w(t+1) = x(t+1) - B x(t) - u - C c(t+1)}
#' The estimated state residuals \eqn{\hat{\mathbf{w}}_{t+1}}{hatw(t+1)} are
#' the difference between estimate of \eqn{\mathbf{x}_{t+1}}{x(t+1)} minus the
#' estimate using \eqn{\mathbf{x}_{t}}{x(t)}.
#' \deqn{ \hat{\mathbf{w}}_{t+1} = \mathbf{x}_{t+1}^T - \mathbf{B}\mathbf{x}_{t}^T - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ hatw(t+1) = xtT(t+1) - B xtT(t) - u - C c(t+1)}
#' where \eqn{\mathbf{x}_{t+1}^T}{xtT(t+1)} is the Kalman smoother estimate of
#' the states at time \eqn{t+1} and \eqn{\mathbf{x}_{t}^T}{xtT(t)} is the
#' Kalman smoother estimate of the states at time \eqn{t}. The estimated state
#' residuals \eqn{\mathbf{w}_{t+1}}{w(t+1)} are returned in `state.residuals`
#' and rows \eqn{n+1} to \eqn{n+m} of `residuals`. `state.residuals[,t]` is
#' \eqn{\mathbf{w}_{t+1}}{w(t+1)} (notice time subscript difference). There
#' are no NAs in the estimated state residuals as an estimate of the state
#' exists whether or not there are associated data.
#'
#' `res1` and `res2` in the code below will be the same.
#' ```
#' dat <- t(harborSeal)[2:3,]
#' TT <- ncol(dat)
#' fit <- MARSS(dat)
#' B <- coef(fit, type="matrix")$B
#' U <- coef(fit, type="matrix")$U
#' statestp1 <- MARSSkf(fit)$xtT[,2:TT]
#' statest <- MARSSkf(fit)$xtT[,1:(TT-1)]
#' res1 <- statestp1 - B %*% statest - U %*% matrix(1,1,TT-1)
#' res2 <- MARSSresiduals(fit, type="tT")$state.residuals[,1:(TT-1)]
#' ```
#' Note that the state residual at the last time step (not shown) will be NA
#' because it is the residual associated with \eqn{\mathbf{x}_T}{x(T)} to
#' \eqn{\mathbf{x}_{T+1}}{x(T+1)} and \eqn{T+1} is beyond the data. Similarly,
#' the variance matrix at the last time step will have NAs for the same reason.
#'
#' **Variance of the residuals**
#'
#' In a state-space model, \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} are
#' stochastic, and the model and state residuals are random variables
#' \eqn{\hat{\mathbf{V}}_{t}}{hatV(t)} and
#' \eqn{\hat{\mathbf{W}}_{t+1}}{hatW(t+1)}. To evaluate the residuals we
#' observed (with \eqn{\mathbf{y}^{(1)}}{y(1)}), we use the joint distribution
#' of \eqn{\hat{\mathbf{V}}_{t}, \hat{\mathbf{W}}_{t+1}}{hatV(t), hatW(t+1)}
#' across all the different possible data sets that our MARSS equations with
#' parameters \eqn{\Theta}{Theta} might generate. Denote the matrix of
#' \eqn{\hat{\mathbf{V}}_{t}, \hat{\mathbf{W}}_{t+1}}{hatV(t), hatW(t+1)}, as
#' \eqn{\widehat{\mathcal{E}}_{t}}{Epsilon(t)}. That distribution has an
#' expected value (mean) and variance:
#' \deqn{ \textrm{E}[\widehat{\mathcal{E}}_{t}] = 0; \textrm{var}[\widehat{\mathcal{E}}_{t}] = \hat{\Sigma}_{t} }{ E[Epsilon(t)] = 0; var[Epsilon(t)] = hatSigma(t)}
#' Our observed residuals (returned in `residuals`) are one sample from this
#' distribution. To standardize the observed residuals, we will use
#' \eqn{ \hat{\Sigma}_{t} }{ hatSigma(t) }. \eqn{ \hat{\Sigma}_{t} }{ hatSigma(t) }
#' is returned in `var.residuals`. Rows/columns 1 to \eqn{n} are the
#' conditional variances of the model residuals and rows/columns \eqn{n+1} to
#' \eqn{n+m} are the conditional variances of the state residuals. The
#' off-diagonal blocks are the covariances between the two types of residuals.
#'
#' **Standardized residuals**
#'
#' `MARSSresiduals` will return the Cholesky standardized residuals sensu Harvey
#' et al. (1998) in `std.residuals` for outlier and shock detection. These are
#' the model and state residuals multiplied by the inverse of the lower triangle
#' of the Cholesky decomposition of `var.residuals` (note `chol()` in R returns
#' the upper triangle thus a transpose is needed). The standardized model
#' residuals are set to NA when there are missing data. The standardized state
#' residuals however always exist since the expected value of the states exist
#' without data. The calculation of the standardized residuals for both the
#' observations and states requires the full residuals variance matrix. Since
#' the state residuals variance is NA at the last time step, the standardized
#' residual in the last time step will be all NA (for both model and state
#' residuals).
#'
#' The interpretation of the Cholesky standardized residuals is not
#' straight-forward when the \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R}
#' variance-covariance matrices are non-diagonal. The residuals which were
#' generated by a non-diagonal variance-covariance matrices are transformed into
#' orthogonal residuals in \eqn{\textrm{MVN}(0,\mathbf{I})}{MVN(0,I)} space.
#' For example, if v is 2x2 correlated errors with variance-covariance matrix R.
#' The transformed residuals (from this function) for the i-th row of
#' \eqn{\mathbf{v}}{v} is a combination of the row 1 effect and the row 1
#' effect plus the row 2 effect. So in this case, row 2 of the transformed
#' residuals would not be regarded as solely the row 2 residual but rather how
#' different row 2 is from row 1, relative to expected. If the errors are highly
#' correlated, then the transformed residuals can look rather non-intuitive.
#'
#' The marginal standardized residuals are returned in `mar.residuals`. These
#' are the model and state residuals multiplied by the inverse of the diagonal
#' matrix formed by the square root of the diagonal of `var.residuals`. These
#' residuals will be correlated (across the residuals at time \eqn{t}) but are
#' easier to interpret when \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} are
#' non-diagonal.
#'
#' The Block Cholesky standardized residuals are like the Cholesky standardized
#' residuals except that the full variance-covariance matrix is not used, only
#' the variance-covariance matrix for the model or state residuals (respectively)
#' is used for standardization. For the model residuals, the Block Cholesky
#' standardized residuals will be the same as the Cholesky standardized
#' residuals because the upper triangle of the lower triangle of the Cholesky
#' decomposition (which is what we standardize by) is all zero. For the state
#' residuals, the Block Cholesky standardization will be different because Block
#' Cholesky standardization treats the model and state residuals as independent
#' (which they are not in the smoothations case).
#'
#' **Normalized residuals**
#'
#' If `normalize=FALSE`, the unconditional variance of \eqn{\mathbf{V}_t}{V(t)}
#' and \eqn{\mathbf{W}_t}{W(t)} are \eqn{\mathbf{R}}{R} and \eqn{\mathbf{Q}}{Q}
#' and the model is assumed to be written as
#' \deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{v}_t}{ y(t) = Z x(t) + a + v(t)}
#' \deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{w}_t}{ x(t) = B x(t-1) + u + w(t)}
#' If normalize=TRUE, the model is assumed to be written
#' \deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{H}\mathbf{v}_t}{ y(t) = Z x(t) + a + Hv(t)}
#' \deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{G}\mathbf{w}_t}{ x(t) = B x(t-1) + u + Gw(t)}
#' with the variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)}
#' equal to \eqn{\mathbf{I}}{I} (identity).
#'
#' `MARSSresiduals.tT` returns the residuals defined as in the first equations.
#' To get the residuals defined as Harvey et al. (1998) define them (second
#' equations), then use `normalize=TRUE`. In that case the unconditional
#' variance of residuals will be \eqn{\mathbf{I}}{I} instead of
#' \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R}.
#'
#' **Missing or left-out data**
#'
#' \eqn{ \textrm{E}[\widehat{\mathcal{E}}_{t}] }{ E[Epsilon(t)] } and
#' \eqn{ \textrm{var}[\widehat{\mathcal{E}}_{t}] }{ var[Epsilon(t)] } are for
#' the distribution across all possible \eqn{\mathbf{X}}{X} and
#' \eqn{\mathbf{Y}}{Y}. We can also compute the expected value and variance
#' conditioned on a specific value of \eqn{\mathbf{Y}}{Y}, the one we observed
#' \eqn{\mathbf{y}^{(1)}}{y(1)} (Holmes 2014). If there are no missing values,
#' this is not very interesting as
#' \eqn{\textrm{E}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}]=\hat{\mathbf{v}}_{t}}{E[hatV(t)|y(1)] = hatv(t)}
#' and
#' \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}] = 0}{var[hatV(t)|y(1)] = 0}.
#' If we have data that are missing because we left them out, however,
#' \eqn{\textrm{E}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}]}{E[hatV(t)|y(1)]}
#' and
#' \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}]}{var[hatV(t)|y(1)]}
#' are the values we need to evaluate whether the left-out data are unusual
#' relative to what you expect given the data you did collect.
#'
#' `E.obs.residuals` is the conditional expected value
#' \eqn{\textrm{E}[\hat{\mathbf{V}}|\mathbf{y}^{(1)}]}{E[hatV(t)|y(1)]}
#' (notice small \eqn{\mathbf{y}}{y}). It is
#' \deqn{\textrm{E}[\mathbf{Y}_{t}|\mathbf{y}^{(1)}] - \mathbf{Z}\mathbf{x}_t^T - \mathbf{a} }{ E[Y(t)|y(1)] - Z xtT(t) - a}
#' It is similar to \eqn{\hat{\mathbf{v}}_{t}}{hatv(t)}. The difference is the
#' \eqn{\mathbf{y}}{y} term.
#' \eqn{\textrm{E}[\mathbf{Y}^{(1)}_{t}|\mathbf{y}^{(1)}] }{ E[Y(1)(t)|y(1)] }
#' is \eqn{\mathbf{y}^{(1)}_{t}}{y(1)(t)} for the non-missing values. For the
#' missing values, the value depends on \eqn{\mathbf{R}}{R}. If
#' \eqn{\mathbf{R}}{R} is diagonal,
#' \eqn{\textrm{E}[\mathbf{Y}^{(2)}_{t}|\mathbf{y}^{(1)}] }{ E[Y(2)(t)|y(1)] }
#' is \eqn{\mathbf{Z}\mathbf{x}_t^T + \mathbf{a}}{Z xtT(t) + a} and the
#' expected residual value is 0. If \eqn{\mathbf{R}}{R} is non-diagonal
#' however, it will be non-zero.
#'
#' `var.obs.residuals` is the conditional variance
#' \eqn{\textrm{var}[\hat{\mathbf{V}}|\mathbf{y}^{(1)}]}{var[hatV(t)|y(1)]}
#' (eqn 24 in Holmes (2014)). For the non-missing values, this variance is 0
#' since
#' \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV(t)|y(1)} is a fixed value. For
#' the missing values,
#' \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV(t)|y(1)} is not fixed because
#' \eqn{\mathbf{Y}^{(2)}}{Y(2)} is a random variable. For these values, the
#' variance of \eqn{\hat{\mathbf{V}}|\mathbf{y}^{(1)}}{hatV(t)|y(1)} is
#' determined by the variance of \eqn{\mathbf{Y}^{(2)}}{Y(2)} conditioned on
#' \eqn{\mathbf{Y}^{(1)}=\mathbf{y}^{(1)}}{Y(1)=y(1)}. This variance matrix is
#' returned in `var.obs.residuals`. The variance of
#' \eqn{\hat{\mathbf{W}}|\mathbf{y}^{(1)}}{hatW(t)|y(1)} is 0 and thus is not
#' included.
#'
#' The variance
#' \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{Y}^{(1)}] }{ var[hatV(t)|Y(1)] }
#' (uppercase \eqn{ \mathbf{Y} }{Y}) returned in the 1 to \eqn{n}
#' rows/columns of `var.residuals` may also be of interest depending on what
#' you are investigating with regards to missing values. For example, it may be
#' of interest in a simulation study or cases where you have multiple replicated
#' \eqn{\mathbf{Y}}{Y} data sets. `var.residuals` would allow you to determine
#' if the left-out residuals are unusual with regards to what you would expect
#' for left-out data in that location of the \eqn{\mathbf{Y}}{Y} matrix but not
#' specifically relative to the data you did collect. If \eqn{\mathbf{R}}{R} is
#' non-diagonal and the \eqn{\mathbf{y}^{(1)}}{y(1)} and
#' \eqn{\mathbf{y}^{(2)}}{y(2)} are highly correlated, the variance of
#' \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{Y}^{(1)}] }{ var[hatV(t)|Y(1)] }
#' and variance of
#' \eqn{\textrm{var}[\hat{\mathbf{V}}_{t}|\mathbf{y}^{(1)}] }{ var[hatV(t)|y(1)] }
#' for the left-out data would be quite different. In the latter, the variance
#' is low because \eqn{\mathbf{y}^{(1)} }{ y(1) } has strong information about
#' \eqn{\mathbf{y}^{(2)} }{ y(2) }. In the former, we integrate over
#' \eqn{\mathbf{Y}^{(1)} }{ Y(1) } and the variance could be high (depending on
#' the parameters).
#'
#' Note, if `Harvey=TRUE` then the rows and columns of `var.residuals`
#' corresponding to missing values will be NA. This is because the Harvey et al.
#' algorithm does not compute the residual variance for missing values.
#'
#' @author Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [MARSSresiduals()], [MARSSresiduals.tt1()], [fitted.marssMLE()], [plot.marssMLE()]
#'
#' @examples
#' dat <- t(harborSeal)
#' dat <- dat[c(2, 11), ]
#' fit <- MARSS(dat)
#'
#' # state residuals
#' state.resids1 <- MARSSresiduals(fit, type = "tT")$state.residuals
#' # this is the same as hatx_t-(hatx_{t-1}+u)
#' states <- fit$states
#' state.resids2 <- states[, 2:30] - states[, 1:29] - matrix(coef(fit, type = "matrix")$U, 2, 29)
#' # compare the two
#' cbind(t(state.resids1[, -30]), t(state.resids2))
#'
#' # normalize the state residuals to a variance of 1
#' Q <- coef(fit, type = "matrix")$Q
#' state.resids1 <- MARSSresiduals(fit, type = "tT", normalize = TRUE)$state.residuals
#' state.resids2 <- (solve(t(chol(Q))) %*% state.resids2)
#' cbind(t(state.resids1[, -30]), t(state.resids2))
#'
#' # Cholesky standardized (by joint variance) model & state residuals
#' MARSSresiduals(fit, type = "tT")$std.residuals
#'
#' # Returns residuals in a data frame in long form
#' residuals(fit, type = "tT")
#'
#' @references
#' Harvey, A., S. J. Koopman, and J. Penzer. 1998. Messy time series: a unified
#' approach. Advances in Econometrics 13: 103-144 (see page 112-113). Equation
#' 21 is the Kalman eqns. Eqn 23 and 24 is the backward recursion to compute
#' the smoothations. This function uses the MARSSkf output for eqn 21 and then
#' implements the backwards recursion in equation 23 and equation 24. Pages
#' 120-134 discuss the use of standardized residuals for outlier and structural
#' break detection.
#'
#' de Jong, P. and J. Penzer. 1998. Diagnosing shocks in time series. Journal
#' of the American Statistical Association 93: 796-806. This one shows the same
#' equations; see eqn 6. This paper mentions the scaling based on the inverse of
#' the sqrt (Cholesky decomposition) of the variance-covariance matrix for the
#' residuals (model and state together). This is in the right column, half-way
#' down on page 800.
#'
#' Koopman, S. J., N. Shephard, and J. A. Doornik. 1999. Statistical algorithms
#' for models in state space using SsfPack 2.2. Econometrics Journal 2: 113-166.
#' (see pages 147-148).
#'
#' Harvey, A. and S. J. Koopman. 1992. Diagnostic checking of
#' unobserved-components time series models. Journal of Business & Economic
#' Statistics 4: 377-389.
#'
#' Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models.
#' Technical Report. arXiv:1411.0045.
#'
#' @export
MARSSresiduals.tT <- function(object, Harvey = FALSE, normalize = FALSE, silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  # These are the residuals and their variance conditioned on all the data
  # Harvey=TRUE uses Harvey et al (1998) algorithm to compute these
  # Harvey=FALSE uses the straight smoother output
  # model.residuals y(t|yT)-Zx(t|yT)-a
  # state.residuals x(t|yT)-Bx(t-1|yT)-u
  # var.residuals variance of above conditioned on y(1)
  # for missing values, Harvey=TRUE returns 0 for var for y_i missing and Harvey=FALSE returns R + Z VtT t(Z)
  # Note, I think there is a problem with the Harvey algorithm when the variance of the state residuals (Q)
  # is non-diagonal and there are missing values; it can become non-invertible

  ######################################
  # Set up variables
  MLEobj <- object
  if (missing(fun.kf)) {
    fun.kf <- MLEobj$fun.kf
  } else {
    MLEobj$fun.kf <- fun.kf
    # to ensure that MARSSkf, MARSShatyt and fitted() use the spec'd fun
  }
  if (fun.kf == "MARSSkfas" && Harvey == TRUE) stop("MARSSresiduals.tT: Harvey=TRUE requires the Kalman gain thus MARSSkfss must be used. Pass in fun.kf='MARSSkfss'.\n", call. = FALSE)

  model.dims <- attr(MLEobj$marss, "model.dims")
  TT <- model.dims[["x"]][2]
  m <- model.dims[["x"]][1]
  n <- model.dims[["y"]][1]
  y <- MLEobj$marss$data
  # set up holders
  et <- st.et <- mar.st.et <- bchol.st.et <- matrix(0, n + m, TT)
  var.et <- array(0, dim = c(n + m, n + m, TT))
  msg <- NULL

  #### list of time-varying parameters
  time.varying <- is.timevarying(MLEobj)

  kf <- MARSSkf(MLEobj)
  Ey <- MARSShatyt(MLEobj)
  Rt <- parmat(MLEobj, "R", t = 1)$R # returns matrix
  Ht <- parmat(MLEobj, "H", t = 1)$H
  Rt <- Ht %*% tcrossprod(Rt, Ht)
  Zt <- parmat(MLEobj, "Z", t = 1)$Z
  Qtp <- parmat(MLEobj, "Q", t = 2)$Q
  Gtp <- parmat(MLEobj, "G", t = 2)$G
  Qtp <- Gtp %*% tcrossprod(Qtp, Gtp)
  Btp <- parmat(MLEobj, "B", t = 2)$B
  utp <- parmat(MLEobj, "U", t = 2)$U
  ######################################

  ######################################
  # Compute residuals via Holmes algorithm
  if (!Harvey) {
    # model.et will be 0 where no data E(y)-modeled(y)
    model.et <- Ey$ytT - fitted(MLEobj, type = "ytT", output = "matrix") # model residuals
    et[1:n, ] <- model.et

    for (t in 1:TT) {
      # model residuals
      if (time.varying$R) Rt <- parmat(MLEobj, "R", t = t)$R # returns matrix
      if (time.varying$H) Ht <- parmat(MLEobj, "H", t = t)$H
      if (time.varying$R || time.varying$H) Rt <- Ht %*% tcrossprod(Rt, Ht)
      if (time.varying$Z) Zt <- parmat(MLEobj, "Z", t = t)$Z

      # model.et defined outside for loop

      # compute the variance of the residuals and state.et
      St <- Ey$yxtT[, , t] - tcrossprod(Ey$ytT[, t, drop = FALSE], kf$xtT[, t, drop = FALSE])
      tmpvar.et <- Rt - Zt %*% tcrossprod(kf$VtT[, , t], Zt) + tcrossprod(St, Zt) + tcrossprod(Zt, St)

      if (t < TT) { # fill in var.et for t (model resid)
        if (time.varying$Q) Qtp <- parmat(MLEobj, "Q", t = t + 1)$Q
        if (time.varying$G) Gtp <- parmat(MLEobj, "G", t = t + 1)$G
        if (time.varying$Q || time.varying$G) Qtp <- Gtp %*% tcrossprod(Qtp, Gtp)

        if (time.varying$B) Btp <- parmat(MLEobj, "B", t = t + 1)$B
        if (time.varying$U) utp <- parmat(MLEobj, "U", t = t + 1)$U

        Sttp <- Ey$yxttpT[, , t] - tcrossprod(Ey$ytT[, t, drop = FALSE], kf$xtT[, t + 1, drop = FALSE])
        cov.et <- tcrossprod(Zt, kf$Vtt1T[, , t + 1]) - Zt %*% tcrossprod(kf$VtT[, , t], Btp) - Sttp + tcrossprod(St, Btp)
        tmpvar.state.et <- Qtp - kf$VtT[, , t + 1] - Btp %*% tcrossprod(kf$VtT[, , t], Btp) + tcrossprod(kf$Vtt1T[, , t + 1], Btp) + tcrossprod(Btp, kf$Vtt1T[, , t + 1])

        et[(n + 1):(n + m), t] <- kf$xtT[, t + 1] - Btp %*% kf$xtT[, t] - utp
      } else {
        cov.et <- matrix(0, n, m)
        tmpvar.state.et <- matrix(0, m, m)
      }
      var.et[1:n, , t] <- cbind(tmpvar.et, cov.et)
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
        var.et[, , t] <- RQinv %*% tcrossprod(var.et[, , t], RQinv)
      }
    }
  } else {
    ######################################
    # use Harvey algorithm
    # Reference page 112-133 in Messy Time Series
    # Reference de Jong and Penzer 1998; with model transformed so sigma^2 = 1
    # refs in man file
    # NOTATION here uses that in Messy Time Series (Harvey et al 1998) with some differences
    # MARSS uses Koopman's terminology where w_t = G%*%w, Harvey uses H%*%w
    # By definition residual = 0 at time t, where x_0 is defined at t=0 or t=1, **when x_0 is estimated**
    # because x_0 is estimated and the max L will be when residual is 0
    # x_TT = Bx_{TT-1}+u_{TT-1}+w_{TT-1}, so w_TT cannot be computed (you'd need X_{TT+1}
    # so state residual at t=TT is NA
    # The state residuals are diff(states) but scaled by t(chol(Q))
    rt <- matrix(0, m, TT)
    ut <- matrix(0, n, TT)
    Nt <- array(0, dim = c(m, m, TT))
    Mt <- array(0, dim = c(n, n, TT))
    Jt <- matrix(0, m, TT)
    vt <- kf$Innov
    # If they are time-varying, Q, G and B at t=T will not appear (cancelled out by r_T and N_T = 0).
    # Set for t=1. Will update in for loop if time-varying.
    pari <- parmat(MLEobj, t = 1)
    Qtp <- pari[["Q"]]
    Gtp <- pari[["G"]]
    Ttp <- pari[["B"]]
    Z <- pari[["Z"]] # base, will be modified if missing values
    R <- pari[["R"]] # base, will be modified if missing values
    Ht <- pari[["H"]]
    for (t in seq(TT, 1, -1)) {
      # define all the, potential time-varying, parameters
      if (t < TT) {
        if (time.varying[["Q"]]) Qtp <- parmat(MLEobj, "Q", t = t + 1)$Q
        if (time.varying[["G"]]) Gtp <- parmat(MLEobj, "G", t = t + 1)$G
        if (time.varying[["B"]]) Ttp <- parmat(MLEobj, "B", t = t + 1)$B
      }
      # Zt and Rt modified in missing values case below so reset even
      # if time-varying
      if (time.varying[["Z"]]) {
        Zt <- parmat(MLEobj, "Z", t = t)$Z
      } else {
        Zt <- Z
      }
      if (time.varying[["R"]]) {
        Rt <- parmat(MLEobj, "R", t = t)$R
      } else {
        Rt <- R
      }
      if (time.varying[["H"]]) Ht <- parmat(MLEobj, "H", t = t)$H

      # implement missing values modifications per Shumway and Stoffer
      diag.Rt <- diag(Rt)
      Rt[is.na(y[, t]), ] <- 0
      Rt[, is.na(y[, t])] <- 0
      # diag(Rt)=diag.Rt
      Zt[is.na(y[, t]), ] <- 0

      # create the m x n+m and n x n+m matrices
      Rstar <- matrix(0, n, n + m)
      if (normalize) {
        Rstar[, 1:n] <- tcrossprod(Ht, pchol(Rt))
      } else {
        Rstar[, 1:n] <- Ht %*% tcrossprod(Rt, Ht)
      }
      Qpstar <- matrix(0, m, n + m)
      # MARSS uses Koopman's terminology where w_t = G%*%w, Harvey uses H%*%w
      if (normalize) {
        Qpstar[, (n + 1):(n + m)] <- tcrossprod(Gtp, pchol(Qtp))
      } else {
        Qpstar[, (n + 1):(n + m)] <- Gtp %*% tcrossprod(Qtp, Gtp)
      }

      # Don't use kf$Sigma since that has (i,i)=1
      Ftinv <- psolve(Zt %*% tcrossprod(kf$Vtt1[, , t], Zt) + Rt)

      # Harvey algorithm modified to return non-normalized errors
      Kt <- Ttp %*% matrix(kf$Kt[, , t], m, n) # R is dropping the dims so we force it to be mxn
      Lt <- Ttp - Kt %*% Zt
      Jt <- Qpstar - Kt %*% Rstar
      ut[, t] <- Ftinv %*% vt[, t, drop = FALSE] - t(Kt) %*% rt[, t, drop = FALSE]
      if (t > 1) rt[, t - 1] <- t(Zt) %*% ut[, t, drop = FALSE] + t(Ttp) %*% rt[, t, drop = FALSE]
      # Mt[,,t] = Ftinv + t(Kt)%*%Nt[,,t]%*%Kt #not used
      if (t > 1) Nt[, , t - 1] <- t(Zt) %*% Ftinv %*% Zt + t(Lt) %*% Nt[, , t] %*% Lt
      et[, t] <- t(Rstar) %*% ut[, t, drop = FALSE] + t(Qpstar) %*% rt[, t, drop = FALSE]
      # see deJong and Penzer 1998, page 800, right column, halfway down
      var.et[, , t] <- t(Rstar) %*% Ftinv %*% Rstar + t(Jt) %*% Nt[, , t] %*% Jt
    }
  }
  ######################################

  ######################################
  # prepare standardized residuals
  for (t in 1:TT) {
    tmpvar <- sub3D(var.et, t = t)
    resids <- et[, t, drop = FALSE]
    # don't includ values for resids if there is no residual (no data)
    # replace NAs with 0s
    is.miss <- c(is.na(y[, t]), rep(FALSE, m))
    resids[is.miss] <- 0

    tmpvar[abs(tmpvar) < sqrt(.Machine$double.eps)] <- 0

    # Marginal
    # inverse of diagonal of variance matrix for marginal standardization
    # psolve deals with 0s on diagonal
    tmpvarinv <- try(psolve(makediag(takediag(tmpvar))), silent = TRUE)
    if (inherits(tmpvarinv, "try-error")) {
      mar.st.et[, t] <- NA
      msg <- c(msg, paste('MARSSresiduals.tT warning: the diagonal matrix of the variance of the residuals at t =", t, "is not invertible.  NAs returned for mar.residuals at t =", t, "\n'))
    } else { # inverse of the diagonal is ok
      mar.st.et[, t] <- sqrt(tmpvarinv) %*% resids
      mar.st.et[is.miss, t] <- NA
    }

    # Block Cholesky
    # psolve and pchol deal with 0s on diagonal
    tmpchol <- try(pchol(tmpvar[(n + 1):(n + m), (n + 1):(n + m), drop = FALSE]), silent = TRUE)
    if (inherits(tmpchol, "try-error")) {
      bchol.st.et[(n + 1):(n + m), t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tT warning: the chol of the variance of the state residuals at t =", t, "returned errors.  NAs returned for bchol.std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
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
      msg <- c(msg, paste("MARSSresiduals.tT warning: the chol of the variance of the residuals at t =", t, "returned errors.  NAs returned for std.residuals at t =", t, ". See MARSSinfo(\"residvarinv\")\n"))
      next
    }
    # chol() returns the upper triangle. We need the lower triangle so t()
    tmpcholinv <- try(psolve(t(tmpchol)), silent = TRUE)
    if (inherits(tmpcholinv, "try-error")) {
      st.et[, t] <- NA
      msg <- c(msg, paste("MARSSresiduals.tT warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =", t, ". See MARSSinfo('residvarinv')\n"))
      next
    }
    st.et[, t] <- tmpcholinv %*% resids
    st.et[is.miss, t] <- NA
  }
  bchol.st.et[1:n, ] <- st.et[1:n, ] # because the upper right block of the lower tri is 0
  ######################################

  ######################################
  # NA at last time step
  # the state.residual at the last time step is NA because it is x(T+1) - f(x(T)) and T+1 does not exist.  For the same reason, the var.residuals at TT will have NAs
  et[(n + 1):(n + m), TT] <- NA
  var.et[, (n + 1):(n + m), TT] <- NA
  var.et[(n + 1):(n + m), , TT] <- NA
  st.et[, TT] <- NA
  mar.st.et[(n + 1):(n + m), TT] <- NA
  bchol.st.et[(n + 1):(n + m), TT] <- NA
  if (Harvey == TRUE) {
    # Harvey algorithm doesn't calculate var for missing data
    for (t in 1:TT) {
      is.miss <- c(is.na(y[, t]), rep(FALSE, m))
      var.et[is.miss, 1:n, t] <- NA
      var.et[1:n, is.miss, t] <- NA
    }
  }
  ######################################

  ######################################
  # Variance of missing values conditioned on the data
  # et is the expected value of the residuals conditioned on y(1)-the observed data
  E.obs.v <- et[1:n, , drop = FALSE]
  var.obs.v <- array(0, dim = c(n, n, TT))
  for (t in 1:TT) var.obs.v[, , t] <- Ey$OtT[, , t] - tcrossprod(Ey$ytT[, t])
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
  if (!is.null(msg) && object[["control"]][["trace"]] >= 0 && !silent) cat("MARSSresiduals.tT reported warnings. See msg element or attribute of returned residuals object.\n")
  ######################################

  return(list(
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
  ))
}
