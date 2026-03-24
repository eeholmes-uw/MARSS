## MARSSinits
## Set up inits
## These will be checked by the MLE object checker.
## Will return a par list that looks just like MLEobj par list
## Wants either a scalar (dim=NULL) or a matrix the same size as $par[[elem]] or a marssMLE object with the par element

#' Initial Values for MLE
#'
#' @description
#' Sets up generic starting values for parameters for maximum-likelihood estimation algorithms that use an iterative maximization routine needing starting values. Examples of such algorithms are the EM algorithm in [MARSSkem()] and Newton methods in [MARSSoptim()]. This is a utility function in the [MARSS-package]. It is not exported to the user. Users looking for information on specifying initial conditions should look at the help file for [MARSS()] and the User Guide section on initial conditions.
#'
#' The function assumes that the user passed in the inits list using the parameter names in whatever form was specified in the [MARSS()] call. The default is form="marxss". The [MARSSinits()] function calls MARSSinits_foo, where foo is the form specified in the [MARSS()] call. MARSSinits_foo translates the inits list in form foo into form marss.
#'
#' @param MLEobj An object of class [marssMLE].
#' @param inits A list of column vectors (matrices with one column) of the estimated values in each parameter matrix.
#'
#' @details
#' Creates an `inits` parameter list for use by iterative maximization algorithms.
#'
#' Default values for `inits` is supplied in `MARSSsettings.R`. The user can alter these and supply any of the following (m is the dim of X and n is the dim of Y in the MARSS model):
#'
#' * `elem=A,U` A numeric vector or matrix which will be constructed into `inits$elem` by the command `array(inits$elem),dim=c(n or m,1))`. If elem is fixed in the model, any `inits$elem` values will be overridden and replaced with the fixed value. Default is `array(0,dim=c(n or m,1))`.
#' * `elem=Q,R,B` A numeric vector or matrix. If length equals the length `MODELobj$fixed$elem` then `inits$elem` will be constructed by `array(inits$elem),dim=dim(MODELobj$fixed$elem))`. If length is 1 or equals dim of `Q` or dim of `R` then `inits$elem` will be constructed into a diagonal matrix by the command `diag(inits$elem)`. If elem is fixed in the model, any `inits$elem` values will be overridden and replaced with the fixed value. Default is `diag(0.05, dim of Q or R)` for `Q` and `R`. Default is `diag(1,m)` for `B`.
#' * `x0` If `inits$x0=-99`, then starting values for `x0` are estimated by a linear regression through the count data assuming `A` is all zero. This will be a poor start if `inits$A` is not 0. If `inits$x0` is a numeric vector or matrix, `inits$x0` will be constructed by the command `array(inits$x0),dim=c(m,1))`. If `x0` is fixed in the model, any `inits$x0` values will be overridden and replaced with the fixed value. Default is `inits$x0=-99`.
#' * `Z` If `Z` is fixed in the model, `inits$Z` set to the fixed value. If `Z` is not fixed, then the user must supply `inits$Z`. There is no default.
#' * `elem=V0` `V0` is never estimated, so this is never used.
#'
#' @return
#' A list with initial values for the estimated values for each parameter matrix in a MARSS model in marss form. So this will be a list with elements `B`, `U`, `Q`, `Z`, `A`, `R`, `x0`, `V0`, `G`, `H`, `L`.
#'
#' @note
#' Within the base code, a form-specific internal `MARSSinits` function is called to allow the output to vary based on form: `MARSSinits_dfa`, `MARSSinits_marss`, `MARSSinits_marxss`.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [marssMODEL], [MARSSkem()], [MARSSoptim()]
#' @export
MARSSinits <- function(MLEobj, inits = list(B = 1, U = 0, Q = 0.05, Z = 1, A = 0, R = 0.05, x0 = -99, V0 = 5, G = 0, H = 0, L = 0)) {
  MODELobj <- MLEobj[["marss"]]
  method <- MLEobj[["method"]]
  if (is.null(inits)) inits <- list()
  if (inherits(inits, "marssMLE")) {
    if (is.null(inits[["par"]])) {
      stop("Stopped in MARSSinits() because inits must have the par element if class marssMLE.\n", call. = FALSE)
    } else {
      inits <- inits$par
      if (!is.list(inits)) stop("Stopped in MARSSinits() because par element of inits$par must be a list if inits is marssMLE object.\n", call. = FALSE)
    }
  } else {
    if (!is.list(inits)) stop("Stopped in MARSSinits() because inits must be a list.\n", call. = FALSE)

    # MARSSinits needs a valid $par element and presumably the inits were passed in in form of $model not marss
    # so call form specific inits function to change those inits to form=marss inits
    inits.fun <- paste("MARSSinits_", attr(MLEobj$model, "form")[1], sep = "")
    tmp <- try(exists(inits.fun, mode = "function"), silent = TRUE)
    if (isTRUE(tmp)) {
      inits <- eval(call(inits.fun, MLEobj, inits))
    } else {
      stop(paste("Stopped in MARSSinits().  You need a MARSSinits_", attr(MLEobj$model, "form")[1], " function to tell MARSS how to interpret the inits list", sep = ""), call. = FALSE)
    }
  }

  alldefaults <- get("alldefaults", envir = pkg_globals)
  default <- alldefaults[[method]][["inits"]]
  for (elem in names(default)) {
    if (is.null(inits[[elem]])) inits[[elem]] <- default[[elem]]
  }
  y <- MODELobj$data
  m <- dim(MODELobj$fixed$x0)[1]
  n <- dim(MODELobj$data)[1]
  d <- MODELobj$free
  f <- MODELobj$fixed
  parlist <- list()

  g1 <- dim(MODELobj$fixed$G)[1] / m
  h1 <- dim(MODELobj$fixed$H)[1] / n
  l1 <- dim(MODELobj$fixed$L)[1] / m
  par.dims <- list(Z = c(n, m), A = c(n, 1), R = c(h1, h1), B = c(m, m), U = c(m, 1), Q = c(g1, g1), x0 = c(m, 1), V0 = c(l1, l1), G = c(m, g1), H = c(n, h1), L = c(m, l1))

  for (elem in names(par.dims)) {
    if (is.fixed(MODELobj$free[[elem]])) {
      parlist[[elem]] <- matrix(0, 0, 1) # always this when fixed
    } else { # not fixed
      # must be numeric
      if (!is.numeric(inits[[elem]])) {
        stop(paste("Stopped in MARSSinits(): ", elem, " inits must be numeric.", sep = ""), call. = FALSE)
      }
      # must be either length 1 or same length as the number of estimated values for elem
      if (!((is.null(dim(inits[[elem]])) & length(inits[[elem]]) == 1) | isTRUE(all.equal(dim(inits[[elem]]), c(dim(MODELobj$free[[elem]])[2], 1))))) {
        stop(paste("Stopped in MARSSinits(): ", elem, " inits must be either a scalar (dim=NULL) or the same size as the par$", elem, " element.", sep = ""), call. = FALSE)
      }
      parlist[[elem]] <- matrix(inits[[elem]], dim(MODELobj$free[[elem]])[2], 1)

      if (elem %in% c("B", "Q", "R", "V0") & is.null(dim(inits[[elem]]))) {
        # if inits is a scalar, make init a diagonal matrix
        # this is a debuging line; this should have been caught earlier
        tmp <- vec(makediag(inits[[elem]], nrow = par.dims[[elem]][1]))

        # replace any fixed elements with their fixed values
        fixed.row <- apply(d[[elem]] == 0, 1, all) # fixed over all t
        tmp[fixed.row] <- f[[elem]][fixed.row, 1, 1] # replace with fix value at time t
        # The funky colSum code sums a 3D matrix over the 3rd dim
        # I want to apply this tmp to all the variances and use an average over the d and f if they are time-varying
        # otherwise I could end up with 0s on the diagonal
        numvals <- colSums(aperm(d[[elem]], c(3, 1, 2)) != 0, dims = 1)
        delem <- colSums(aperm(d[[elem]], c(3, 1, 2)), dims = 1) / numvals
        delem[numvals == 0] <- 0
        numvals <- colSums(aperm(f[[elem]], c(3, 1, 2)) != 0, dims = 1)
        felem <- colSums(aperm(f[[elem]], c(3, 1, 2)), dims = 1) / numvals
        felem[numvals == 0] <- 0
        # use a pseudoinverse here so D's with 0 columns don't fail
        parlist[[elem]] <- pinv(t(delem) %*% delem) %*% t(delem) %*% (tmp - felem)
      } # c("Q","R","B","V0")

      if (elem == "x0") {
        dx0 <- sub3D(d$x0, t = 1)
        fx0 <- sub3D(f$x0, t = 1)
        if (identical(unname(inits$x0), -99)) { # get estimate of x0
          y1 <- y[, 1, drop = FALSE]
          # replace NAs (missing vals) with 0s
          y1[is.na(y1)] <- 0
          Zmat <- sub3D(f$Z, t = 1) + sub3D(d$Z, t = 1) %*% parlist$Z
          Zmat <- unvec(Zmat, dim = c(n, m))
          Amat <- sub3D(f$A, t = 1) + sub3D(d$A, t = 1) %*% parlist$A
          if (MODELobj$tinitx == 0) { # y=Z*(B*pi+U)+A
            Bmat <- sub3D(f$B, t = 1) + sub3D(d$B, t = 1) %*% parlist$B
            Bmat <- unvec(Bmat, dim = c(m, m))
            Umat <- sub3D(f$U, t = 1) + sub3D(d$U, t = 1) %*% parlist$U
          } else { # y=Z*pi + A
            Bmat <- diag(1, m)
            Umat <- matrix(0, m, 1)
          }
          # the following is by solving for pipi using
          # y1=Z1*(D*pipi+f)+a if tinit=1 or y1=Z1*(B(D*pipi+f)+U)+A if tinit=0
          tmp <- Zmat %*% Bmat %*% dx0
          if (is.solvable(tmp) == "underconstrained") {
            if (MODELobj$tinitx == 0) {
              stop("Stopped in MARSSinits(): Z B d_x0 is underconstrained and inits for x0 cannot be computed.  \n Pass in inits$x0 manually using inits=list(x0=...). \n This is not an error. It is simply a constraint of the method used to compute x0 inits in MARSSinits().")
            } else {
              stop("Stopped in MARSSinits(): Z d_x0 is underconstrained and inits for x0 cannot be computed.  \n Pass in inits$x0 manually using inits=list(x0=...).  \n This is not an error. It is simply a constraint of the method used to compute x0 inits in MARSSinits().")
            }
          }
          parlist$x0 <- pinv(tmp) %*% (y1 - Zmat %*% Bmat %*% fx0 - Zmat %*% Umat - Amat)
        } # inits$x0=-99 means solve for x0
      } # elem == x0
      rownames(parlist[[elem]]) <- colnames(d[[elem]])
    } # if elem not fixed
  } # for across elems

  parlist
}
