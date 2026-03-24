#' Names for marssMLE Object Components
#'
#' @description
#' Puts names on the par, start, par.se, init components of [marssMLE]
#' objects. This is a utility function in the **MARSS** package and is not
#' exported.
#'
#' @param MLEobj An object of class [marssMLE].
#'
#' @details
#' The X.names and Y.names are attributes of [marssMODEL] objects (which
#' would be in `$marss` and `$model` in the [marssMLE] object). These names
#' are applied to the par elements in the [marssMLE] object.
#'
#' @return
#' The object passed in, with row and column names on matrices as specified.
#'
#' @author
#' Eli Holmes, NOAA, Seattle, USA.
#'
#' @seealso [marssMLE], [marssMODEL]
#'
#' @keywords internal
MARSSapplynames <- function(MLEobj) {
  ## Helper function to put names on the elements in a marssMLE object
  if (!inherits(MLEobj, "marssMLE")) {
    stop("Stopped in MARSSapplynames() because this function is for marssMLE objects only.\n", call. = FALSE)
  }

  MODELobj <- MLEobj[["marss"]]
  par.names <- attr(MODELobj, "par.names")
  X.names <- attr(MODELobj, "X.names")
  Y.names <- attr(MODELobj, "Y.names")

  # The par element has rownames that come from the column names of the free matrix
  for (elem in par.names) {
    if (!is.null(MLEobj[["par"]][[elem]]) & is.null(rownames(MLEobj[["par"]][[elem]]))) rownames(MLEobj[["par"]][[elem]]) <- colnames(MODELobj$free[[elem]])
  }

  rownames(MLEobj[["model"]][["data"]]) <- Y.names
  rownames(MLEobj[["marss"]][["data"]]) <- Y.names
  if (!is.null(MLEobj[["Ey"]][["ytT"]])) rownames(MLEobj[["Ey"]][["ytT"]]) <- Y.names
  if (!is.null(MLEobj[["kf"]][["xtT"]])) rownames(MLEobj[["kf"]][["xtT"]]) <- X.names
  if (!is.null(MLEobj[["kf"]][["xtt1"]])) rownames(MLEobj[["kf"]][["xtt1"]]) <- X.names
  if (!is.null(MLEobj[["kf"]][["xtt"]]) && !is.character(MLEobj[["kf"]][["xtt"]])) rownames(MLEobj[["kf"]][["xtt"]]) <- X.names
  if (!is.null(MLEobj[["kf"]][["x0T"]])) rownames(MLEobj[["kf"]][["x0T"]]) <- X.names
  if (!is.null(MLEobj[["kf"]][["Innov"]]) && !is.character(MLEobj[["kf"]][["Innov"]])) rownames(MLEobj[["kf"]][["Innov"]]) <- Y.names
  if (!is.null(MLEobj[["states.se"]])) rownames(MLEobj[["states.se"]]) <- X.names
  if (!is.null(MLEobj[["states"]])) rownames(MLEobj[["states"]]) <- X.names

  return(MLEobj)
}
