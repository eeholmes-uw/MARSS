# NAMESPACE bootstrap for the incremental roxygen2 migration.
#
# This file mirrors the hand-written NAMESPACE so that devtools::document()
# can be safely run at any point during the migration without breaking the
# package.
#
# HOW TO MAINTAIN during the migration:
#   - When an R file gains @export for a function, remove that function's
#     @rawNamespace export() line from the "Exported functions" section below.
#   - When an R file gains @method generic class, remove the matching
#     @rawNamespace S3method() line from the "S3 method registrations" section.
#   - When a function's R file gains @importFrom, remove the matching line
#     from the "Package-level imports" section (or leave it — duplicates are
#     harmless but noisy).
#   - Delete this file entirely when all R files have been converted.
#
# NOTE: export(MARSScv) is intentionally absent — it is already produced by
#       the @export tag in R/MARSScv.R.


# ---- Exported functions -------------------------------------------------------
# Remove each line below once the corresponding R file gains an @export tag.
#
#' @rawNamespace export(accuracy)
#' @rawNamespace export(autoplot.marssMLE)
#' @rawNamespace export(autoplot.marssPredict)
#' @rawNamespace export(CSEGriskfigure)
#' @rawNamespace export(CSEGtmufigure)
#' @rawNamespace export(forecast)
#' @rawNamespace export(glance)
#' @rawNamespace export(MARSS)
#' @rawNamespace export(MARSSaic)
#' @rawNamespace export(MARSSboot)
#' @rawNamespace export(MARSShessian)
#' @rawNamespace export(MARSSinfo)
#' @rawNamespace export(MARSSinits)
#' @rawNamespace export(MARSSinnovationsboot)
#' @rawNamespace export(MARSSkem)
#' @rawNamespace export(MARSSkemcheck)
#' @rawNamespace export(MARSSkf)
#' @rawNamespace export(MARSShatyt)
#' @rawNamespace export(MARSSkfss)
#' @rawNamespace export(MARSSkfas)
#' @rawNamespace export(MARSSoptim)
#' @rawNamespace export(MARSSparamCIs)
#' @rawNamespace export(MARSSresiduals)
#' @rawNamespace export(MARSSsimulate)
#' @rawNamespace export(MARSSFisherI)
#' @rawNamespace export(MARSSvectorizeparam)
#' @rawNamespace export(tidy)
#' @rawNamespace export(zscore)
#' @rawNamespace export(ldiag)
#' @rawNamespace export(MARSSfit)
NULL

# ---- Package-level imports ----------------------------------------------------
# These broad imports are intentional (stats, utils, graphics are base packages
# that a user could theoretically detach, breaking MARSS).  The importFrom lines
# are for specific external-package functions actually used.
# Remove / move these lines when the corresponding @import / @importFrom tags
# are added to the appropriate R source files.
#
#' @rawNamespace import(stats)
#' @rawNamespace import(utils)
#' @rawNamespace import(graphics)
#' @rawNamespace importFrom(mvtnorm, rmvnorm)
#' @rawNamespace importFrom(nlme, fdHess)
#' @rawNamespace importFrom(KFAS, SSModel, SSMcustom, KFS)
#' @rawNamespace importFrom("grDevices", "contourLines")
#' @rawNamespace importFrom(generics,forecast)
#' @rawNamespace importFrom(generics,accuracy)
#' @rawNamespace importFrom(generics,glance)
#' @rawNamespace importFrom(generics,tidy)
NULL

# ---- S3 method registrations --------------------------------------------------
# Remove each line once the corresponding R file gains an @method tag.
#
#' @rawNamespace S3method(MARSSfit, default)
#' @rawNamespace S3method(MARSSfit, kem)
#' @rawNamespace S3method(MARSSfit, BFGS)
#' @rawNamespace S3method(accuracy, marssMLE)
#' @rawNamespace S3method(accuracy, marssPredict)
#' @rawNamespace S3method(coef, marssMLE)
#' @rawNamespace S3method(fitted, marssMLE)
#' @rawNamespace S3method(forecast, marssMLE)
#' @rawNamespace S3method(glance, marssMLE)
#' @rawNamespace S3method(model.frame, marssMODEL)
#' @rawNamespace S3method(model.frame, marssMLE)
#' @rawNamespace S3method(print, marssMODEL)
#' @rawNamespace S3method(print, marssMLE)
#' @rawNamespace S3method(print, marssPredict)
#' @rawNamespace S3method(plot, marssMLE)
#' @rawNamespace S3method(plot, marssPredict)
#' @rawNamespace S3method(plot, marssResiduals)
#' @rawNamespace S3method(stats::predict, marssMLE)
#' @rawNamespace S3method(logLik, marssMLE)
#' @rawNamespace S3method(residuals, marssMLE)
#' @rawNamespace S3method(simulate, marssMLE)
#' @rawNamespace S3method(summary, marssMODEL)
#' @rawNamespace S3method(summary, marssMLE)
#' @rawNamespace S3method(tidy, marssMLE)
#' @rawNamespace S3method(toLatex, marssMODEL)
#' @rawNamespace S3method(toLatex, marssMLE)
#' @rawNamespace S3method(stats::tsSmooth, marssMLE)
NULL

# ---- Conditional S3 method registrations (ggplot2) ----------------------------
# Remove once the autoplot methods gain @method tags in their R files.
#
#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(ggplot2::autoplot, marssMLE)
#'   S3method(ggplot2::autoplot, marssPredict)
#'   S3method(ggplot2::autoplot, marssResiduals)
#' }
NULL
