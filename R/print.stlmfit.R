#' Prints the fitted coefficient table of a spatio-temporal fitted linear model.
#'
#' This function uses the object that is output from \code{\link{stlmfit}()} of class \code{stlmfit}.
#'
#' @param x is an object generated from \code{\link{stlmfit}()}
#' @param digits is the number of digits to be displayed in the
#' model output
#' @param ... further arguments passed to or from other methods.
#' @return a list with \itemize{
#'   \item a table of fixed effects estimates and estimated spatiotemporal covariance estimates
#' }
#' @examples
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords",
#' ycoord = "ycoords", tcoord = "times") 
#' print(stlmfit_obj)
#' @import stats
#' @export

print.stlmfit <- function(x,
                         digits = max(3L, getOption("digits") - 3L),
                         ...) {
  
  cat("\nCoefficients:\n")
  print(x$fixed_parms)
  
  cat("\nCovariance Parameters:\n")
  print(x$cov_parms, digits = digits)
}


#' Prints the summary of a fitted spatio-temporal linear model.
#'
#' This function uses the object that is output from \code{\link{summary.stlmfit}()}.
#'
#' @param x is an summary object generated from \code{\link{summary.stlmfit}()}
#' @param digits is the number of digits to be displayed in the
#' model output
#' @param signif.stars is an option to show which predictors
#' are significant.
#' @param ... further arguments passed to or from other methods.
#' @return a list with \itemize{
#'   \item model formula
#'   \item summary statistics for the residuals.
#'   \item a table of fixed effects estimates and associated standard errors.
#'   \item estimated spatial covariance parameter estimates.
#' }
#' @examples
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords",
#' ycoord = "ycoords", tcoord = "times") 
#' print(summary(stlmfit_obj))
#' @import stats
#' @export

print.summary.stlmfit <- function(x,
                                  digits = max(3L,
                                               getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  
  
  cat("\nCall:\n", paste(deparse(x$catCall),
                         sep = "\n", collapse = "\n"),
      "\n", sep = "")
  
  cat("\nResiduals:\n")
  resQ = c(min(x$Residuals), quantile(x$Residuals,
                                      p = c(0.25, 0.5, 0.75),
                                      na.rm = TRUE), max(x$Residuals))
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)
  
  cat("\nCoefficients:\n")
  coefs <- x$FixedEffects
  colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  
  cat("\nCovariance Parameters:\n")
  print(x$CovarianceParms)
  
  cat("\nGeneralized R-squared:", x$GeneralizedR2,"\n")
}
