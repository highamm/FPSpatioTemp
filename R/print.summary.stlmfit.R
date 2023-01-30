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
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#' print(summary(stlmfit_obj))
#' @import stats
#' @export

print.summary.stlmfit <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
  
  
  cat("\nCall:\n", paste(deparse(x$catcall),
                         sep = "\n", collapse = "\n"),
      "\n", sep = "")
  
  cat("\nResiduals:\n")
  resQ = c(min(x$resid_vec), quantile(x$resid_vec,
                                      p = c(0.25, 0.5, 0.75),
                                      na.rm = TRUE), max(x$resid_vec))
  
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)
  
  cat("\nCoefficients:\n")
  coefs <- x$fixed_effects_tab
  colnames(coefs) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  
  cat("\nCovariance Parameters:\n")
  print(x$cov_parms)
  
  cat("\nCorrelation Functions:\n")
  print(x$covariance_forms)
  
}

