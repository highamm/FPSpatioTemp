#' Extract Model Residuals from an \code{stlmfit} object.
#'
#' @param object a \code{stlmfit} object generated from the \code{\link{stlmfit}()} function.
#' @param type residual type: \code{"raw"} (the default) or \code{"normalized"}
#' @param ... further arguments passed to or from other methods.
#' @return a vector of residuals, consisting of each observed response minus the estimated mean
#' @examples
#' stlmfit_obj <- stlmfit(formula = response_na ~ x,
#'  data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#' residuals(stlmfit_obj, type = "raw")
#' residuals(stlmfit_obj, type = "normalized")
#' @export

residuals.stlmfit <- function(object, type = "raw", ...) {
  if (type == "raw") {
    return(object$resids)
  } else if (type == "normalized") {
    cholcov <- t(chol(object$Sigma_ss))
    resids <- forwardsolve(cholcov, object$resids)
    return(resids)
  } else {
    stop("type must be either raw or normalized")
  }
}