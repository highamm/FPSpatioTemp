#' Extract the AIC from a \code{\link{stlmfit}} object for comparing models.
#'
#' @param object a \code{\link{stlmfit}} object
#' @param ... further arguments passed to or from other methods.
#' @return The AIC value of the stlmfit object.  Here, AIC for restricted maximum likelihood (REML) is computed as 2 times the
#' negative log-likelihood plus 2 times the number of covariance model parameters. For REML, AIC should only be used to compare two models with the same fixed effects but different spatial covariance structures.
#' @examples
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times",
#'  cor_model_sp = "gaussian", cor_model_t = "exponential") 
#' AIC(stlmfit_obj)
#' @export

AIC.stlmfit <- function(object, ...) {
  
  return(object$AIC)
  
}

