#' Summarizes a fitted spatio-temporal linear model.
#'
#' In conjunction with \code{print.summary.slmfit()}, the output looks similar
#' to output from \code{R}'s standard \code{lm()} function.
#'
#' @param object is an object generated from \code{\link{stlmfit}()} of class \code{stlmfit}.
#' @param ... further arguments passed to or from other methods.
#' @return a list with \itemize{
#'   \item model formula
#'   \item a table of fixed effects estimates and associated standard errors
#'   \item estimated spatial covariance parameter estimates
#'   \item residuals
#'        }
#' @examples
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#' summary(stlmfit_obj)
#' @import stats
#' @export

summary.stlmfit <- function(object, ...) {
  
  catcall <- object$summary_stlmm$Call
  
  fixed_effects_tab <- object$summary_stlmm$FixedEffects
  
  cov_parms <- object$summary_stlmm$CovarianceParameters[c("s_de", "s_ie",
                                              "s_range", 
                                              "t_de", "t_ie",
                                              "t_range",
                                              "st_de", "st_ie")]

  names(cov_parms) <- c("sp_de", "sp_ie", "sp_range",
                        "t_de", "t_ie", "t_range", 
                        "spt_de", "spt_ie")
  
  
  covariance_forms <- object$summary_stlmm$CovarianceForms[ ,2:3]
  names(covariance_forms) <- c("sp_cor", "t_cor")
  
  resid_vec <- object$summary_stlmm$Residuals

  
  outpt <- list(catcall = catcall,
                fixed_effects_tab = fixed_effects_tab,
                cov_parms = cov_parms,
                covariance_forms = covariance_forms,
                resid_vec = resid_vec)

  class(outpt) <- "summary.stlmfit"
  return(outpt)
  
}
