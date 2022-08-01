#' Prints the fitted coefficient table of a spatiotemporal fitted spatial linear model.
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
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'        sp_de = 0.5, sp_ie = 0.5,
#'        t_de = 0.5, t_ie = 0.5,
#'        spt_de = 0.5, spt_ie = 0.5,
#'        sp_range = 4, t_range = 2)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data |>
#'  dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'   true = 1, false = 0))
#' stlmfit_obj <- stlmfit(formula = response_na ~ 1, data = samp_data,
#'  xcoord = "xcoords",
#' ycoord = "ycoords", tcoord = "times") 
#' print(stlmfit_obj)
#' @import stats
#' @export

print.stlmfit <- function(x,
                         digits = max(3L, getOption("digits") - 3L),
                         ...) {
  
  cat("\nParameters Estimates:\n")
  
  print(c(unlist(x$fixed_parms), x$cov_parms), digits = digits)
}
