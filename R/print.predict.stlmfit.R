#' Prints a short summary for the \code{\link{predict.stlmfit}()} function.
#'
#' This function uses the object that is output from \code{\link{predict.stlmfit}()} of class \code{predict.stlmfit}.
#' @param x is a prediction object generated from \code{\link{predict.stlmfit}()}
#' @param digits is the number of digits to be displayed in the
#' model output
#' @param ... further arguments passed to or from other methods.
#' @examples 
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'       sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'       sigma_nugget_spacetime = 0.5)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data |> 
#' dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'  true = 1, false = 0))
#' stlmfit_obj <- stlmfit(formula = response_na ~ 1, data = samp_data,
#'  xcoordcol = "xcoords",
#' ycoordcol = "ycoords", tcol = "times") 
#' pred_obj <- predict(object = stlmfit_obj, wtscol = "predwts")
#' print(pred_obj)
#' @import stats
#' @export

print.predict.stlmfit <- function(x, digits =
                                   max(3L, getOption("digits") - 3L),
                                 ...) {
  
  if (inherits(x, "predict.stlmfit") == FALSE) {
    stop("x must be of class predict.stlmfit")
  }
  
  totalpred <- x$totalpred
  predvar <- x$predvar
  pred_level <- x$pred_level
  
  
  simptab <- t(matrix(c(totalpred, sqrt(predvar))))
  colnames(simptab) <- c("Prediction", "SE")
  
  confbounds <- matrix(c(as.numeric(totalpred) + c(1, -1) *
                           stats::qnorm((1 - pred_level) / 2) *
                           sqrt(as.numeric(predvar))),
                       nrow = 1)
  labs <- c(paste(pred_level * 100, "% LB", sep = ""),
            paste(pred_level * 100,
                  "% UB", sep = ""))
  colnames(confbounds) <- labs
  
  basic_tab <- cbind(simptab, confbounds)
  rownames(basic_tab) <- base::all.vars(formula)[1]
  
  cat("Prediction Info:\n")
  print(basic_tab, digits = digits)
  
  
}
