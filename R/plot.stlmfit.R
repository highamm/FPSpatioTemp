#' Plot a histogram of the normalized residuals and a normalized residuals vs. fitted values scatterplot
#'
#' Provides a quick way to assess model assumptions of normality and constant variance. Note that, when predicting a total, normality is only required when the number of sites is small.
#'
#' @param x is an object of class \code{\link{stlmfit}}.
#' @param ... further arguments passed to or from other methods.
#' @return plots of normalized residuals and normalized residuals vs. fitted values
#' @import stats
#' @import ggplot2
#' @examples
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times")
#' plot(stlmfit_obj)
#' @export

plot.stlmfit <- function(x, ...) {
  
  if (inherits(x, "stlmfit") == FALSE) {
    stop("x must be of class `stlmfit`")
  }
  
  norm_res <- residuals(x, type = "normalized")

  resid_df <- data.frame(resids = norm_res, fitted = x$fitted_samp)
  
  p1 <- ggplot(data = resid_df, aes(x = .data$resids)) +
    geom_histogram(colour = "black", fill = "white", 
                   bins = 20) +
    labs(x = "Normalized Residuals")
  
  p2 <- ggplot(data = resid_df, aes(x = .data$yo, y = .data$resids)) +
    geom_point() +
    labs(y = "Normalized Residuals",
         x = "Fitted Values")
  
  gridExtra::grid.arrange(p1, p2)
  
}