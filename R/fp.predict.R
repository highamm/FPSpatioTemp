#' Finite Population Spatio-temporal Prediction
#'
#' \code{fp.predict} predicts the population total for the current year of sampling from sampled sites in the current year as well as sampled sites in past years.
#'
#' @param spatiotempdata is an object from the `give.nas` function with simulated counts, spatial coordinates, and time points.
#' @return A list with (1) the data frame with the predicted site-by-site predictions, (2) the true total for the current year, (3) the estimated total for the current year, and (4) the prediction variance.
#' @export fp.predict

fp.predict <- function(spatiotempdata) {
  
  ## separable models: can separate the spatial and temporal covariance 
  ## components so that they are multiplied in the errors
  ## To start, assume that we have a separable model
  
  ## grab the object from the give.nas function. This is a data frame
  ## with coordianates, true counts, observed counts, a sampling indicator,
  ## and a prediction indicator (a 1 indicates a site in the present that
  ## we want to include in our prediction for the current population total)
  
  
  return(pred_obj)
}