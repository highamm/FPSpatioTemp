#' Select Observations that Have Missing Values
#'
#' The purpose of \code{give.nas} select certain values in the data frame output of \code{sim.spat.temp} to be missing.
#'
#' @param sim_obj is an object from \code{sim.spat.temp}
#' @param nsamp is the number of sites to be sampled 
#' @param method is the how to select which sites will be sampled (for now, the only \code{method} is "Random"). In the future, method might be different if we know that some sites tend to get sampled over and over again through time.
#' @return a data frame with the spatial coordinates, time points, counts, counts with missing values. an indicator for which sites were sampled (1 = sampled, 0 = unsampled), and an indicator for which sites are in the current smapling year (1 = current, 0 = past).
#' @export give.nas



give.nas <- function(sim_obj, nsamp, method) {
  index <- sample(1:nrow(sim_obj), size = nsamp, replace = FALSE)
  sim_obj$obscounts <- sim_obj$counts
  sim_obj$obscounts[-index] <- NA
  sim_obj$sampind <- as.integer(is.na(sim_obj$obscounts) == FALSE)
  
  pred_index <- sim_obj$ts == max(sim_obj$ts)
  sim_obj$predind <- as.integer(pred_index)
  
  return(sim_obj)

}
