#' Simulate Spatio-temporal Data.
#'
#' The purpose of \code{m2ll.spatiotemp.mL.R} is to simulate a set of spatio-temporal data for a separable model with some fixed spatial coordinates and a few fixed time points.
#'
#' @param coords_times is a data frame with x-coordinates in the first column, y-coordinates in the second column, and time points in the third column
#' @param spatvec is the vector of spatial parameters (partial sill, range)
#' @param rho is the autocorrelation parameter for an AR(1) time series.
#' @param betavec is the vector of coefficients for fixed effects.
#' @param XDesign is the design matrix for the fixed effects.
#' @return a list with (1) a data frame with the \code{coords_times} input binded together with the simulated counts, (2) the spatio-temporal covariance matrix
#' @export sim.spat.temp

sim.spat.temp <- function(coords_times, spatvec, rho, betavec, XDesign) {
  
  spatcoords <- unique(coords_times[ ,c(1, 2)])
  nspat <- nrow(spatcoords)
  
  distancemat <- as.matrix(stats::dist(spatcoords))
  
  Sigma <- spatvec[1] * exp(-distancemat / spatvec[2])
  
  times <- unique(coords_times[ ,3])
  ntime <- length(times)
  
  Sigmatime <- matrix(0, nrow = ntime, ncol = ntime)
  H <- abs(outer(times, times, "-")) 
  
  Sigmatime <- rho ^ H 
  Sigmaboth <- kronecker(Sigmatime, Sigma)
  
  ntotal <- nrow(coords_times)

  Dchol <- chol(Sigmaboth)
  mu <- XDesign %*% betavec
  
  z <- as.vector(round(mu + t(Dchol) %*%
      stats::rnorm(ntotal, mean = 0, sd = 1)))
  z[z < 0] <- 0
  count_vec <- z
  
  coords_times$counts <- count_vec
  
  return(list(df = coords_times, Sigmaboth))
}
