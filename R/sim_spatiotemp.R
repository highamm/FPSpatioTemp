#' Simulation of Spatiotemporal Data
#'
#' The primary purpose of \code{sim_spatiotemp.R} is to simulate
#' spatiotemporal data from (for now) a linear sum-with-error model
#' that has both spatial covariance, temporal covariance, and 
#' completely independent random error.
#' 
#' @param nx is the number of x coordinates for a spatial grid
#' @param ny is the number of y coordinates for a temporal grid
#' @param ntime is the number of time points
#' @param betavec is the parameter vector for the fixed effects
#' @param XDesign is the fixed effects design matrix
#' @param sigma_parsil_spat is the spatial partial sill
#' @param range is the spatial range for an exponential covariance
#' @param sigma_nugget_spat is the spatial nugget
#' @param sigma_parsil_time is the temporal partial sill
#' @param rho is the temporal autocorrelation
#' @param sigma_nugget_time is the temporal nugget
#' @param sigma_nugget_spacetime is the independent variance parameter
#' @param seed a seed
#' @return a list with \itemize{
#'     \item a data frame \code{out_df} containing the \code{response}, spatial coordinates \code{xcoords} and \code{ycoords}, and time points \code{times}
#'     }
#' @examples 
#' sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'        sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'        sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'        sigma_nugget_spacetime = 0.5)
#' @import dplyr
#' @importFrom tidyr expand_grid
#' @export sim_spatiotemp

sim_spatiotemp <- function(nx = 10, ny = 10, ntime = 5, betavec = 0,
                           XDesign = matrix(1, nrow = nx * ny * ntime),
                           sigma_parsil_spat = 0.9, range = 5, 
                           sigma_nugget_spat = 0.1, 
                           sigma_parsil_time = 0.7, rho = 0.8,
                           sigma_nugget_time = 0.3,
                           sigma_nugget_spacetime = 0.4,
                           seed = round(runif(1, min = 1, max = 10000000))) {
  
  set.seed(seed)
  ##XDesign <- matrix(1, nrow = nx * ny * ntime)
  
  ## construct the distance matrix
  xcoords <- 1:nx; ycoords <- 1:ny
  allcoords <- expand.grid(xcoords, ycoords)
  names(allcoords) <- c("xcoords", "ycoords")
  distancemat <- as.matrix(dist(allcoords))
  
  ## calculate N, the total number of sites
  nspat <- nrow(distancemat)
  times <- 1:ntime
  ntime <- length(times)
  N <- ntime * nspat

  ## construct spatial correlation matrix
  Rs <- exp(-distancemat / range)
  
  ## build Zs, spatial random effects design matrix.
  onetime <- diag(1, nspat) %>% as.data.frame()
  Zs <- onetime %>% slice(rep(row_number(), ntime)) %>% as.matrix()
  
  ## build spatial components of overall variance
  comp_1 <- sigma_parsil_spat * Zs %*% Rs %*% t(Zs)
  comp_2 <- sigma_nugget_spat * Zs %*% t(Zs)
  
  
  ## construct temporal correlation matrix
  H <- abs(outer(times, times, "-")) 
  Rt <- rho ^ H 
  
  ## build Zt, temporal random effects design matrix
  Zt <- lapply(1:ntime, matrix, data = 0, nrow = nspat, ncol = ntime)
  for (i in 1:ntime) {
    Zt[[i]][ ,i] <- 1
  }
  Zt <- do.call(rbind, Zt)
  
  ## build temporal components of overall variance
  comp_3 <- sigma_parsil_time * Zt %*% Rt %*% t(Zt)
  comp_4 <- sigma_nugget_time * Zt %*% t(Zt)
  
  ## build independent error component of overall variance
  comp_5 <- diag(sigma_nugget_spacetime, nrow = N)
  
  Sigma <- comp_1 + comp_2 + comp_3 + comp_4 + comp_5
  
  epsilon <- t(chol(Sigma)) %*% rnorm(N)
  
  response <- as.vector(XDesign %*% betavec + epsilon)
  
  ## reminder: base R's expand.grid doesn't work with matrices or
  ## data frames
  space_time_info <- tidyr::expand_grid(times, allcoords)
  
  out_df <- dplyr::tibble(space_time_info, response)
  out_obj <- list(out_df = out_df, seed = seed)
  return(out_obj)
  
}

