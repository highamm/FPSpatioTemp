#' Simulation of Spatiotemporal Data
#'
#' The primary purpose of \code{sim_spatiotemp.R} is to simulate
#' spatiotemporal data from a product-sum model.
#' 
#' @param resp_type is the distribution of the response variable (either
#' \code{"normal"} or \code{"lognormal"})
#' @param nx is the number of x coordinates for a spatial grid on the unit square
#' @param ny is the number of y coordinates for a spatial grid on the unit square
#' @param ntime is the number of equally spaced time points on the [0, 1] interval.
#' @param betavec is the parameter vector for the fixed effects
#' @param XDesign is the fixed effects design matrix
#' @param sp_de is the spatial partial sill
#' @param sp_ie is the spatial nugget
#' @param t_de is the temporal partial sill
#' @param t_ie is the temporal nugget
#' @param spt_de is the spatiotemporal partial sill
#' @param spt_ie is the independent variance parameter
#' @param sp_range is the spatial range for an exponential covariance
#' @param t_range is the temporal range for exponential covariance
#' @param seed a seed
#' @return a list with \itemize{
#'     \item a data frame \code{out_df} containing the \code{response}, spatial coordinates \code{xcoords} and \code{ycoords}, and time points \code{times}
#'     }
#' @examples 
#' sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'        sp_de = 0.5, sp_ie = 0.5,
#'        t_de = 0.5, t_ie = 0.5,
#'        spt_de = 0.5, spt_ie = 0.5,
#'        sp_range = 4, t_range = 2)
#' @importFrom dplyr slice row_number
#' @importFrom tidyr expand_grid
#' @importFrom tibble tibble
#' @export sim_spatiotemp

sim_spatiotemp <- function(resp_type = "normal", 
                           nx = 10, ny = 10, ntime = 5, betavec = 0,
                           XDesign = matrix(1, nrow = nx * ny * ntime),
                           sp_de = 0.9, 
                           sp_ie = 0.1, 
                           t_de = 0.7,
                           t_ie = 0.3,
                           spt_de = 0.3,
                           spt_ie = 0.4,
                           sp_range = sqrt(2) / 3, t_range = 1 / 3,
                           seed = round(runif(1, min = 1, max = 10000000))) {
  
  set.seed(seed)
  ##XDesign <- matrix(1, nrow = nx * ny * ntime)
  
  ## construct the distance matrix
  xcoords <- seq(from = 0, to = 1, length.out = nx)
  ycoords <- seq(from = 0, to = 1, length.out = ny)
  
  allcoords <- expand.grid(xcoords, ycoords)
  names(allcoords) <- c("xcoords", "ycoords")
  distancemat <- dist(allcoords) |> as.matrix()
  
  ## calculate N, the total number of sites
  nspat <- nrow(distancemat)
  times <- seq(from = 0, to = 1, length.out = ntime)
  ntime <- length(times)
  N <- ntime * nspat

  ## construct spatial correlation matrix
  Rs <- build_r("exponential", range = sp_range, dist_mat = distancemat)
  
  ## build Zs, spatial random effects design matrix.
  onetime <- diag(1, nspat) |> as.data.frame()
  Zs <- onetime |> dplyr::slice(rep(dplyr::row_number(), ntime)) |>
    as.matrix()
  
  ## build spatial components of overall variance
  comp_1 <- sp_de * Zs %*% Rs %*% t(Zs)
  comp_2 <- sp_ie * Zs %*% t(Zs)
  
  
  ## construct temporal correlation matrix
  dist_time <- dist(times, upper = TRUE, diag = TRUE) |> as.matrix()
  
  Rt <- build_r("exponential", range = t_range, dist_mat = dist_time)
  
  
  ## build Zt, temporal random effects design matrix
  Zt <- lapply(1:ntime, matrix, data = 0, nrow = nspat, ncol = ntime)
  for (i in 1:ntime) {
    Zt[[i]][ ,i] <- 1
  }
  Zt <- do.call(rbind, Zt)
  
  ## build temporal components of overall variance
  comp_3 <- t_de * Zt %*% Rt %*% t(Zt)
  comp_4 <- t_ie * Zt %*% t(Zt)
  
  comp_5 <- spt_de * (Zs %*% Rs %*% t(Zs)) *
    (Zt %*% Rt %*% t(Zt))
  comp_6 <- diag(spt_ie, nrow = N)
  
  Sigma <- comp_1 + comp_2 + comp_3 + comp_4 + comp_5 + comp_6
  
  epsilon <- t(chol(Sigma)) %*% rnorm(N)
  
  response <- as.vector(XDesign %*% betavec + epsilon)
  
  if (resp_type == "lognormal") {
    response <- exp(response)
  }
  ## reminder: base R's expand.grid doesn't work with matrices or
  ## data frames
  space_time_info <- tidyr::expand_grid(times, allcoords)
  
  out_df <- tibble::tibble(space_time_info, response)
  out_obj <- list(out_df = out_df, seed = seed)
  
  class(out_obj) <- "sim_spatiotemp"
  return(out_obj)
  
}

