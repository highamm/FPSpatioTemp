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
#'       sp_de = 0.5, sp_range = 4, sp_ie = 0.5,
#'       t_de = 0.5, t_range = 0.7, t_ie = 0.5,
#'       spt_ie = 0.5)
#' @importFrom tibble tibble
#' @export sim_spatiotemp

sim_spatiotemp <- function(resp_type = "normal", 
                           nx = 10, ny = 10, ntime = 5, betavec = 0,
                           XDesign = matrix(1, nrow = nx * ny * ntime),
                           sp_de = 0.9, sp_ie = 0.1, 
                           t_de = 0.7, t_ie = 0.3,
                           spt_de = 0.3, spt_ie = 0.4,
                           sp_range = sqrt(2) / 3, t_range = 1 / 3,
                           seed = round(runif(1, min = 1,
                                              max = 1e7))) {
  
  set.seed(seed)

  xcoords <- seq(from = 0, to = 1, length.out = nx)
  ycoords <- seq(from = 0, to = 1, length.out = ny)
  times <- seq(from = 0, to = 1, length.out = ntime)
  
  allcoords <- expand.grid(xcoords, ycoords, times)
  names(allcoords) <- c("xcoords", "ycoords", "times")
  
  order_spt_obj <- order_spt(allcoords,
                             xcoord = "xcoords",
                             ycoord = "ycoords",
                             tcoord = "times")
  
  ## construct spatial correlation matrix
  ## * 3 to convert regular range to effective range
  R_sp <- build_r("exponential", range = sp_range * 3,
                dist_mat = order_spt_obj$h_sp_small)
  
  Z_mats <- build_z_mats(order_spt_obj$full_data)
  Z_sp <- Z_mats$Z_sp
  
  ## construct temporal correlation matrix

  R_t <- build_r("exponential", range = t_range * 3,
                dist_mat = order_spt_obj$h_t_small)
  Z_t <- Z_mats$Z_t
  
  Sigma <- build_sigma(sp_de = sp_de, sp_ie = sp_ie, 
                       t_de = t_de, t_ie = t_ie,
                       spt_de = spt_de, spt_ie = spt_ie,
                       model_type = "product_sum",
                       R_sp = R_sp, R_t = R_t,
                       Z_sp = Z_sp, Z_t = Z_t)
  
  epsilon <- t(chol(Sigma)) %*% rnorm(nrow(Sigma))
  
  response <- as.vector(XDesign %*% betavec + epsilon)
  
  if (resp_type == "lognormal") {
    response <- exp(response)
  }

  out_df <- tibble::tibble(allcoords, response)
  out_obj <- list(out_df = out_df, seed = seed)
  
  class(out_obj) <- "sim_spatiotemp"
  return(out_obj)
  
}

