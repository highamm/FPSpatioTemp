#' Build Random Effect Matrices
#' 
#' Create Z random effect matrices for spatial and temporal effects.
#' 
#' @param data_ord is the name of a data.frame or tibble where the data has already been ordered to space within time..
#' 
#' @return a list with \itemize{
#'   \item the spatial Z matrix.
#'   \item the temporal Z matrix
#'   }
#' @import stats
#' @importFrom purrr map_dfr
#' @export build_z_mats


build_z_mats <- function(data_ord) {
  
  n_sp <- max(data_ord$spindex)
  n_t <- max(data_ord$tindex)
  
  one_time_zsp <- diag(1, n_sp) |> as.data.frame()
  
  ## stack one_time_zs for the number of unique time points
  vec <- 1:n_t
  Zsp <- purrr::map_dfr(vec, ~one_time_zsp) |>
    as.matrix()
  Zsp_obs <- Zsp[data_ord$.observed, ]
  
  one_spat_zt <- diag(1, n_t) |> as.matrix()
  
  Zt <- rep(one_spat_zt, each = n_sp) |>
    matrix(nrow = n_sp * n_t, ncol = n_t)
  Zt_obs <- Zt[data_ord$.observed, ]
  
  
  Z_mats <- list(Z_sp = Zsp_obs, Z_t = Zt_obs)
  
  return(Z_mats)
}
