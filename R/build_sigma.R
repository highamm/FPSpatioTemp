#' Build Spatiotemporal Covariance Matrix
#' 
#' Construct a spatiotemporal covariance matrix.
#' 
#' @param parsil_sp is the spatial partial sill.
#' @param nugget_sp is the spatial nugget.
#' @param parsil_t is the temporal partial sill.
#' @param nugget_t is the temporal nugget.
#' @param parsil_spt is the spatiotemporal partial sill (for a product-sum model).
#' @param parsil_spt is the spatiotemporal nugget. 
#' 
#' 
#' @return a spatiotemporal covariance matrix.
#' @examples 
#' example_df <- tibble::tibble(tcoord = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4),
#' xcoord = c(1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2),
#' ycoord = c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1))
#' data_unordered <- dplyr::sample_n(example_df, size = nrow(example_df))
#' data_ord <- order_spt(data = data_unordered, xcoord = xcoord, ycoord = ycoord, tcoord = tcoord)$full_data
#'
#' z_mats <- build_z_mats(data_ord = data_ord)
#' example_r_sp <- build_r(cov_type = "exponential",
#'  range = 2, dist_mat = h_sp_small) 
#' example_r_t <- build_r(cov_type = "exponential", range = 2, dist_mat = h_t_small)

#' build_sigma(model_type = "product_sum", R_sp = example_r_sp, R_t = example_r_t,
#'            Z_sp = z_mats$Z_sp, Z_t = z_mats$Z_t)
#' @import stats
#' @export build_sigma

build_sigma <- function(parsil_sp = 0.9, nugget_sp = 0.1,
                        parsil_t = 0.7, nugget_t = 0.3,
                        parsil_spt = 0.4, nugget_spt = 0.5,
                        model_type = "product_sum",
                        R_sp, R_t,
                        Z_sp, Z_t) {
  
  cov_sp <- parsil_sp * Z_sp %*% R_sp %*% t(Z_sp) + nugget_sp * Z_sp %*% t(Z_sp)
  
  cov_t <- parsil_t * Z_t %*% R_t %*% t(Z_t) + nugget_t * Z_t %*% t(Z_t)
  
  if (model_type == "product_sum") {
    
    cov_spt <- parsil_spt * (Z_sp %*% R_sp %*% t(Z_sp)) *
      (Z_t %*% R_t %*% t(Z_t)) +
      diag(nugget_spt, nrow = nrow(cov_sp))
    
  } else if (model_type == "sum_with_error") {
    
    cov_spt <- diag(nugget_spt, nrow = nrow(cov_sp))
    
  }
  
  sigma <- cov_sp + cov_t + cov_spt
  
  return(sigma)
}