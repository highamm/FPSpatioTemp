#' Build Spatiotemporal Covariance Matrix
#' 
#' Construct a spatiotemporal covariance matrix.
#' 
#' @param parsil_sp is the spatial dependent error variance (spatial partial sill).
#' @param nugget_sp is the spatial independent error variance (spatial nugget).
#' @param parsil_t is the temporal dependent error variance (temporal partial sill).
#' @param nugget_t is the temporal independent error variance (temporal nugget).
#' @param parsil_spt is the spatio-temporal dependent error variance (spatio-temporal partial sill).
#' @param nugget_spt is the spatio-temporal independent error variance (spatio-temporal nugget). 
#' @param R_sp is the spatial correlation matrix
#' @param R_t is the temporal correlation matrix
#' @param Z_sp is the spatial random effects matrix
#' @param Z_t is the temporal random effects matrix
#' @param model_type is either \code{"product_sum"} (by default) or \code{"sum_with_error"}
#' @return a spatio-temporal covariance matrix.
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