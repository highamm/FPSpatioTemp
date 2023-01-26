#' Build Spatio-temporal Covariance Matrix
#' 
#' Construct a spatio-temporal covariance matrix.
#' 
#' @param sp_de is the spatial dependent error variance (spatial partial sill).
#' @param sp_ie is the spatial independent error variance (spatial nugget).
#' @param t_de is the temporal dependent error variance (temporal partial sill).
#' @param t_ie is the temporal independent error variance (temporal nugget).
#' @param spt_de is the spatio-temporal dependent error variance (spatio-temporal partial sill).
#' @param spt_ie is the spatio-temporal independent error variance (spatio-temporal nugget). 
#' @param R_sp is the spatial correlation matrix
#' @param R_t is the temporal correlation matrix
#' @param Z_sp is the spatial random effects matrix
#' @param Z_t is the temporal random effects matrix
#' @param model_type is either \code{"product_sum"} (by default) or \code{"sum_with_error"}
#' @return a spatio-temporal covariance matrix.
#' @import stats
#' @export build_sigma

build_sigma <- function(sp_de = 0.9, sp_ie = 0.1,
                        t_de = 0.7, t_ie = 0.3,
                        spt_de = 0.4, spt_ie = 0.5,
                        R_sp, R_t,
                        Z_sp, Z_t,
                        model_type = "product_sum") {
  
  cov_sp <- sp_de * Z_sp %*% R_sp %*% t(Z_sp) + sp_ie * Z_sp %*% t(Z_sp)
  
  cov_t <- t_de * Z_t %*% R_t %*% t(Z_t) + t_ie * Z_t %*% t(Z_t)
  
  if (model_type == "product_sum") {
    
    cov_spt <- spt_de * (Z_sp %*% R_sp %*% t(Z_sp)) *
      (Z_t %*% R_t %*% t(Z_t)) +
      diag(spt_ie, nrow = nrow(cov_sp))
    
  } else if (model_type == "sum_with_error") {
    
    cov_spt <- diag(spt_ie, nrow = nrow(cov_sp))
    
  }
  
  sigma <- cov_sp + cov_t + cov_spt
  
  return(sigma)
}