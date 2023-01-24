#' Build Correlation Matrix
#' 
#' Builds a correlation matrix from a covariance type, range, and distance matrix.
#' 
#' @param cov_type is a covariance model type (either "exponential" or "gaussian").
#' @param range is the range parameter.
#' @param dist_mat is the distance matrix.
#' @return r, a correlation matrix.
#' @examples 
#' 
#' coords <- cbind(c(1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2), c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1))
#' dist_mat <- stats::dist(coords, diag = TRUE, upper = TRUE) |> as.matrix()
#' 
#' build_r(cov_type = "exponential", range = 4, dist_mat = dist_mat)
#' build_r(cov_type = "gaussian", range = 0.1, dist_mat = dist_mat)   

#' @import stats
#' @export build_r

build_r <- function(cov_type = "exponential", range, dist_mat) {
  
  r <- switch(
    cov_type,
    exponential = exp(-(dist_mat / range)),
    gaussian = exp(-(dist_mat / range) ^ 2),
    triangular = (1 - dist_mat / range) * (dist_mat <= range),
    cosine = cos(dist_mat / range),
    spherical = (1 - (3 / 2) * (dist_mat / range) + (1 / 2) * (dist_mat / range) ^ 3) * (dist_mat <= range),
    none = diag(1, nrow = nrow(dist_mat)),
    stop("Choose exponential, gaussian, triangular, cosine, spherical, or none as the covariance structure")
  )
  
  
  return(r)
}
