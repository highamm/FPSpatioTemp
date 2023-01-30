#' Spatiotemporal Model Fit
#' 
#' Estimates regression coefficients, spatial autocorrelation
#' parameters, and temporal autocorrelation parameters, given spatial
#' coordinates, time points and a model formula.
#' 
#' @param formula is an \code{R} linear model formula specifying the
#' response variable as well as covariates for predicting the response on the unsampled sites.
#' @param data is a data frame or tibble with the response column, the covariates to be used for the block kriging, the spatio-temporal coordinates for all of the sites.
#' @param xcoord is the name of the column in the data frame with x coordinates or longitudinal coordinates.
#' @param ycoord is the name of the column in the data frame with y coordinates or latitudinal coordinates.
#' @param tcoord is the name of the column in the data frame with 
#' the time points.
#' @param cor_model_sp a correlation model for the spatial correlation. Options are \code{"exponential"}, \code{"gaussian"}, \code{"spherical"}, or \code{"tent"}.
#' @param cor_model_t a correlation model for the temporal correlation Options are \code{"exponential"}, \code{"gaussian"}, \code{"spherical"}, or \code{"tent"}.
#' @return a list with \itemize{
#'   \item the original data set, appended with a few variables
#'   \item the estimated spatiotemporal covariance matrix
#'   \item the design matrix for the fixed effects
#'   \item a vector of the sampled response values
#'   \item a list with the 8 spatio-temporal parameter
#'    estimates for the fixed effects, 
#'    the Exponential spatial covariance and AR(1) temporal
#'    covariance and a sum-with-error linear mixed model
#'     (vector of fixed effects parameter estimates, spatial
#'    partial sill, range, spatial nugget, temporal partial sill,
#'    autocorrelation parameter, temporal nugget, and spatiotemporal
#'    nugget).
#'  }
#' @examples 
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#' @import stats
#' @export stlmfit
             
# stlmfit(formula = response_na ~ 1, data = samp_data,
# xcoord = "xcoords", ycoord = "ycoords", tcoord = "times")
                     
stlmfit <- function(formula, data, xcoord, ycoord, tcoord,
                    cor_model_sp = "exponential",
                    cor_model_t = "exponential") {

  ## order the data so that sites are in the same order within each time
  ## need to mess with this some more for NSE
  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)
  tcoord <- substitute(tcoord)
  

  order_spt_obj <- order_spt(data, xcoord, ycoord, tcoord)
  
  ## rows with missing coordinates (space or time) have been removed
  data_ordered <- order_spt_obj$full_data
  
  ## data_obs contains rows in the original sampling frame of data.
  data_og <- data_ordered |> dplyr::filter(.data$.observed == TRUE)
  
  
  ## make data frame of only predictors
  data_preds_only <- data.frame(data_og[ ,all.vars(formula)[-1]])
  colnames(data_preds_only) <- all.vars(formula)[-1]
  
  
  ## design matrix for all sites
  X_all <- model.matrix(formula, model.frame(formula, data_og,
                                             na.action = stats::na.pass))
  
  missing_ind <- base::apply(is.na(X_all), MARGIN = 1, FUN = sum)
  
  X_all_miss_cov <- X_all[missing_ind == 0, ]
  
  n_missing <- sum(missing_ind >= 1)
  
  ## data_obs has all observations except for....
  ## missing covariates and has already removed the "extra" 
  ## rows from order_spt()
  
  data_obs <- data_og[missing_ind == 0, ]
  
  full_mf <- stats::model.frame(formula, na.action =
                                 stats::na.pass, data = data_obs)
  
  resp_var <- stats::model.response(full_mf, "numeric")
  resp_density <- resp_var
  
  ind_sa <- !is.na(resp_var)
  data_obs$ind_sa <- ind_sa
  ind_un <- is.na(resp_var)
  data_obs$ind_un <- ind_un
  
  resp_density_samp <- resp_density[ind_sa]
  
  
  ## create Zs and Zt matrices
  ## build_z_mats filters out the rows where
  ## .observed == FALSE at the end
  z_mats <- build_z_mats(data_ord = data_ordered)
  
  Z_sp <- z_mats$Z_sp[missing_ind == 0, ] ## keep only rows where covariates were observed
  Z_t <- z_mats$Z_t[missing_ind == 0, ]
  
  total_var <- var(resp_density_samp)
  
  ## replace with grid search eventually
  sp_de_initial <- total_var / 10
  sp_ie_initial <- total_var / 10
  t_de_initial <- total_var / 10
  t_ie_initial <- total_var / 10
  spt_de_initial <- total_var / 10
  spt_ie_initial <- total_var / 10
  
  total_var_initial <- sum(
    sp_de_initial, sp_ie_initial, t_de_initial,
    t_ie_initial, spt_de_initial, spt_ie_initial
  )

  
  h_sp_small <- order_spt_obj$h_sp_small
  h_t_small <- order_spt_obj$h_t_small
  
  max_dist_sp <- max(h_sp_small, na.rm = TRUE)
  min_dist_sp <- min(h_sp_small[h_sp_small != 0], na.rm = TRUE)
  
  max_dist_t <- max(h_t_small, na.rm = TRUE)
  min_dist_t <- min(h_t_small[h_t_small != 0], na.rm = TRUE)
  
  minimax_vec <- c(max_dist_sp = max_dist_sp,
                   min_dist_sp = min_dist_sp,
                   max_dist_t = max_dist_t,
                   min_dist_t = min_dist_t)
  
  sp_range_initial <- median(h_sp_small)
  t_range_initial <- median(h_t_small)
  
  sim_initial <-  make_covparam_object(
    s_de = sp_de_initial, s_ie = sp_ie_initial,
    t_ie = t_ie_initial, t_de = t_de_initial,
    st_de = spt_de_initial, st_ie = spt_ie_initial,
    s_range = sp_range_initial, t_range = t_range_initial,
    stcov = "productsum", estmethod = "reml")
  
  data_sa <- data_obs[ind_sa, ]
  data_un <- data_obs[ind_un, ]
  
  # data_test <- data |> mutate(resp2 = rnorm(12, 0, 1))

  ## CANNOT HAVE A COLUMN NAMED index, spindex, or tindex
  fast_est <- stlmm(
    formula = formula,
    data = data_sa |> dplyr::select(-dplyr::any_of(c("`_index`", "`_spindex`", "`_tindex`"))),
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    stcov = "productsum",
    estmethod = "reml",
    s_cor = cor_model_sp,
    t_cor = cor_model_t,
    initial = sim_initial,
    condition = 1e-4
  )
  
  summary_output <- summary(fast_est)
  
  m2ll <- fast_est$Objective[1]

  parm_est <- c(unclass(fast_est$CovarianceParameters), use.names = FALSE)
  
  fitted_samp <- fast_est$model[["FixedDesignMatrix"]] %*% fast_est$Coefficients
  AIC <- m2ll + 2 * length(parm_est)
  
  resid_samp <- fast_est$Residuals
  
  sp_de_hat <- parm_est[1]
  sp_range_hat <- parm_est[7] 
  sp_ie_hat <- parm_est[2]
  
  t_de_hat <- parm_est[3]
  t_range_hat <- parm_est[8]
  t_ie_hat  <- parm_est[4]
  
  spt_de_hat <- parm_est[5]
  spt_ie_hat <- parm_est[6]

  ## overwrite with correct estimates
   parm_est[7] <- sp_range_hat
   parm_est[8] <- t_range_hat
  
  cov_parms <- c(sp_de = sp_de_hat, sp_ie = sp_ie_hat, sp_range = sp_range_hat,
                 t_de = t_de_hat, t_ie = t_ie_hat, t_range = t_range_hat,
                 spt_de = spt_de_hat, spt_ie = spt_ie_hat)
  
  R_sp_hat <- build_r(cov_type = cor_model_sp, range = sp_range_hat,
                      dist_mat = h_sp_small)
  
  R_t_hat <- build_r(cov_type = cor_model_t, range = t_range_hat,
                     dist_mat = h_t_small)

  Sigma_hat <- build_sigma(sp_de = sp_de_hat, sp_ie = sp_ie_hat,
                           t_de = t_de_hat, t_ie = t_ie_hat,
                           spt_de = spt_de_hat, spt_ie = spt_ie_hat,
                           model_type = "product_sum",
                           R_sp = R_sp_hat, R_t = R_t_hat,
                           Z_sp = Z_sp, Z_t = Z_t)
    

  Sigma_ss <- Sigma_hat[ind_sa, ind_sa, drop = FALSE]
  Sigma_ssi <- solve(Sigma_ss) ## can output this so it doesn't
  ## also get inverted in predict.stlmfit()
  
  beta_hat <- as.vector(fast_est$Coefficients)
  names(beta_hat) <- c("Intercept", all.vars(formula)[-1])
  
  stlmfit_obj <- list(data = data_obs, Sigma_hat = Sigma_hat, X_all = X_all, 
                      resp_density = resp_density_samp, cov_parms = cov_parms,
                      fixed_parms = beta_hat,
                      resids = resid_samp,
                      Sigma_ss = Sigma_ss,
                      AIC = AIC,
                      fitted_samp = fitted_samp,
                      summary_stlmm = summary_output,
                      minimax_vec = minimax_vec)
  
  class(stlmfit_obj) <- "stlmfit"
  return(stlmfit_obj)
}
