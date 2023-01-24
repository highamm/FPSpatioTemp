#' Spatiotemporal Model Fit
#' 
#' Estimates regression coefficients, spatial autocorrelation
#' parameters, and temporal autocorrelation parameters, given spatial
#' coordinates, time points and a model formula.
#' NOTE: Current version has not been tested with any fixed effects
#' in the design matrix X
#' 
#' @param formula is an \code{R} linear model formula specifying the
#' response variable as well as covariates for predicting the response on the unsampled sites.
#' @param data is a data frame or tibble with the response column, the covariates to be used for the block kriging, the spatio-temporal coordinates for all of the sites.
#' @param xcoord is the name of the column in the data frame with x coordinates or longitudinal coordinates.
#' @param ycoord is the name of the column in the data frame with y coordinates or latitudinal coordinates.
#' @param tcoord is the name of the column in the data frame with 
#' the time points.
#' @param areacol is the name of the column in the data frame 
#' with the site areas
#' NOTE: package has not yet been tested for use of this argument.
#' @param cor_model_sp a correlation model for the spatial correlation. Options are \code{"exponential"}, \code{"gaussian"}, \code{"triangular"}, or \code{"cosine"}.
#' @param cor_model_t a correlation model for the temporal correlation Options are \code{"exponential"}, \code{"gaussian"}, \code{"triangular"}, or \code{"cosine"}.
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
#' NOTE: This function will eventually be split into stlmfit() to fit
#' the model and estcov() as a helper in model fitting
#' @examples 
#' set.seed(07262022)
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sp_de = 0.5, sp_range = 4, sp_ie = 0.5,
#'       t_de = 0.5, t_range = 0.7, t_ie = 0.5,
#'       spt_ie = 0.5)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data |>
#'  dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'   true = 1, false = 0))
#' samp_data$x <- rnorm(nrow(samp_data), 0, 1)
#' samp_data <- samp_data
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#' @import stats
#' @export stlmfit
             
# stlmfit(formula = response_na ~ 1, data = samp_data,
# xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
                     
stlmfit <- function(formula, data, xcoord, ycoord, tcoord,
                    areacol = NULL,
                    cor_model_sp = "exponential",
                    cor_model_t = "exponential") {

  ## order the data so that sites are in the same order within each time
  ## need to mess with this some more for NSE
  xcoord <- substitute(xcoord)
  ycoord <- substitute(ycoord)
  tcoord <- substitute(tcoord)
  

  order_spt_obj <- order_spt(data, xcoord, ycoord, tcoord)
  
  data_ordered <- order_spt_obj$full_data
  
  ## data_obs contains rows in the original sampling frame of data.
  data_og <- data_ordered |> dplyr::filter(.data$.observed == TRUE)
  
  
  ## make data frame of only predictors
  data_preds_only <- data.frame(data_og[ ,all.vars(formula)[-1]])
  colnames(data_preds_only) <- all.vars(formula)[-1]
  
  
  ## design matrix for all sites
  X_all <- model.matrix(formula, model.frame(formula, data_og,
                                             na.action = stats::na.pass))
  
  ## get rid of observations with any missing predictors
  missing_ind <- base::apply(is.na(X_all), MARGIN = 1, FUN = sum)
  n_missing <- sum(missing_ind >= 1)
  
  ## data_obs has all observations except for....
  ## missing covariates and has already removed the "extra" 
  ## rows from order_spt()
  
  data_obs <- data_og[missing_ind == 0, ]
  
  

  
  
  if (is.null(areacol) == TRUE) {
    area_var <- rep(1, nrow(data_obs))
  } else {
    area_var <- data_obs[ ,areacol]
  }
  
  full_mf <- stats::model.frame(formula, na.action =
                                 stats::na.pass, data = data_obs)
  
  resp_var <- stats::model.response(full_mf, "numeric")
  resp_density <- resp_var / area_var
  
  ind_sa <- !is.na(resp_var)
  data_obs$ind_sa <- ind_sa
  ind_un <- is.na(resp_var)
  data_obs$ind_un <- ind_un
  
  resp_density_samp <- resp_density[ind_sa]
  
  
  ## build_z_mats filters out the rows where
  ## .observed == FALSE at the end
  z_mats <- build_z_mats(data_ord = data_ordered)
  
  ## create Zs and Zt matrices
  # N <- nrow(data)
  # nspat <- N / length(uniquetimes)
  
  Z_sp <- z_mats$Z_sp
  Z_t <- z_mats$Z_t
  
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
  
  # sp_dist_mat <- dist(cbind(data[[xcoord]], data[[ycoord]]),
  #                     diag = TRUE, upper = TRUE) |>
  #   as.matrix()
  # 
  # t_dist_mat <- dist(data[[tcoord]],
  #                    diag = TRUE, upper = TRUE) |>
  #   as.matrix()
  
  h_sp_small <- order_spt_obj$h_sp_small
  h_t_small <- order_spt_obj$h_t_small
  
  sp_range_initial <- median(h_sp_small)
  t_range_initial <- median(h_t_small)
  
  sim_initial <-  DumelleEtAl2021STLMM::make_covparam_object(
    s_de = sp_de_initial, s_ie = sp_ie_initial,
    t_ie = t_ie_initial, t_de = t_de_initial,
    st_de = spt_de_initial, st_ie = spt_ie_initial,
    s_range = sp_range_initial, t_range = t_range_initial,
    stcov = "productsum", estmethod = "reml")
  
  data_sa <- data_obs[ind_sa, ]
  data_un <- data_obs[ind_un, ]
  
  ## mike uses the -3h / range parameterization in stlmm()
  ## want to convert this to -h / range_1 after estimation
  
  # data_test <- data |> mutate(resp2 = rnorm(12, 0, 1))

  ## CANNOT HAVE A COLUMN NAMED index, spindex, or tindex
  fast_est <- DumelleEtAl2021STLMM::stlmm(
    formula = formula,
    data = data_sa |> dplyr::select(-.data$index, -.data$spindex, -.data$tindex),
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

  
  parm_est <- c(unclass(fast_est$CovarianceParameters), use.names = FALSE)
  sp_de_hat <- parm_est[1]
  sp_range_hat <- parm_est[7] / 3  ## dividing by 3 to cancel dumelle alternative parameterization
  sp_ie_hat <- parm_est[2]
  
  t_de_hat <- parm_est[3]
  t_range_hat <- parm_est[8] / 3
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
  

  comp_1 <- sp_de_hat * Z_sp %*% R_sp_hat %*% t(Z_sp)
  comp_2 <- sp_ie_hat * Z_sp %*% t(Z_sp)
  
  R_t_hat <- build_r(cov_type = cor_model_t, range = t_range_hat,
                     dist_mat = h_t_small)

  comp_3 <- t_de_hat * Z_t %*% R_t_hat %*% t(Z_t)
  comp_4 <- t_ie_hat * Z_t %*% t(Z_t)
  
  comp_5 <- spt_de_hat * (Z_sp %*% R_sp_hat %*% t(Z_sp)) * (Z_t %*% R_t_hat %*% t(Z_t))
  
  comp_6 <- diag(spt_ie_hat, nrow = nrow(data_obs))

  
  Sigma_hat <- comp_1 + comp_2 + comp_3 + comp_4 + comp_5 + comp_6
  
  Sigma_ss <- Sigma_hat[ind_sa, ind_sa, drop = FALSE]
  Sigma_ssi <- solve(Sigma_ss) ## can output this so it doesn't
  ## also get inverted in predict.stlmfit()
  
  beta_hat <- as.vector(fast_est$Coefficients)
  names(beta_hat) <- c("Intercept", all.vars(formula)[-1])
  
  stlmfit_obj <- list(data = data_obs, Sigma_hat = Sigma_hat, X_all = X_all, 
                      resp_density = resp_density_samp, cov_parms = cov_parms,
                      fixed_parms = beta_hat)
  
  class(stlmfit_obj) <- "stlmfit"
  return(stlmfit_obj)
}
