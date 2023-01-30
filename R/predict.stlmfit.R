#' Abundance Prediction
#' 
#' Predicts the abundance given an object from \code{stlmfit()}.
#' 
#' @param object is an object of class `stlmfit`
#' @param wts is the name of the column in the data frame with
#' the weights for which sites we want the prediction for. By default, the function will return the prediction for the total in the most recent time point given.
#' @param pred_level is the level used for a prediction interval
#' @param ... further arguments passed to or from other methods.
#' @return a list with the following \itemize{

#'   \item the prediction of the total abundance 
#'   \item the prediction variance
#'   \item a 90% prediction interval lower bound
#'   \item a 90% prediction interval upper bound
#'   \item a vector of site-by-site predictions for the unsampled sites
#'   \item a vector of site-by-site prediction variances for the unsampled sites
#'   \item the original data set appended with predictions
#'   \item the prediction interval level
#'  }
#'  
#' @examples 
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sp_de = 0.5, sp_range = 4, sp_ie = 0.5,
#'       t_de = 0.5, t_range = 0.7, t_ie = 0.5,
#'       spt_ie = 0.5)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data |> 
#' dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'  true = 1, false = 0))
#' samp_data$x <- rnorm(nrow(samp_data), 0, 1)
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords",
#' ycoord = "ycoords", tcoord = "times") 
#' predict(object = stlmfit_obj, wts = "predwts")
#' @import stats
#' @export

predict.stlmfit <- function(object, wts = NULL, pred_level = 0.90, ...) {
  
  if (inherits(object, "stlmfit") == FALSE) {
    stop("object must be of class stlmfit")
  }
  
  data <- object$data
  
  if (is.null(wts) == TRUE) {
    data$`_pred.wts` <- as.numeric((data$`_tindex` == max(data$`_tindex`,
                                               na.rm = TRUE)))
  } else if(is.character(wts)) {
    data$`_pred.wts` <- data[[wts]]
  } else if(length(wts) == nrow(data)) {
    data$`_pred.wts` <- wts
  } else {
    stop("if specified, wts must either be a string giving the name of a column in the data argument in stlmfit or must be a vector with the same number of rows as the data argument in stlmfit.")
  }
  
  X_all <- object$X_all
  ind_sa <- data$ind_sa
  ind_un <- data$ind_un
  Sigma_hat <- object$Sigma_hat
  beta_hat <- object$fixed_parms
  
  Sigma_ss <- Sigma_hat[ind_sa, ind_sa, drop = FALSE]
  Sigma_us <- Sigma_hat[ind_un, ind_sa, drop = FALSE]
  Sigma_su <- t(Sigma_us)
  Sigma_uu <- Sigma_hat[ind_un, ind_un, drop = FALSE]
  
  Sigma_ssi <- solve(Sigma_ss)
  
  nonzero_weight <- (!is.na(data$`_pred.wts`) & data$`_pred.wts` > 0)
  
  ## indicator for unsampled sites that are in the current year of 
  ## interest
  unsampcurrind <- ind_sa == 0 & nonzero_weight == TRUE
  
  
  X_ucurr <- X_all[unsampcurrind, , drop = FALSE]
  X_s <- X_all[ind_sa, , drop = FALSE]
  X_u <- X_all[ind_un, , drop = FALSE]
  ## matrix of covariance between all sites in the current year and
  ## all sites that were sampled
  Sigma_cs <- Sigma_hat[nonzero_weight == TRUE, ind_sa == 1]

    ## covariance matrix of sites in the current year
  Sigma_cc  <- Sigma_hat[nonzero_weight == TRUE,
                        nonzero_weight == TRUE]
  
  ## prediction indicator vector for all sampled sites
  bs_all <- data$`_pred.wts`[ind_sa == 1]
  
  ## prediction indicator vector for sampled sites in the current year
  b_s <- bs_all[bs_all > 0]
  
  ## prediction indicator vector for all unsampled sites 
  bu_all <- data$`_pred.wts`[ind_sa == 0]
  
  ## prediction indicator vector for
  ## the unsampled sites in the current year
  b_u <- bu_all[bu_all > 0]
  
  ## part 1 of the predictor
  p1 <- bs_all 
  
  ## covariance of the unsampled sites in the current year with all of the
  ## sampled sites
  Sigma_ucurrs <- Sigma_hat[unsampcurrind == 1, ind_sa == 1]
  Sigma_ssi <- solve(Sigma_hat[ind_sa == 1, ind_sa == 1])
  
  p2 <- t(b_u) %*% Sigma_ucurrs %*% Sigma_ssi
  
  p3 <- -t(b_u) %*% (Sigma_ucurrs %*% Sigma_ssi %*% X_s %*%
                      solve(t(X_s) %*% Sigma_ssi %*% X_s) %*% t(X_s) %*% Sigma_ssi)
  p4 <- t(b_u) %*% (X_ucurr %*%  solve(t(X_s) %*% Sigma_ssi %*% X_s) %*%
                     t(X_s) %*% Sigma_ssi)
  
  ## kriging weights
  tlambda <- (p1 + p2 + p3 + p4)
  
  resp_density <- object$resp_density
  
  totalpred <- tlambda %*% resp_density |> as.vector()
  
  ## prediction weights for sites in the current year (usually just a 
  ## vector of 1's if we want to predict the total for the current year).
  b_c <- c(b_s, b_u)
  
  predvar <- (tlambda %*% Sigma_ss %*% t(tlambda) -
    2 * t(b_c) %*% Sigma_cs %*% t(tlambda) +
    t(b_c) %*% Sigma_cc %*% b_c) |>
    as.vector()
  
  lb <- totalpred + 1 * stats::qnorm((1 - pred_level) / 2) * sqrt(predvar)
  ub <- totalpred + -1 * stats::qnorm((1 - pred_level) / 2) * sqrt(predvar)
  
  
  ## an equivalent calculation for the total:
  muhat_s <- X_s %*% beta_hat; muhat_u <- X_ucurr %*% beta_hat
  ## the predicted values for the sites that were not sampled
  zhatu <- Sigma_ucurrs %*% Sigma_ssi %*% (resp_density -
                                             muhat_s) + muhat_u
  ##totalpred_equiv <- sum(zhatu) + sum(density[ind_sa == 1 & data[[wts]] == 1])
  W <- t(X_u) - t(X_s) %*% Sigma_ssi %*% Sigma_su
  Vmat <- solve(t(X_s) %*% Sigma_ssi %*% X_s)
  sitecov <- Sigma_uu - Sigma_us %*% Sigma_ssi %*% Sigma_su +
    t(W) %*% Vmat %*% W
  sitevar <- diag(sitecov)
  
  
  data$predictions_ <- rep(NA, nrow(data))
  
  ## site by site predictions for past years are left out
  data$predictions_[unsampcurrind == 1] <- zhatu
  data$predictions_[ind_sa == 1] <- resp_density
  
  pred_obj <- list(totalpred = totalpred, predvar = predvar,
                   lb = lb, ub = ub, zhatu = zhatu, sitevar = sitevar,
                   data = data, pred_level = pred_level,
                   formula = object$summary_stlmm$Call)
  class(pred_obj) <- "predict.stlmfit"
  
  return(pred_obj)
  
}
