#' Abundance Prediction
#' 
#' Predicts the abundance given an object from \code{stlmfit()}.
#' 
#' @param object is an object of class `stlmfit`
#' @param wtscol is the name of the column in the data frame with
#' the weights for which sites we want the prediction for.
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
#'       sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'       sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'       sigma_nugget_spacetime = 0.5)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data %>% 
#' dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'  true = 1, false = 0))
#' samp_data$x <- rnorm(nrow(samp_data), 0, 1)
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoordcol = "xcoords",
#' ycoordcol = "ycoords", tcol = "times") 
#' predict(object = stlmfit_obj, wtscol = "predwts")
#' @import stats
#' @export

predict.stlmfit <- function(object, wtscol, pred_level = 0.90, ...) {
  
  if (inherits(object, "stlmfit") == FALSE) {
    stop("object must be of class stlmfit")
  }
  
  data <- object$data
  Xall <- object$Xall
  ind.sa <- data$ind.sa
  ind.un <- data$ind.un
  Sigmaest <- object$Sigmaest
  betahat <- object$parms$betahat
  
  Sigma.ss <- Sigmaest[ind.sa, ind.sa, drop = FALSE]
  Sigma.us <- Sigmaest[ind.un, ind.sa, drop = FALSE]
  Sigma.su <- t(Sigma.us)
  Sigma.uu <- Sigmaest[ind.un, ind.un, drop = FALSE]
  
  Sigma.ssi <- solve(Sigma.ss)
  
  ## indicator for unsampled sites that are in the current year of 
  ## interest
  unsampcurrind <- ind.sa == 0 & data[[wtscol]] == 1
  
  Xucurr <- Xall[unsampcurrind, , drop = FALSE]
  Xs <- Xall[ind.sa, , drop = FALSE]
  Xu <- Xall[ind.un, , drop = FALSE]
  ## matrix of covariance between all sites in the current year and
  ## all sites that were sampled
  Sigma.cs <- Sigmaest[data[[wtscol]] == 1, ind.sa == 1]
  
  ## covariance matrix of sites in the current year
  Sigma.cc <- Sigmaest[data[[wtscol]] == 1, data[[wtscol]] == 1]
  
  ## prediction indicator vector for all sampled sites
  bsall <-  data[[wtscol]][ind.sa == 1]
  
  ## prediction indicator vector for sampled sites in the current year
  bs <- bsall[bsall == 1]
  
  ## prediction indicator vector for all unsampled sites 
  buall <- data[[wtscol]][ind.sa == 0]
  
  ## prediction indicator vector for
  ## the unsampled sites in the current year
  bu <- buall[buall == 1]
  
  ## part 1 of the predictor
  p1 <- bsall 
  
  ## covariance of the unsampled sites in the current year with all of the
  ## sampled sites
  Sigma.ucurrs <- Sigmaest[unsampcurrind == 1, ind.sa == 1]
  Sigma.ssi <- solve(Sigmaest[ind.sa == 1, ind.sa == 1])
  
  p2 <- t(bu) %*% Sigma.ucurrs %*% Sigma.ssi
  
  p3 <- -t(bu) %*% (Sigma.ucurrs %*% Sigma.ssi %*% Xs %*%
                      solve(t(Xs) %*% Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi)
  p4 <- t(bu) %*% (Xucurr %*%  solve(t(Xs) %*% Sigma.ssi %*% Xs) %*%
                     t(Xs) %*% Sigma.ssi)
  
  ## kriging weights
  tlambda <- (p1 + p2 + p3 + p4)
  
  z.density <- object$z.density
  
  totalpred <- tlambda %*% z.density
  ## need to add in back area for this, if sites have unequal areas
  
  ## prediction weights for sites in the current year (usually just a 
  ## vector of 1's if we want to predict the total for the current year).
  bc <- c(bs, bu)
  
  predvar <- tlambda %*% Sigma.ss %*% t(tlambda) -
    2 * t(bc) %*% Sigma.cs %*% t(tlambda) +
    t(bc) %*% Sigma.cc %*% bc
  
  lb <- as.vector(totalpred) + 1 * stats::qnorm((1 - pred_level) / 2) * sqrt(as.vector(predvar))
  ub <- as.vector(totalpred) + -1 * stats::qnorm((1 - pred_level) / 2) * sqrt(as.vector(predvar))
  
  
  ## an equivalent calculation for the total:
  muhats <- Xs %*% betahat; muhatu <- Xucurr %*% betahat
  ## the predicted values for the sites that were not sampled
  zhatu <- Sigma.ucurrs %*% Sigma.ssi %*% (z.density -
                                             muhats) + muhatu
  ##totalpred_equiv <- sum(zhatu) + sum(density[ind.sa == 1 & data[[wtscol]] == 1])
  
  W <- t(Xu) - t(Xs) %*% Sigma.ssi %*% Sigma.su
  Vmat <- solve(t(Xs) %*% Sigma.ssi %*% Xs)
  sitecov <- Sigma.uu - Sigma.us %*% Sigma.ssi %*% Sigma.su +
    t(W) %*% Vmat %*% W
  sitevar <- diag(sitecov)
  
  
  data$predictions_ <- rep(NA, nrow(data))
  
  ## site by site predictions for past years are left out
  data$predictions_[unsampcurrind == 1] <- zhatu
  data$predictions_[ind.sa == 1] <- z.density
  
  pred_obj <- list(totalpred = totalpred, predvar = predvar,
                   lb = lb, ub = ub, zhatu = zhatu, sitevar = sitevar,
                   data = data, pred_level = pred_level)
  class(pred_obj) <- "predict.stlmfit"
  
  return(pred_obj)
  
}
