#' Finite Population Spatio-temporal Prediction
#'
#' \code{fp.predict} predicts the population total for the current year of sampling from sampled sites in the current year as well as sampled sites in past years.
#'
#' @param spatiotempdata is an object from the `give.nas` function with simulated counts, spatial coordinates, and time points.
#' @return A list with (1) the data frame with the predicted site-by-site predictions, (2) the true total for the current year, (3) the estimated total for the current year, and (4) the prediction variance.
#' @export fp.predict

fp.predict <- function(spatiotempdata) {
  
  ## separable models: can separate the spatial and temporal covariance 
  ## components so that they are multiplied in the errors
  ## To start, assume that we have a separable model
  
  ## grab the object from the give.nas function. This is a data frame
  ## with coordianates, true counts, observed counts, a sampling indicator,
  ## and a prediction indicator (a 1 indicates a site in the present that
  ## we want to include in our prediction for the current population total)
  
  fulldf <- spatiotempdata
  
  ## total number of sites sampled, past and present
  nsamp <- sum(fulldf$sampind)
  
  ## vector of only the sampled counts
  countssamp <- with(fulldf,
    counts[sampind == 1])
  
  uniquecoords <- unique(cbind(fulldf$xcoords, fulldf$ycoords))
  distancemat <- as.matrix(stats::dist(uniquecoords))
  
  ##  make sure the function is working
  ##  can replace this step with a grid search to give `parmest` 
  ##  a better starting point than (1, 1, 1, 1).
  m2LL.spatiotemp.ML(theta = c(1, 1, 1, 1), zcol = fulldf$obscounts,
    XDesign = matrix(1, nrow = nrow(fulldf)), xcoord = uniquecoords[ ,1],
    ycoord = uniquecoords[ ,2],
    timepoints = unique(fulldf$ts), CorModel = "Exponential")
  
  ## find ML estimates or REML estimates
  parmest <- stats::optim(c(1, 1, 1, 1), m2LL.spatiotemp.REML,
    zcol = fulldf$obscounts,
    XDesign = matrix(1, nrow = nrow(fulldf)), xcoord = uniquecoords[ ,1],
    ycoord = uniquecoords[ ,2],
    timepoints = unique(fulldf$ts), CorModel = "Exponential")
  
  ## fitted spatial and temporal parameters
  nugget_hat <- exp(parmest$par[1])
  parsil_hat <- exp(parmest$par[2])
  range_hat <- exp(parmest$par[3])
  rho_hat <- exp(parmest$par[4]) / (1 + exp(parmest$par[4]))
  
  ## estimated spatial covariance matrix
  Sigmaspatest <- diag(nugget_hat, nrow = nrow(uniquecoords)) +
    parsil_hat * exp(-distancemat / range_hat)
  
  ## estimated AR(1) temporal covariance matrix
  times <- min(fulldf$ts):max(fulldf$ts)
  ntime <- length(times)
  H <- abs(outer(times, times, "-")) 
  Sigmatimeest <- rho_hat ^ H 
  
  ## including an additional term in the temporal model (sigma_time)
  ## would be redundant since the model is separable?
  
  ## full covariance matrix of separable model
  Sigmaest <- kronecker(Sigmaspatest, Sigmatimeest)
  
  
  Sigma.ss <- Sigmaest[fulldf$sampind == 1, fulldf$sampind == 1]
  Sigma.us <- Sigmaest[fulldf$sampind == 0, fulldf$sampind == 1]
  Sigma.su <- t(Sigma.us)
  Sigma.uu <- Sigmaest[fulldf$sampind == 0, -fulldf$sampind == 0]
  
  Sigma.ssi <- solve(Sigma.ss)
  
  Xs <- matrix(1, sum(fulldf$sampind))
  
  ## define an indicator that denotes sites that are __unsampled__ and
  ## in the __current__ year (the sites that we want predictions for).
  
  fulldf$unsampcurrind <- fulldf$sampind == 0 & fulldf$predind == 1
  
  Xucurr <- matrix(1, sum(fulldf$unsampcurrind))
  
  ## the generalized least squares regression coefficient estimates
  betahat <- solve(t(Xs) %*% Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi %*%
    countssamp
  
  ## matrix of covariance between all sites in the current year and
  ## all sites that were sampled
  Sigma.cs <- Sigmaest[fulldf$predind == 1, fulldf$sampind == 1]
  
  ## covariance matrix of sites in the current year
  Sigma.cc <- Sigmaest[fulldf$predind == 1, fulldf$predind == 1]
  
  ## prediction indicator vector for all sampled sites
  bsall <-  fulldf$predind[fulldf$sampind == 1]
  
  ## prediction indicator vector for sampled sites in the current year
  bs <- bsall[bsall == 1]
  
  ## prediction indicator vector for all unsampled sites 
  buall <- fulldf$predind[fulldf$sampind == 0]
  
  ## prediction indicator vector for the unsampled sites in the current year
  bu <- buall[buall == 1]
  
  ## part 1 of the predictor
  p1 <- bsall 
  
  ## covariance of the unsampled sites in the current year with all of the
  ## sampled sites
  Sigma.ucurrs <- Sigmaest[fulldf$unsampcurrind == 1, fulldf$sampind == 1]
  Sigma.ssi <- solve(Sigmaest[fulldf$sampind == 1, fulldf$sampind == 1])
  
  p2 <- t(bu) %*% Sigma.ucurrs %*% Sigma.ssi
  
  p3 <- -t(bu) %*% (Sigma.ucurrs %*% Sigma.ssi %*% Xs %*%
      solve(t(Xs) %*% Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi)
  p4 <- t(bu) %*% (Xucurr %*%  solve(t(Xs) %*% Sigma.ssi %*% Xs) %*%
      t(Xs) %*% Sigma.ssi)
  
  ## kriging weights
  tlambda <- (p1 + p2 + p3 + p4)
  
  ## prediction for the total in the current year
  totalpred <- tlambda %*% countssamp
  
  ## prediction weights for sites in the current year (usually just a 
  ## vector of 1's if we want to predict the total for the current year).
  bc <- c(bs, bu)
  
  predvar <- tlambda %*% Sigma.ss %*% t(tlambda) -
    2 * t(bc) %*% Sigma.cs %*% t(tlambda) +
    t(bc) %*% Sigma.cc %*% bc
  
  
  ## an equivalent calculation for the total:
  muhats <- Xs %*% betahat; muhatu <- Xucurr %*% betahat
  ## the predicted values for the sites that were not sampled
  zhatu <- Sigma.ucurrs %*% Sigma.ssi %*% (countssamp -
      muhats) + muhatu
  totalpred_equiv <- sum(zhatu) + sum(fulldf$obscounts[fulldf$sampind == 1 & fulldf$predind == 1])
  
  truetotal <- sum(fulldf$counts[fulldf$predind == 1])
  
  fulldf$predictions <- rep(NA, nrow(fulldf))
  
  fulldf$predictions[fulldf$unsampcurrind == 1] <- zhatu
  fulldf$predictions[fulldf$sampind == 1] <- countssamp
  
  pred_obj <- list(fulldf, truetotal, totalpred, predvar)
  return(pred_obj)
}