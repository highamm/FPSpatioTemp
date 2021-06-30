library(tidyverse)
moose_df <- read_csv("inst/moose_04_06.csv")
moose_samp <- moose_df %>% filter(!is.na(totalmoosena))

counts <- moose_df$totalmoosena
xcoords <- moose_df$centrlat ## need to change to TM
ycoords <- moose_df$centrlon

library(sptotal)
tm_obj <- LLtoTM(mean(moose_df$centrlon), lat = moose_df$centrlat,
                 lon = moose_df$centrlon)
moose_df$xTM <- tm_obj$xy[ ,1]
moose_df$yTM <- tm_obj$xy[ ,2]

cormodeltype <- "Exponential"
uniquecoords <- unique(cbind(moose_df$xTM, moose_df$yTM))
xcoord <- uniquecoords[ ,1]
ycoord <- uniquecoords[ ,2]
allcoord <- cbind(xcoord, ycoord)
distmat <- as.matrix(dist(allcoord, diag = TRUE, upper = TRUE))

N <- nrow(moose_df); n <- nrow(moose_samp)


moose_df


fulldf <- spatiotempdata

## vector of only the sampled counts
countssamp <- moose_samp$totalmoosena

uniquecoords <- unique(cbind(fulldf$xTM, fulldf$yTM))
distancemat <- as.matrix(stats::dist(uniquecoords))




##  make sure the function is working
##  can replace this step with a grid search to give `parmest` 
##  a better starting point than (1, 1, 1, 1).
m2LL.spatiotemp.ML(theta = c(1, 1, 1, 1), zcol = fulldf$totalmoosena,
                   XDesign = matrix(1, nrow = nrow(fulldf)), xcoord = uniquecoords[ ,1],
                   ycoord = uniquecoords[ ,2],
                   timepoints = unique(fulldf$Surveyyear), CorModel = "Exponential")

## find ML estimates or REML estimates
parmest <- stats::optim(c(1, 1, 1, 1), m2LL.spatiotemp.ML,
                        zcol = fulldf$totalmoosena,
                        XDesign = matrix(1, nrow = nrow(fulldf)),
                        xcoord = uniquecoords[ ,1],
                        ycoord = uniquecoords[ ,2],
                        timepoints = unique(fulldf$Surveyyear),
                        CorModel = "Exponential")

## fitted spatial and temporal parameters
nugget_hat <- exp(parmest$par[1])
parsil_hat <- exp(parmest$par[2])
range_hat <- exp(parmest$par[3])
rho_hat <- exp(parmest$par[4]) / (1 + exp(parmest$par[4]))

## estimated spatial covariance matrix
Sigmaspatest <- diag(nugget_hat, nrow = nrow(uniquecoords)) +
  parsil_hat * exp(-distancemat / range_hat)

## estimated AR(1) temporal covariance matrix
times <- min(fulldf$Surveyyear):max(fulldf$Surveyyear)
ntime <- length(times)
H <- abs(outer(times, times, "-")) 
Sigmatimeest <- rho_hat ^ H 

## including an additional term in the temporal model (sigma_time)
## would be redundant since the model is separable?

## full covariance matrix of separable model
Sigmaest <- kronecker(Sigmaspatest, Sigmatimeest)

fulldf$sampind

fulldf <- fulldf %>% mutate(sampind = if_else(!is.na(fulldf$totalmoosena), true = 1, false = 0))
Sigma.ss <- Sigmaest[fulldf$sampind == 1, fulldf$sampind == 1]
Sigma.us <- Sigmaest[fulldf$sampind == 0, fulldf$sampind == 1]
Sigma.su <- t(Sigma.us)
Sigma.uu <- Sigmaest[fulldf$sampind == 0, -fulldf$sampind == 0]

Sigma.ssi <- solve(Sigma.ss)

Xs <- matrix(1, sum(fulldf$sampind))

## define an indicator that denotes sites that are __unsampled__ and
## in the __current__ year (the sites that we want predictions for).

fulldf <- fulldf %>%
  mutate(predind = if_else(Surveyyear == 2006, true = 1, false = 0))
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
totalpred_equiv <- sum(zhatu) + sum(fulldf$totalmoosena[fulldf$sampind == 1 & fulldf$predind == 1])

fulldf$predictions <- rep(NA, nrow(fulldf))

fulldf$predictions[fulldf$unsampcurrind == 1] <- zhatu
fulldf$predictions[fulldf$sampind == 1] <- countssamp

as.vector(totalpred) + c(-1, 1) * 1.645 * sqrt(as.vector(predvar))

## 1729
## (1653, 1804)

## comparison to sptotal

moose2006 <- fulldf %>% filter(Surveyyear == 2006)
mod <- slmfit(totalmoosena ~ 1, data = moose2006, xcoordcol = "xTM",
       ycoordcol = "yTM")
predict(mod)
## 1793
## (1553, 2033)