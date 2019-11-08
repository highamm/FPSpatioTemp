## separable models: can separate the spatial and temporal covariance 
## components so that they are multiplied in the errors

## Suppose that we have a separable model and we have data on a 20 x 20
## grid over 5 years (for 2000 observations)









## get rid of some of the observations
nsamp <- 500
indxkeep <- sample(1:ntotal, nsamp, replace = FALSE)

countswithna <- counts
countswithna[-indxkeep] <- NA

countssamp <- counts[indxkeep]
Sigmasamp <- Sigmaboth[indxkeep, indxkeep]

countstemporary <- countswithna[1:((ntime - 1) * nspat)]
countsallnew_compold <- c(countstemporary[!is.na(countstemporary)],
  countswithna[((ntime - 1) * nspat + 1): ntotal])

Sigmaallnew_compold <- Sigmaboth[c(countswithna[1:((ntime - 1) * nspat)][!is.na(countswithna[1:((ntime - 1) * nspat)])], countswithna[((ntime - 1) * nspat + 1): ntotal]), 
  c(countswithna[1:((ntime - 1) * nspat)][!is.na(countswithna[1:((ntime - 1) * nspat)])], countswithna[((ntime - 1) * nspat + 1): ntotal])]
nrow(Sigmaallnew_compold)

source("/Desktop/FPSpatioTemp/R/m2ll.spatiotemp.ml.R")
source("mginv.R")

m2LL.spatiotemp.ML(theta = c(1, 1, 1, 1), zcol = countswithna,
  XDesign = matrix(1, nrow = length(counts)), xcoord = allcoords[ ,1],
  ycoord = allcoords[ ,2],
  timepoints = 1:5, CorModel = "Exponential")

parmest <- optim(c(1, 1, 1, 1), m2LL.spatiotemp.ML,
  zcol = countswithna,
  XDesign = matrix(1, nrow = length(counts)), xcoord = allcoords[ ,1],
  ycoord = allcoords[ ,2],
  timepoints = 1:5, CorModel = "Exponential")

nugget_hat <- exp(parmest$par[1])
parsil_hat <- exp(parmest$par[2])
range_hat <- exp(parmest$par[3])
rho_hat <- exp(parmest$par[4]) / (1 + exp(parmest$par[4]))


Sigmaspatest <- diag(nugget_hat, nrow = nrow(distancemat)) +
  parsil_hat * exp(-distancemat / range_hat)
Sigmatimeest <- rho_hat ^ H 
Sigmaest <- kronecker(Sigmaspatest, Sigmatimeest)


Sigma.ss <- Sigmaest[indxkeep, indxkeep]
Sigma.us <- Sigmaest[-indxkeep, indxkeep]
Sigma.su <- t(Sigma.us)
Sigma.uu <- Sigma[-indxkeep, -indxkeep]


Sigma.ssi <- solve(Sigma.ss)

Xs <- matrix(1, length(counts[indxkeep]))
Xucurr <- matrix(1, length(counts[ucurrind]))

## the generalized least squares regression coefficient estimates
betahat <- solve(t(Xs) %*% Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi %*% countssamp
## estimator for the mean vector



ntotal - nspat
yearind <- rep(0, ntotal)
yearind[(ntotal - nspat + 1):ntotal] <- 1

Sigma.cs <- Sigmaest[yearind == 1, indxkeep]
Sigma.cc <- Sigmaest[yearind == 1, yearind == 1]

bsall <- yearind[is.na(countswithna) == FALSE]
bs <- bsall[bsall == 1]
bucurr <- rep(1, length(yearind[is.na(countswithna) == TRUE & 
    yearind == 1]))

p1 <- bsall 

sampind <- matrix(is.na(countswithna) == FALSE)
ucurrind <- matrix(yearind == 1 & sampind == 0)

Sigma.ucurrs <- Sigmaest[ucurrind, indxkeep]
Sigma.ssi <- solve(Sigmaest[indxkeep, indxkeep])

p2 <- t(bucurr) %*% Sigma.ucurrs %*% Sigma.ssi



p3 <- -t(bucurr) %*% (Sigma.ucurrs %*% Sigma.ssi %*% Xs %*%  solve(t(Xs) %*% Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi)
p4 <- t(bucurr) %*% (Xucurr %*%  solve(t(Xs) %*% Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi)
tlambda <- (p1 + p2 + p3 + p4)
tlambda %*% countssamp

bc <- c(bs, bucurr)
predvar <- tlambda %*% Sigma.ss %*% t(tlambda) - 2 * t(bc) %*% Sigma.cs %*% bc +
  t(bc) %*% Sigma.cc %*% bc
sqrt(predvar)

sum(counts[1601:2000])
 


## this is right but something is off with the above equation...
muhats <- Xs %*% betahat; muhatu <- Xucurr %*% betahat
## the predicted values for the sites that were not sampled
zhatu <- Sigma.ucurrs %*% Sigma.ssi %*% (countssamp -
    muhats) + muhatu



bucurr %*% Sigma.ucurrs %*% Sigma.ssi %*% countssamp -
  bucurr %*% Sigma.ucurrs %*% Sigma.ssi %*% muhats +
  bucurr %*% muhatu +
  sum(countssamp[yearind[is.na(countswithna) == FALSE] == 1])
sum(zhatu) + sum(countssamp[indxkeep == TRUE])
sum(counts[1601:2000])

## next step: code this function so that it is within the body of the package
## and then document this function appropriately
## next step: clean up notation in write-up 
## next step: clean up the function so that it reads better and isn't just awful