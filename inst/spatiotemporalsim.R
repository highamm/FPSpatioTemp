## separable models: can separate the spatial and temporal covariance 
## components so that they are multiplied in the errors

## Suppose that we have a separable model and we have data on a 20 x 20
## grid over 5 years (for 2000 observations)

library(MASS)

beta0 <- 10

xcoords <- 1:20; ycoords <- 1:20
allcoords <- expand.grid(xcoords, ycoords)
distancemat <- as.matrix(dist(allcoords))
nspat <- nrow(distancemat)
range <- 10; parsil <- 4
Sigma <- parsil * exp(-distancemat / range)

ntime <- 5
times <- 1:ntime

Sigmatime <- matrix(0, nrow = ntime, ncol = ntime)
H <- abs(outer(times, times, "-")) 
##assume sigma is 1: we have a separable model
rhotime <- 0.7; sigmatime <- 1

Sigmatime <- sigmatime * rhotime ^ H 
Sigmaboth <- kronecker(Sigma, Sigmatime)
Sigmaboth[1:5, 1:5]

ntotal <- ntime * nspat
epsilon <- mvrnorm(n = 1, mu = rep(0, ntotal), Sigma = Sigmaboth)

counts <- round(beta0 + epsilon, 0)

## get rid of some of the observations
nsamp <- 400
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

source("m2ll.spatiotemp.ml.R")
source("mginv.R")

m2LL.spatiotemp.ML(theta = c(1, 1, 1, 2, 1), zcol = countswithna,
  XDesign = matrix(1, nrow = length(counts)), xcoord = allcoords[ ,1],
  ycoord = allcoords[ ,2],
  timepoints = 1:5, CorModel = "Exponential")

parmest <- optim(c(1, 1, 1, 2, 1), m2LL.spatiotemp.ML,
  zcol = countswithna,
  XDesign = matrix(1, nrow = length(counts)), xcoord = allcoords[ ,1],
  ycoord = allcoords[ ,2],
  timepoints = 1:5, CorModel = "Exponential")

nugget_hat <- exp(parmest$par[1])
parsil_hat <- exp(parmest$par[2])
range_hat <- exp(parmest$par[3])
beta_hat <- parmest$par[4]
rho_hat <- exp(parmest$par[5]) / (1 + exp(parmest$par[5]))


Sigmaspatest <- diag(nugget_hat, nrow = nrow(distancemat)) +
  parsil_hat * exp(-distancemat / range_hat)
Sigmatimeest <- rho_hat ^ H 
Sigmaest <- kronecker(Sigmaspatest, Sigmatimeest)


Sigma.ss <- Sigmaest[indxkeep, indxkeep]
Sigma.us <- Sigmaest[-indxkeep, indxkeep]
Sigma.su <- t(Sigma.us)
Sigma.uu <- Sigma[-indxkeep, -indxkeep]


Sigma.ssi <- solve(Sigma.ss)

## the generalized least squares regression coefficient estimates
betahat <- beta_hat

## estimator for the mean vector



ntotal - nspat
yearind <- rep(0, ntotal)
yearind[(ntotal - nspat + 1):ntotal] <- 1


bsall <- yearind[is.na(countswithna) == FALSE]
bucurr <- rep(1, length(yearind[is.na(countswithna) == TRUE & 
    yearind == 1]))

p1 <- bsall

sampind <- matrix(is.na(countswithna) == FALSE)
ucurrind <- matrix(yearind == 1 & sampind == 0)

Sigma.ucurrs <- Sigmaest[ucurrind, indxkeep]
Sigma.ssi <- solve(Sigmaest[indxkeep, indxkeep])

p2 <- t(bucurr) %*% Sigma.ucurrs %*% Sigma.ssi

Xs <- matrix(1, length(counts[indxkeep]))
Xucurr <- matrix(1, length(counts[ucurrind]))

p3 <- t(bucurr) %*% (Sigma.ucurrs %*% Sigma.ssi %*% Xs %*% solve(t(Xs) %*%
    Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi + 
    Xucurr %*% solve(t(Xs) %*%
        Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi)
(p1 + p2 + p3) %*% countssamp
sum(counts[1601:2000])
 
(Sigma.ucurrs %*% Sigma.ssi + (Sigma.ucurrs %*% Sigma.ssi %*% Xs %*% solve(t(Xs) %*%
    Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi + 
    Xucurr %*% solve(t(Xs) %*%
        Sigma.ssi %*% Xs) %*% t(Xs) %*% Sigma.ssi)) %*% matrix(countssamp)


## this is right but something is off with the above equation...
muhats <- Xs %*% betahat; muhatu <- Xucurr %*% betahat
## the predicted values for the sites that were not sampled
zhatu <- Sigma.ucurrs %*% Sigma.ssi %*% (countssamp -
    muhats) + muhatu

sum(zhatu) + sum(countswithna[sampind == 1 & yearind == 1])
sum(counts[1601:2000])

## next step: translate the above zhatu to lambda and find bsall and bu
## so that the prediction comes out accurately 