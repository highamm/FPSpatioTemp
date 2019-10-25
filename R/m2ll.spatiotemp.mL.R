#' Covariance Parameter and Mean Estimation Function through Maximum Likelihood.
#'
#' The primary purpose of \code{m2ll.spatiotemp.mL.R} is to estimate the spatial
#' covariance parameters and fixed effects through maximum likelihood.
#'
#' @param theta is the parameter vector of (nugget, partialsill, range,
#' (vector of mean parameters), rho), where rho is the correlation for 
#' an AR(1) time series
#' @param zcol is the response vector of counts, with \code{NA} values 
#' for the missing values for sites that were not observed at a particular
#' time point.
#' @param XDesign is the design matrix containing the covariates used to predict animal or plant abundance (including a column of 1's for the intercept) for all
#' sites at all time points..
#' @param xcoord is a vector of the x spatial coordinates (in UTM)
#' @param ycoord is a vector of the y spatial coordinates (in UTM)
#' @param timepoints is a vector of timepoints for the model
#' @param CorModel is the geostatistical spatial correlation model to be used. See the \code{corModels} documentation for possible models to use.

#' @return A numeric output of minus 2 times the restricted log likelihood to be minimized by `optim` to obtain spatial parameter estimates.

m2LL.spatiotemp.ML <- function(theta, zcol, XDesign, xcoord, ycoord,
  timepoints, CorModel)
{
  nspat <- length(zcol) / length(timepoints)
  p <- length(XDesign[1,])
  nugget <- as.numeric(exp(theta[1]))
  parsil <- as.numeric(exp(theta[2]))
  range <- as.numeric(exp(theta[3]))
  beta <- matrix(as.numeric(theta[4:(length(theta) - 1)]))

  DM <- matrix(0, nspat, nspat)
  DM[lower.tri(DM)] <- stats::dist(as.matrix(cbind(xcoord, ycoord)))
  Dismat <- DM + t(DM)

  if (CorModel == "Exponential") {
    Sigmat <- parsil * exp(-Dismat / range)
    Cmat.spatial <- diag(nugget, nrow = nrow(Sigmat)) + Sigmat
  } ##else if (CorModel == "Gaussian") {
    ##Sigmat <- parsil * (corModelGaussian(Dismat, range))
    ##Cmat.nodet <- diag(nugget, nrow = nrow(Sigmat)) + Sigmat
  ##} else if (CorModel == "Spherical") {
  ##  Sigmat <- parsil * corModelSpherical(Dismat, range)
  ##  Cmat.nodet <- diag(nugget, nrow = nrow(Sigmat)) +
  ##    Sigmat
  ##}
  
  times <- min(timepoints):max(timepoints)
  ntime <- length(times)
  H <- abs(outer(times, times, "-")) 
  
  Sigmatime <- matrix(0, nrow = ntime, ncol = ntime)

  rhotime <- exp(theta[5]) / (1 + exp(theta[5]))
  
  Sigmatime <- 1 * rhotime ^ H 
  
  sampindx <- is.na(zcol) == FALSE
  
  Cmat.nodet <- kronecker(Cmat.spatial, Sigmatime)[sampindx, sampindx]
  
  Ci <- mginv(Cmat.nodet)

  zcolsamp <- zcol[sampindx]
  
  minus2loglik <- log(det(Cmat.nodet)) +
    (t(as.matrix(zcolsamp) - as.matrix(XDesign[sampindx, ] %*% beta))) %*%
    Ci %*%
    (as.matrix(zcolsamp) - as.matrix(XDesign[sampindx, ] %*% beta)) ##+
   ## n * log(2 * pi)

  return(as.numeric(minus2loglik))
}
