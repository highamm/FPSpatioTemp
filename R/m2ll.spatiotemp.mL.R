#' Covariance Parameter and Mean Estimation Function through Maximum Likelihood.
#'
#' The primary purpose of \code{m2LL.spatiotemp.mL.R} is to estimate the spatial
#' covariance parameters and fixed effects through maximum likelihood.
#'
#' @param theta is the parameter vector of (spatial partial sill,
#' range, spatial nugget, range, temporal partial sill, temporal autocorrelation, temporal nugget, and spatiotemporal nugget)
#' @param zcol is the response vector of counts, with \code{NA} values 
#' for the missing values for sites that were not observed at a particular
#' time point.
#' @param XDesign is the design matrix containing the covariates used to predict animal or plant abundance (including a column of 1's for the intercept) for all
#' sites at all time points..
#' @param xcoord is a vector of the x spatial coordinates (in UTM)
#' @param ycoord is a vector of the y spatial coordinates (in UTM)
#' @param timepoints is a vector of timepoints for the model
#' @param CorModel is the geostatistical spatial correlation model to be used. See the \code{corModels} documentation for possible models to use.
#' @param Zs is the design matrix of spatial random effects.
#' @param Zt is the design matrix of temporal random effects.
#' @param H is a matrix used to generate Rt.
#' @param Dismat is the spatial distance matrix for unique spatial coordinates.

#' @return A numeric output of minus 2 times the restricted log likelihood to be minimized by `optim` to obtain spatial parameter estimates.
#' @export m2LL.spatiotemp.ML

m2LL.spatiotemp.ML <- function(theta, zcol, XDesign, xcoord, ycoord,
                               timepoints, CorModel, Zs, Zt, H, Dismat)
{
  
  p <- length(XDesign[1, ])
  sigma_parsil_spat <- as.numeric(exp(theta[1]))
  range <- as.numeric(exp(theta[2]))
  sigma_nugget_spat  <- as.numeric(exp(theta[3]))
  ##beta <- matrix(as.numeric(theta[4:(length(theta) - 1)]))
  
  n <- nrow(XDesign)
  
  if (CorModel == "Exponential") {
    Rs <- exp(-Dismat / range)
  }
  
  comp_1 <- sigma_parsil_spat * Zs %*% Rs %*% t(Zs)
  comp_2 <- sigma_nugget_spat * Zs %*% t(Zs)
  
  sigma_parsil_time <- as.numeric(exp(theta[4]))
  rhotime <- as.numeric(exp(theta[5]) / (1 + exp(theta[5])))
  sigma_nugget_time  <- as.numeric(exp(theta[6]))
  
  Rt <- rhotime ^ H 
  comp_3 <- sigma_parsil_time * Zt %*% Rt %*% t(Zt)
  comp_4 <- sigma_nugget_time * Zt %*% t(Zt)
  
  sigma_nugget_spacetime <- as.numeric(exp(theta[7]))
  comp_5 <- diag(sigma_nugget_spacetime, nrow = n)
  
  Sigma <- comp_1 + comp_2 + comp_3 + comp_4 + comp_5
  
  sampindx <- is.na(zcol) == FALSE
  Sigma_samponly <- Sigma[sampindx, sampindx]
  
  Sigma_samponly_i <- mginv(Sigma_samponly)
  
  zcolsamp <- zcol[sampindx]
  
  beta <- matrix(mginv(t(XDesign[sampindx, ]) %*% Sigma_samponly_i %*%
                         XDesign[sampindx, ]) %*%
                   t(XDesign[sampindx, ]) %*% Sigma_samponly_i %*%
                   as.matrix(zcolsamp))
  minus2loglik <- log(det(Sigma_samponly)) +
    (t(as.matrix(zcolsamp) - as.matrix(XDesign[sampindx, ] %*% beta))) %*%
    Sigma_samponly_i %*%
    (as.matrix(zcolsamp) - as.matrix(XDesign[sampindx, ] %*% beta)) ##+
  ## n * log(2 * pi)
  
  return(as.numeric(minus2loglik))
}
