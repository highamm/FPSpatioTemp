#' Spatio-Temporal Model Fit
#' 
#' Estimates regression coefficients, spatial autocorrelation
#' parameters, and temporal autocorrelation parameters, given spatial
#' coordinates, time points and a model formula.
#' NOTE: Current version has not been tested with any fixed effects
#' in the design matrix X
#' NOTE: Current version only does Maximum Likelihood, though REML 
#' will be added.
#' 
#' @param formula is an \code{R} linear model formula specifying the
#' response variable as well as covariates for predicting the response on the unsampled sites.
#' @param data is a data frame or tibble with the response column, the covariates to be used for the block kriging, the spatio-temporal coordinates for all of the sites.
#' @param xcoordcol is the name of the column in the data frame with x coordinates or longitudinal coordinates.
#' @param ycoordcol is the name of the column in the data frame with y coordinates or latitudinal coordinates.
#' @param tcol is the name of the column in the data frame with 
#' the time points.
#' @param areacol is the name of the column in the data frame 
#' with the site areas
#' NOTE: package has not yet been tested for use of this argument.
#' @param CorModel a correlation model for the spatial correlation ("Exponential" is the only choice for now).
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
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'       sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'       sigma_nugget_spacetime = 0.5)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data %>%
#'  dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'   true = 1, false = 0))
#' samp_data$x <- rnorm(nrow(samp_data), 0, 1)
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoordcol = "xcoords",
#' ycoordcol = "ycoords", tcol = "times") 
#' @import stats
#' @export stlmfit
                                               
stlmfit <- function(formula, data, xcoordcol, ycoordcol, tcol,
                    areacol = NULL,
                    CorModel = "Exponential") {
  
  ## order the data so that sites are in the same order within each time
  
  data <- data %>% dplyr::ungroup() %>%
    dplyr::arrange_(tcol, xcoordcol, ycoordcol) 
  ## data %>% arrange(!! rlang::sym(c("tcol")))
  
  ## make data frame of only predictors
  datapredsonly <- data.frame(data[ ,all.vars(formula)[-1]])
  colnames(datapredsonly) <- all.vars(formula)[-1]
  
  ## design matrix for all sites
  Xall <- model.matrix(formula, model.frame(formula, data,
                                            na.action = stats::na.pass))
  
  ## get rid of observations with any missing predictors
  missingind <- base::apply(is.na(Xall), MARGIN = 1, FUN = sum)
  nmissing <- sum(missingind >= 1)
  
  datanomiss <- data[missingind == 0, ]
  
  xcoordsTM <- datanomiss[[xcoordcol]]
  ycoordsTM <- datanomiss[[ycoordcol]]
  
  if (is.null(areacol) == TRUE) {
    areavar <- rep(1, nrow(datanomiss))
  } else {
    areavar <- datanomiss[ ,areacol]
  }
  
  fullmf <- stats::model.frame(formula, na.action =
                                 stats::na.pass, data = datanomiss)
  
  yvar <- stats::model.response(fullmf, "numeric")
  density <- yvar / areavar
  
  ## remove any rows with missing values in any of the predictors
  formula.onlypreds <- formula[-2]
  
  X <- model.matrix(formula.onlypreds,
                    model.frame(formula.onlypreds, datanomiss,
                                na.action = stats::na.omit))
  
  ## divide data set into sampled sites and unsampled sites based
  ## on whether the response variable has a number (for sampled sites)
  ## or NA (for unsampled sites)
  
  ind.sa <- !is.na(yvar)
  data$ind.sa <- ind.sa
  ind.un <- is.na(yvar)
  data$ind.un <- ind.un
  
  data.sa <- datanomiss[ind.sa, ]
  data.un <- datanomiss[ind.un, ]
  
  m.un <- stats::model.frame(formula, data.un, na.action =
                               stats::na.pass)
  
  Xu <- model.matrix(formula.onlypreds,
                     model.frame(formula.onlypreds, data.un,
                                 na.action = stats::na.omit))
  
  ## sampled response values and design matrix
  m.sa <- stats::model.frame(formula, data.sa, na.action =
                               stats::na.omit)
  z.sa <- stats::model.response(m.sa)
  Xs <- stats::model.matrix(formula, m.sa)
  
  if (abs(det(t(Xs) %*% Xs)) < 1e-10) {
    stop("There are collinearity issues in the predictors. Remove collinear predictors and re-fit the model.")
  }
  
  
  z.density <- z.sa / areavar[ind.sa]
  n <- nrow(Xs)
  
  
  prednames <- colnames(Xs)
  
  ## x and y coordinates for sampled and unsampled sites
  x.sa <- xcoordsTM[ind.sa]
  y.sa <- ycoordsTM[ind.sa]
  x.un <- xcoordsTM[ind.un]
  y.un <- ycoordsTM[ind.un]
  
  ## number of sites that were sampled
  n.sa <- nrow(Xs)
  ## number of sites that were not sampled
  n.un <- nrow(Xu)
  
  ###########################
  ## stuff pertaining to time
  ###########################
  
  ## spatial sites must have __exactly__ the same coordinates
  ## issue with this approach: it's very sensitive to ordering.
  ## If the sites aren't in the same order within year, then there are 
  ## problems using the kronecker and [sampind, sampind] approach
  ## in m2ll.spatiotemp.ML
  
  uniquecoords <- unique(cbind(data[[xcoordcol]], data[[ycoordcol]]))
  distancemat <- as.matrix(stats::dist(uniquecoords))
  
  uniquetimes <- unique(data[[tcol]])
  times <- min(data[[tcol]]):max(data[[tcol]])
  ntime <- length(times)
  
  ## create Zs and Zt matrices
  N <- nrow(data)
  nspat <- N / length(uniquetimes)
  onetime <- diag(1, nspat) %>% as.data.frame()
  Zs <- onetime %>% slice(rep(row_number(), ntime)) %>% as.matrix()
  

  Zt <- lapply(1:ntime, matrix, data = 0, nrow = nspat, ncol = ntime)
  
  for (i in 1:ntime) {
    Zt[[i]][ ,i] <- 1
  }
  Zt <- do.call(rbind, Zt)
  
  H <- abs(outer(times, times, "-")) 
  
  DM <- matrix(0, nspat, nspat)
  DM[lower.tri(DM)] <- stats::dist(as.matrix(cbind(uniquecoords[ ,1], uniquecoords[ ,2])))
  Dismat <- DM + t(DM)
  
 
  
  varstart <- log(var(density, na.rm = TRUE) / 10)
  rangestart <- median(Dismat) / 2
  rhostart <- 0
  
  m2LL.spatiotemp.ML(theta = c(varstart, rangestart, varstart,
                              varstart, rhostart, varstart, varstart),
                     zcol = density,
                     XDesign = as.matrix(X),
                     xcoord = uniquecoords[ ,1],
                     ycoord = uniquecoords[ ,2],
                     timepoints = uniquetimes,
                     CorModel = "Exponential",
                     Zs = Zs, Zt = Zt, H = H, Dismat = Dismat)

  
  parmest <- stats::optim(c(varstart, rangestart, varstart,
                            varstart, rhostart, varstart, varstart),
                          m2LL.spatiotemp.ML,
                          zcol = density,
                          XDesign = as.matrix(X),
                          xcoord = uniquecoords[ ,1],
                          ycoord = uniquecoords[ ,2],
                          timepoints = uniquetimes,
                          CorModel = "Exponential",
                          Zs = Zs, Zt = Zt, H = H, Dismat = Dismat)
  
  sigma_parsil_spat_hat <- exp(parmest$par[1])
  range_hat <- exp(parmest$par[2])
  sigma_nugget_spat_hat <- exp(parmest$par[3])
  
  sigma_parsil_time_hat <- exp(parmest$par[4])
  rhotime_hat <- exp(parmest$par[5]) / (1 + exp(parmest$par[5]))
  sigma_nugget_time_hat  <- exp(parmest$par[6])
  
  sigma_nugget_spacetime_hat <- exp(parmest$par[7])
  
  if (CorModel == "Exponential") {
    Rs_hat <- exp(-Dismat / range_hat)
  }
  
  comp_1 <- sigma_parsil_spat_hat * Zs %*% Rs_hat %*% t(Zs)
  comp_2 <- sigma_nugget_spat_hat * Zs %*% t(Zs)
  
  Rt_hat <- rhotime_hat ^ H 
  comp_3 <- sigma_parsil_time_hat * Zt %*% Rt_hat %*% t(Zt)
  comp_4 <- sigma_nugget_time_hat * Zt %*% t(Zt)
  
  comp_5 <- diag(sigma_nugget_spacetime_hat, nrow = N)
  
  Sigmaest <- comp_1 + comp_2 + comp_3 + comp_4 + comp_5
  Sigma.ss <- Sigmaest[ind.sa, ind.sa, drop = FALSE]
  Sigma.ssi <- solve(Sigma.ss) ## can output this so it doesn't
  ## also get inverted in predict.stlmfit()
  
  betahat <- solve(t(Xs) %*% Sigma.ssi %*% Xs) %*%
    t(Xs) %*% Sigma.ssi %*%
    z.density
  
  stlmfit_obj <- list(data = data, Sigmaest = Sigmaest, Xall = Xall, 
                      z.density = z.density, parms = list(betahat = betahat, sigma_parsil_spat = sigma_parsil_spat_hat,
                               range = range_hat,
                               sigma_nugget_spat = sigma_nugget_spat_hat,
                               sigma_parsil_time = sigma_parsil_time_hat,
                               rho = rhotime_hat,
                               sigma_nugget_time = sigma_nugget_time_hat,
                               sigma_nugget_spacetime = sigma_nugget_spacetime_hat))
  
  class(stlmfit_obj) <- "stlmfit"
  return(stlmfit_obj)
}
