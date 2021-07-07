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
#' @param wtscol is the name of the column in the data frame with
#' the weights for which sites we want the prediction for.
#' @param areacol is the name of the column in the data frame 
#' with the site areas
#' NOTE: package has not yet been tested for use of this argument.
#' a list of class \code{slmfit} with
#' @return a list with 
#'   \item the prediction of the total abundance
#'   \item the prediction variance
#'   \item a 90% prediction interval lower bound
#'   \item a 90% prediction interval upper bound
#'   \item a vector of site-by-site predictions for the unsampled sites
#'   \item a vector of site-by-site prediction variances for the unsampled sites
#'   \item a list with the 7 spatio-temporal covariance parameter
#'    estimates for Exponential spatial covariance and AR(1) temporal
#'    covariance and a sum-with-error linear mixed model (spatial
#'    partial sill, range, spatial nugget, temporal partial sill,
#'    autocorrelation parameter, temporal nugget, and spatiotemporal
#'    nugget).
#' NOTE: This function will eventually be split into stlmfit() to fit
#' the model, estcov() as a helper in model fitting, 
#' and predict.stlmfit(), in much the same way that sptotal is split.
#' @import stats
#' @import dplyr
#' @export stlmfit


# moose_df2 <- moose_df %>% sample_n(nrow(moose_df))
# formula <- totalmoosena ~ 1
# data <- moose_df2
# xcoordcol <- "xTM"
# ycoordcol <- "yTM"
# tcol <- "Surveyyear"
# wtscol <- "yearind"
# 
# 
# ## 2004 through 2012 data
# moose_df <- readr::read_csv("inst/moose_04_12_all.csv")
# 
# tm_obj <- sptotal::LLtoTM(mean(moose_df$centrlon),
#                           lat = moose_df$centrlat,
#                           lon = moose_df$centrlon)
# moose_df$xTM <- tm_obj$xy[ ,1]
# moose_df$yTM <- tm_obj$xy[ ,2]
# 
# moose_df <- moose_df %>%
#   mutate(yearind = if_else(Surveyyear == 2012, true = 1, false = 0))
# 
# formula <- totalmoosena ~ 1
# data <- moose_df2
# xcoordcol <- "xTM"
# ycoordcol <- "yTM"
# tcol <- "Surveyyear"
# wtscol <- "yearind"

stlmfit <- function(formula, data, xcoordcol, ycoordcol, tcol,
                    wtscol, areacol = NULL,
                    CorModel = "Exponential") {
  
  ## order the data so that sites are in the same order within each time
  
  data <- data %>% arrange_(tcol, xcoordcol, ycoordcol) 
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
  ind.un <- is.na(yvar)
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
  
  ##  make sure the function is working
  ##  can replace this step with a grid search to give `parmest` 
  ##  a better starting point than (1, 1, 1, 1).
  m2LL.spatiotemp.ML(theta = c(0, 0, 0, 1, 0, 0, 0),
                     zcol = density,
                     XDesign = as.matrix(X),
                     xcoord = uniquecoords[ ,1],
                     ycoord = uniquecoords[ ,2],
                     timepoints = uniquetimes,
                     CorModel = "Exponential",
                     Zs = Zs, Zt = Zt, H = H, Dismat = Dismat)
  
  parmest <- stats::optim(c(0, 0, 0, 1, 0, 0, 0), m2LL.spatiotemp.ML,
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

  

  ###########################################################
  ## break this next part out into its own predict() function
  ###########################################################

  
  Sigma.ss <- Sigmaest[ind.sa, ind.sa, drop = FALSE]
  Sigma.us <- Sigmaest[ind.un, ind.sa, drop = FALSE]
  Sigma.su <- t(Sigma.us)
  Sigma.uu <- Sigmaest[ind.un, ind.un, drop = FALSE]
  
  Sigma.ssi <- solve(Sigma.ss)
  
  ## indicator for unsampled sites that are in the current year of 
  ## interest
  unsampcurrind <- ind.sa == 0 & data[[wtscol]] == 1
  
  Xucurr <- Xall[unsampcurrind, , drop = FALSE]
  
  betahat <- solve(t(Xs) %*% Sigma.ssi %*% Xs) %*%
    t(Xs) %*% Sigma.ssi %*%
    z.sa
  
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
  
  ## prediction for the total in the current year
  totalpred <- tlambda %*% z.density
  
  ## prediction weights for sites in the current year (usually just a 
  ## vector of 1's if we want to predict the total for the current year).
  bc <- c(bs, bu)
  
  predvar <- tlambda %*% Sigma.ss %*% t(tlambda) -
    2 * t(bc) %*% Sigma.cs %*% t(tlambda) +
    t(bc) %*% Sigma.cc %*% bc
  
  lb <- as.vector(totalpred) + c(-1) * 1.645 * sqrt(as.vector(predvar))
  ub <- as.vector(totalpred) + c(1) * 1.645 * sqrt(as.vector(predvar))
  
  
  ## an equivalent calculation for the total:
  muhats <- Xs %*% betahat; muhatu <- Xucurr %*% betahat
  ## the predicted values for the sites that were not sampled
  zhatu <- Sigma.ucurrs %*% Sigma.ssi %*% (z.density -
                                             muhats) + muhatu
  totalpred_equiv <- sum(zhatu) + sum(density[ind.sa == 1 & data[[wtscol]] == 1])
  
  W <- t(Xu) - t(Xs) %*% Sigma.ssi %*% Sigma.su
  Vmat <- solve(t(Xs) %*% Sigma.ssi %*% Xs)
  sitecov <- Sigma.uu - Sigma.us %*% Sigma.ssi %*% Sigma.su +
    t(W) %*% Vmat %*% W
  sitevar <- diag(sitecov)
  
  
  data$predictions <- rep(NA, nrow(data))
  
  data$predictions[unsampcurrind == 1] <- zhatu
  data$predictions[ind.sa == 1] <- z.density
  
  return(list(totalpred, predvar, lb, ub, zhatu, sitevar, 
              list(sigma_parsil_spat = sigma_parsil_spat_hat,
                   range = range_hat,
                   sigma_nugget_spat = sigma_nugget_spat_hat,
                   sigma_parsil_time = sigma_parsil_time_hat,
                   rho = rhotime_hat,
                   sigma_nugget_time = sigma_nugget_time_hat,
                   sigma_nugget_spacetime = sigma_nugget_spacetime_hat)))
  
}
