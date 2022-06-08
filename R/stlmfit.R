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
                    CorModel = "exponential") {
  
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
  zdensity <- yvar / areavar
  
  ind_sa <- !is.na(yvar)
  data$ind_sa <- ind_sa
  ind_un <- is.na(yvar)
  data$ind_un <- ind_un
  
  zdensity_samp <- zdensity[ind_sa]
  
  
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
  
  
  
  totalvar <- var(zdensity_samp)
  
  s_de_initial <- totalvar / 10
  s_ie_initial <- totalvar / 10
  t_de_initial <- totalvar / 10
  t_ie_initial <- totalvar / 10
  st_de_initial <- totalvar / 10
  st_ie_initial <- totalvar / 10
  total_var_initial <- sum(
    s_de_initial, s_ie_initial, t_de_initial,
    t_ie_initial, st_de_initial, st_ie_initial
  )
  s_range_initial <- median(Dismat)
  t_range_initial <- median(H)
  
  siminitial <-  DumelleEtAl2021STLMM::make_covparam_object(
    s_de = s_de_initial, s_ie = s_ie_initial,
    t_ie = t_ie_initial, t_de = t_de_initial,
    st_de = st_de_initial, st_ie = st_ie_initial,
    s_range = s_range_initial, t_range = t_range_initial,
    stcov = "productsum", estmethod = "reml"
  )
  
  data_sa <- datanomiss[ind_sa, ]
  data_un <- datanomiss[ind_un, ]
  
  # data_sampled <- data %>% dplyr::filter(!is.na(yvar))

  fast_est <- DumelleEtAl2021STLMM::stlmm(
    formula = formula,
    data = data_sa,
    xcoord = xcoordcol,
    ycoord = ycoordcol,
    tcoord = tcol,
    stcov = "productsum",
    estmethod = "reml",
    s_cor = CorModel,
    t_cor = CorModel,
    initial = siminitial,
    condition = 1e-4
  )
  
  
  # varstart <- log(var(density, na.rm = TRUE) / 10)
  # rangestart <- median(Dismat) / 2
  # rhostart <- 0
  
  # m2LL.spatiotemp.ML(theta = c(varstart, rangestart, varstart,
  #                             varstart, rhostart, varstart, varstart),
  #                    zcol = density,
  #                    XDesign = as.matrix(X),
  #                    xcoord = uniquecoords[ ,1],
  #                    ycoord = uniquecoords[ ,2],
  #                    timepoints = uniquetimes,
  #                    CorModel = "Exponential",
  #                    Zs = Zs, Zt = Zt, H = H, Dismat = Dismat)

  
  # parmest <- stats::optim(c(varstart, rangestart, varstart,
  #                           varstart, rhostart, varstart, varstart),
  #                         m2LL.spatiotemp.ML,
  #                         zcol = density,
  #                         XDesign = as.matrix(X),
  #                         xcoord = uniquecoords[ ,1],
  #                         ycoord = uniquecoords[ ,2],
  #                         timepoints = uniquetimes,
  #                         CorModel = "Exponential",
  #                         Zs = Zs, Zt = Zt, H = H, Dismat = Dismat)
  
  parmest <- unclass(fast_est$CovarianceParameters)
  sigma_parsil_spat_hat <- parmest[1]
  range_hat <- parmest[7]
  sigma_nugget_spat_hat <- parmest[2]
  
  sigma_parsil_time_hat <- parmest[3]
  rhotime_hat <- parmest[8]
  sigma_nugget_time_hat  <- parmest[4]
  
  sigma_parsil_spacetime_hat <- parmest[5]
  sigma_nugget_spacetime_hat <- parmest[6]

  
  
  
  if (CorModel == "exponential") {
    Rs_hat <- exp(-(3 * (Dismat / range_hat)))
  }
  
  comp_1 <- sigma_parsil_spat_hat * Zs %*% Rs_hat %*% t(Zs)
  comp_2 <- sigma_nugget_spat_hat * Zs %*% t(Zs)
  
  if (CorModel == "exponential") {
    Rt_hat <- exp(-(3 * (H / rhotime_hat)))
  }

  comp_3 <- sigma_parsil_time_hat * Zt %*% Rt_hat %*% t(Zt)
  comp_4 <- sigma_nugget_time_hat * Zt %*% t(Zt)
  
  comp_5 <- diag(sigma_nugget_spacetime_hat, nrow = N)
  
  comp_6 <- sigma_parsil_spacetime_hat * (Zs %*% Rs_hat %*% t(Zs)) * (Zt %*% Rt_hat %*% t(Zt))
  
  Sigmaest <- comp_1 + comp_2 + comp_3 + comp_4 + comp_5 + comp_6
  
  Sigma.ss <- Sigmaest[ind_sa, ind_sa, drop = FALSE]
  Sigma.ssi <- solve(Sigma.ss) ## can output this so it doesn't
  ## also get inverted in predict.stlmfit()
  
  betahat <- fast_est$Coefficients
  
  stlmfit_obj <- list(data = data, Sigmaest = Sigmaest, Xall = Xall, 
                      z.density = zdensity_samp, parms = list(betahat = betahat, sigma_parsil_spat = sigma_parsil_spat_hat,
                               range = range_hat,
                               sigma_nugget_spat = sigma_nugget_spat_hat,
                               sigma_parsil_time = sigma_parsil_time_hat,
                               rho = rhotime_hat,
                               sigma_nugget_time = sigma_nugget_time_hat,
                               sigma_nugget_spacetime = sigma_nugget_spacetime_hat,
                               sigma_parsil_spacetime = sigma_parsil_spacetime_hat))
  
  class(stlmfit_obj) <- "stlmfit"
  return(stlmfit_obj)
}
