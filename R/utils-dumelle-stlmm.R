## this file contains .R files from https://github.com/michaeldumelle/DumelleEtAl2021STLMM/tree/main/R

#' Fit a Spatio-Temporal Linear Mixed Model
#'
#' @param data A data object containing all necessary variables.
#'
#' @param formula A formula of the form \code{y ~ x}, where \code{y} is the response variable
#'   and \code{x} are the predictor variables.
#'
#' @param xcoord A character vector specifying the column name of the x-coordinate
#'   variable in \code{data}.
#'
#' @param ycoord A character vector specifying the column name of the y-coordinate
#'   variable in \code{data}.
#'
#' @param tcoord A character vector specifying the column name of the t-coordinate (time)
#'   variable in \code{data}.
#'
#' @param stcov The spatio-temporal covariance type
#'  \describe{
#'    \item{\code{product}}{The product LMM}
#'    \item{\code{sum_with_error}}{The sum-with-error LMM}
#'    \item{\code{productsum}}{The product sum LMM}
#'  }
#'
#' @param estmethod The estimation method
#'  \describe{
#'    \item{\code{reml}}{Restricted Maximum Likelihood}
#'    \item{\code{svwls}}{Semivariogram Weighted Least Squares}
#'  }
#'
#' @param s_cor The spatial correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'   }
#'
#' @param t_cor The temporal correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'     \item{\code{tent}}{The tent (linear with sill) correlation.}
#'   }
#'
#' @param chol Should the Cholesky decomposition be used? If \code{FALSE},
#'   efficient inversion algorithms are implemented. Defaults to \code{FALSE}.
#'
#' @param condition A small number added to the diagonals of matrices before
#'   inverting them to prevent ill-conditioning (defaults to \code{1e-4}).
#'
#' @param logdet Should the log determinant be returned? (defaults to \code{FALSE}).
#'
#' @param weights Weights when \code{estmethod = "svwls"} (defaults to
#'   \code{"cressie"}) for the Cressie's weighted least squares weights).
#'
#' @param initial Initial values for the parameters. Must be made with
#'   \code{make_covparam_object()} (defaults to even spread of across variance
#'   parameters, a total variance matching the sample variance of OLS residuals,
#'   and ranges equaling half the maximum distance in the domain).
#'
#' @param optim_options A list containing additional options to pass to
#'   \link[stats]{optim}.
#'
#' @param h_options A list containing options to compute distances if
#'   \code{response}, \code{xcoord}, \code{ycoord}, and \code{tcoord} are
#'   provided. Named arguments are
#'   \describe{
#'     \item{\code{h_t_distmetric}}{The temporal distance matrix (defaults to
#'     \code{"euclidean"}).}
#'     \item{\code{h_s_distmetric}}{The spatial distance matrix (defaults to
#'     \code{"euclidean"}).}
#'  }
#'
#' @param max_options A list containing additonal options for placing upper
#' bounds on the total variance, spatial range, and temporal range. This can
#' be helpful for numerical stability in optimization.
#'   \describe{
#'     \item{\code{max_v}}{The maximum total variance (defaults to four times the
#'      variance of OLS residuals)}
#'     \item{\code{max_s_range}}{The maximum spatial range variance (defaults to
#'     four times the maximum observed spatial distance)}
#'     \item{\code{max_t_range}}{The maximum temporal range variance (defaults to
#'     four times the maximum observed temporal distance)}
#'   }
#'
#' @param stempsv_options A list containing additional options for the
#'   empirical spatio-temporal semivariogram. Named arguments are
#'   \describe{
#'     \item{\code{n_s_lag}}{The number of spatial distance classes (defaults to 16).}
#'     \item{\code{n_t_lag}}{The number of temporal distance classes (defaults to 16).}
#'     \item{\code{h_s_max}}{The maximum spatial distance. Deafaults to half the
#'       maximum distance in the spatial domain.}
#'     \item{\code{h_t_max}}{The maximum temporal distance. Deafaults to half the
#'       maximum distance in the temporal domain.}
#'   }
#'
#' @param ... Additonal arguments.
#'
#' @return A list containing several objects
#'   \describe{
#'     \item{\code{CovarianceParameters}}{Estimated covariance parameters.}
#'     \item{\code{Coefficients}}{Fixed effect estimates.}
#'     \item{\code{NamesCoefficients}}{Names of the fixed effect estimates.}
#'     \item{\code{CovCoefficients}}{The covariance matrix of the fixed effect estimates.}
#'     \item{\code{Objective}}{A list containing optimization information.}
#'     \item{\code{CovarianceForms}}{The spatial, temopral, and spatio-temporal correlation forms.}
#'     \item{\code{formula}}{The model formula.}
#'     \item{\code{model}}{A list containing the fixed effect design matrix and response vector.}
#'     \item{\code{data_object}}{An ordered data object.}
#'     \item{\code{invert_object}}{An inverse object.}
#'     \item{\code{coord_names}}{The names of the coordinate vectors.}
#'     \item{\code{coords}}{The coordinate vectors.}
#'     \item{\code{h_options}}{Returning the \code{h_options} argument.}
#'     \item{\code{stempsv_options}}{Returning the \code{stempsv_options} argument.}
#'     \item{\code{stempsv}}{The empirical spatio-temporal semivariogram (if \code{estmethod = "svwls"})}
#'     \item{\code{optim_options}}{Returning the \code{stempsv_options} argument.}
#'     \item{\code{max_options}}{Returning the \code{max_options} argument.}
#'     \item{\code{chol}}{Returning the \code{chol} argument.}
#'     \item{\code{condition}}{Returning the \code{condition} argument.}
#'     \item{\code{residuals}}{Raw residuals.}
#'   }
#' @noRd


stlmm <- function(data, formula, ...) {
  UseMethod("stlmm", object = data)
}

#' @name stlmm
#' @method stlmm data.frame
#' @noRd

stlmm.data.frame <- function(data, formula, xcoord, ycoord = NULL, tcoord, stcov,
                             estmethod = "reml", s_cor = "exponential", t_cor = "exponential", chol = FALSE, condition = 1e-4,
                             logdet = FALSE, weights = "cressie", initial = NULL,
                             optim_options = NULL, h_options = NULL,
                             max_options = NULL, stempsv_options = NULL, ...) {
  
  # create the data object
  data_object <- make_data_object(
    formula = formula,
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    data = data,
    h_options = h_options
  )
  
  # create the covest object
  covest_object <- make_covest_object(
    initial = initial,
    estmethod = estmethod,
    stcov = stcov,
    data_object = data_object,
    condition = condition,
    chol = chol,
    s_cor = s_cor,
    t_cor = t_cor,
    weights = weights,
    max_options = max_options,
    optim_options = optim_options,
    stempsv_options = stempsv_options
  )
  
  
  
  
  # estimate the profiled covariance parameters
  covest_output <- covest_wrapper(covest_object = covest_object, data_object = data_object)
  
  # give the estimate parameters the same class as the original covest_object
  class(covest_output$par_r) <- class(covest_object)
  
  # create the invert object
  invert_object <- make_invert_object(
    covparam_object = covest_output$par_r,
    chol = chol,
    condition = condition,
    h_s_large = data_object$h_s_large,
    h_t_large = data_object$h_t_large,
    h_s_small = data_object$h_s_small,
    h_t_small = data_object$h_t_small,
    logdet = logdet,
    m_index = data_object$m_index,
    o_index = data_object$o_index,
    s_cor = s_cor,
    t_cor = t_cor,
    xo = data_object$ordered_xo,
    yo = data_object$ordered_yo
  )
  
  # compute the inverse covariance matrix (times the design matrix and response)
  invert_output <- invert(invert_object)
  
  
  # estimate the fixed effects
  betaest_output <- betaest(
    xo = data_object$ordered_xo,
    sigmainv_xyo = invert_output$sigmainv_o,
    condition = condition,
    return_estlist = FALSE # don't need the extra output required for reml optimization
  )
  
  # return the relevant summary output
  stlmm_object <- structure(
    list(
      CovarianceParameters = covest_output$par_r,
      Coefficients = betaest_output$betahat,
      NamesCoefficients = colnames(data_object$original_xo),
      CovCoefficients = betaest_output$cov_betahat,
      Objective = c(
        value = covest_output$value,
        counts = covest_output$count,
        convergence = covest_output$convergence,
        stempsv_seconds = covest_object$stempsv_seconds,
        optim_seconds = covest_output$optim_seconds
      ),
      CovarianceForms = c(
        stcov = stcov,
        s_cor = s_cor,
        t_cor = t_cor
      ),
      formula = formula,
      model = list(
        FixedDesignMatrix = data_object$original_xo,
        Response = data_object$original_yo
      ),
      data_object = data_object,
      invert_output = invert_output,
      coordnames = list(
        xcoord = xcoord,
        ycoord = ycoord,
        tcoord = tcoord
      ),
      coords = list(
        xcoord = data_object$ordered_data_o[[xcoord]],
        ycoord = data_object$ordered_data_o[[ycoord]],
        tcoord = data_object$ordered_data_o[[tcoord]]
      ),
      h_options = data_object$h_options,
      stempsv_options = covest_object$stempsv_options,
      stempsv = covest_object$stempsv,
      optim_options = covest_object$optim_options,
      max_options = covest_object$max_options,
      chol = chol,
      condition = condition
    ),
    class = "stlmm" # give the object class "stlmm"
  )
  
  # compute the residuals
  stlmm_object$Residuals <- residuals(object = stlmm_object)
  
  # return the stlmm_object
  return(stlmm_object)
}

betaest <- function(xo, sigmainv_xyo, condition, return_estlist = FALSE) {
  
  # number of columns of X (of the observed data)
  ncol_xo <- ncol(xo)
  
  # a sequence from 1 to the number of columns in X
  xo_dims <- seq.int(1, ncol_xo)
  
  # Taking sigma_inverse %*% X
  sigmainv_xo <- sigmainv_xyo[, xo_dims, drop = FALSE]
  
  # Taking sigma_inverse %*% Y
  sigmainv_yo <- sigmainv_xyo[, ncol_xo + 1, drop = FALSE]
  
  # Inverse of cov beta hat
  invcov_betahat <- t(xo) %*% sigmainv_xo
  
  # Adding diagonal stability
  diag(invcov_betahat) <- diag(invcov_betahat) + condition
  
  # The cholesky of this matrix
  chol_invcov_betahat <- chol(invcov_betahat)
  
  # Finding cov beta hat
  cov_betahat <- chol2inv(chol_invcov_betahat)
  
  # Computing beta hat
  betahat <- cov_betahat %*% t(xo) %*% sigmainv_yo
  
  # Returning betahat and the covariance in a list
  betaest_output <- list(betahat = betahat, cov_betahat = cov_betahat)
  
  # returning relevant output necessary for estimation
  if (return_estlist) {
    betaest_output$estlist <- list(
      # the log determinant of the inverse of cov beta hat
      ldet_cicb = 2 * sum(log(diag(chol_invcov_betahat))),
      sigmainv_xo = sigmainv_xo,
      sigmainv_yo = sigmainv_yo,
      p = ncol_xo
    )
  }
  return(betaest_output)
}

covest <- function(par, covest_object, ...) {
  UseMethod("covest", object = covest_object)
}

# semivariogram weighted least squares optimization
covest.svwls <- function(par, covest_object, data_object) {
  
  # transform profiled variance parameters to regular
  plo2r <- plo2r.svwls(par, covest_object)
  
  # copy the class of covest_object to us the appropriate generic
  class(plo2r) <- class(covest_object)
  
  # make the spatio-temopral semivariogram
  theo_sv <- make_stsemivariogram(
    covparam_object = plo2r,
    h_s_large = covest_object$stempsv$h_s_avg,
    h_t_large = covest_object$stempsv$h_t_avg,
    s_cor = covest_object$s_cor,
    t_cor = covest_object$t_cor
  )
  # create the weights used in the weighted least squares optimization
  wts <- switch(
    covest_object$weights,
    "cressie" = weights_cressie(sv = covest_object$stempsv, theo_sv = theo_sv),
    stop("choose valid weights")
  )
  
  # create the objective function
  sumsq <- (covest_object$stempsv$gammahat - theo_sv)^2
  
  # return the objective function
  return(sum(wts * sumsq))
}

# create cressie (1985) weights
weights_cressie <- function(sv, theo_sv) {
  wts <- sv$n / theo_sv^2
  return(wts)
}

# reml optimization
covest.reml <- function(par, covest_object, invert_object) {
  
  # transform profiled variance parameters to regular
  ## the overall variance is 1 because it has been profiled
  plo2r <- plo2r.reml(par, covest_object, ov_var = 1)
  
  # change the covariance parameter vecotr in invert_object
  invert_object$covparams <- plo2r
  
  # compute the inverse
  invert_output <- invert(invert_object)
  
  # compute minus twice the negative log likelihood
  m2ll <- minus2loglik.reml(invert_object = invert_object, invert_output = invert_output)
  
  # return minus twice the negative log likelihood
  return(m2ll)
}

covest_wrapper <- function(covest_object, data_object) {
  UseMethod("covest_wrapper", object = covest_object)
}

# wrapper for the weighted least squarse optimization
covest_wrapper.svwls <- function(covest_object, data_object) {
  
  # performing the optimization
  
  # timing
  optim_start <- Sys.time()
  covest_output <- optim(
    par = covest_object$initial_plo,
    fn = covest.svwls,
    covest_object = covest_object,
    data_object = data_object,
    method = covest_object$optim_options$method,
    control = covest_object$optim_options$control
  )
  optim_end <- Sys.time()
  covest_output$optim_seconds <- as.numeric(optim_end - optim_start, units = "secs")
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  # covest_output$par[covest_output$par > 7] <- 7
  # covest_output$par[covest_output$par < -7] <- -7
  
  # transforming the optimized parameter values to regular values
  covest_output$par_r <- plo2r.svwls(par = covest_output$par, covest_object = covest_object)
  
  # a warning if convergence did not occur
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }
  
  # returning the covariance parameter output
  return(covest_output)
}

covest_wrapper.reml <- function(covest_object, data_object) {
  
  # newly storing the initial values (profiled log odds)
  initial_plo_noclass <- covest_object$initial_plo
  
  # giving this the same class as the covest_object
  class(covest_object$initial_plo) <- class(covest_object)
  
  # making the inverse object - many are from the data object
  invert_object <- make_invert_object(
    covparam_object = covest_object$initial_plo,
    chol = covest_object$chol,
    co = NULL,
    condition = covest_object$condition,
    h_s_large = data_object$h_s_large,
    h_t_large = data_object$h_t_large,
    h_s_small = data_object$h_s_small,
    h_t_small = data_object$h_t_small,
    logdet = covest_object$logdet,
    m_index = data_object$m_index,
    o_index = data_object$o_index,
    s_cor = covest_object$s_cor,
    t_cor = covest_object$t_cor,
    xo = data_object$ordered_xo,
    yo = data_object$ordered_yo
  )
  
  # performing the optimization
  
  # timing
  optim_start <- Sys.time()
  covest_output <- optim(
    par = initial_plo_noclass,
    fn = covest.reml,
    covest_object = covest_object,
    invert_object = invert_object,
    method = covest_object$optim_options$method,
    control = covest_object$optim_options$control
  )
  optim_end <- Sys.time()
  covest_output$optim_seconds <- as.numeric(optim_end - optim_start, units = "secs")
  
  
  # saving the profiled covariance parameter output
  invert_object$covparams <- plo2r.reml(covest_output$par, covest_object = covest_object, ov_var = 1)
  
  # computing the inverse
  invert_output <- invert(invert_object)
  
  # computing the overall variance
  ov_var <- varest.reml(invert_object = invert_object, invert_output = invert_output)
  
  # storing the regular covariance parameter output
  covest_output$par_r <- plo2r.reml(covest_output$par, covest_object = covest_object, ov_var = ov_var)
  
  # a warning if convergence did not occur
  if (covest_output$convergence != 0) {
    warning("covariance parameter convergence may not have been achieved - consider
            setting new initial values, lowering the relative tolerance, or increasing
            the maximum iterations")
  }
  
  # returning the covariance parameter output
  return(covest_output)
}

#' Create Initial Values
#'
#' A wrapper around \code{make_covparam_object} that requires specification
#' of \code{estmethod}.
#'
#' @inheritParams make_covparam_object
#'
#' @return A named vector with covariance parameters having class equal to
#' the \code{estmethod} argument and the \code{stcov} argument.
#' @noRd

initial <- function(s_de,
                    s_ie,
                    t_de,
                    t_ie,
                    st_de,
                    st_ie,
                    v_s,
                    v_t,
                    s_range,
                    t_range,
                    estmethod,
                    stcov) {
  
  # make the initial object covariance parameter object
  # it is the same as make_covparam_object() but forces
  # the user to specify an estimation method
  initialcov <- make_covparam_object(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie,
    v_s = v_s,
    v_t = v_t,
    s_range = s_range,
    t_range = t_range,
    estmethod = estmethod,
    stcov = stcov
  )
  
  # returning the initial covaiance parameter vector
  return(initialcov)
}

#' Covariance Matrix Inversion
#'
#' @param invert_object An inversion object made from \code{make_invert_object}.
#'
#' @return Relevant covariance matrix inversion information.
#' @noRd

invert <- function(invert_object) {
  UseMethod("invert", object = invert_object)
}

#' @name invert
#'
#' @method invert productsum
#' @noRd

invert.productsum <- function(invert_object) {
  # invert a product sum covariance matrix
  # invert using a the standard cholesky decomposition approach
  if (invert_object$chol) {
    
    # make the covariance matrix
    sigma <- make_stcovariance.productsum(
      covparam_object = invert_object$covparams,
      h_s_large = invert_object$h_s_large,
      h_t_large = invert_object$h_t_large,
      s_cor = invert_object$s_cor,
      t_cor = invert_object$t_cor
    )
    
    # adding condition number stability
    diag(sigma) <- diag(sigma) + invert_object$condition
    
    # finding the Cholesky decomposition
    chol_sigma <- chol(sigma)
    
    # computing the inverse
    siginv <- chol2inv(chol_sigma)
    
    # multiplying this inverse on the right
    siginv_o <- siginv %*% invert_object$xyc_o
    
    # computing the log determinant if requested
    if (invert_object$logdet) {
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {
    
    # saving the required objects to invert the matrix - this needs to be cleaned up
    # at some point
    r_s_small <- make_r(
      h = invert_object$h_s_small,
      range = invert_object$covparams[["s_range"]],
      structure = invert_object$s_cor
    )
    r_t_small <- make_r(
      h = invert_object$h_t_small,
      range = invert_object$covparams[["t_range"]],
      structure = invert_object$t_cor
    )
    s_de <- invert_object$covparams[["s_de"]]
    s_ie <- invert_object$covparams[["s_ie"]]
    t_de <- invert_object$covparams[["t_de"]]
    t_ie <- invert_object$covparams[["t_ie"]]
    st_de <- invert_object$covparams[["st_de"]]
    st_ie <- invert_object$covparams[["st_ie"]]
    xyc_o <- invert_object$xyc_o
    condition <- invert_object$condition
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index
    n_s <- invert_object$n_s
    n_t <- invert_object$n_t
    n_st <- n_s * n_t
    
    dense <- length(o_index) == n_st
    
    
    
    
    ## adding diagonal tolerances for invertibility stability - for the correlation matrices, they are
    ## also rescaled so the diagonal is 1
    r_s_small <- r_s_small / (1 + condition)
    diag(r_s_small) <- 1
    r_t_small <- r_t_small / (1 + condition)
    diag(r_t_small) <- 1
    s_ie <- s_ie + condition
    t_ie <- t_ie + condition
    st_ie <- st_ie + condition
    
    ## finding eigendecompositions
    r_s_small_eigen <- eigen(r_s_small)
    r_t_small_eigen <- eigen(r_t_small)
    
    ## creating w matrix
    w <- kronecker(r_t_small_eigen$vectors, r_s_small_eigen$vectors)
    # creating v matrix
    v <- st_de * kronecker(r_t_small_eigen$values, r_s_small_eigen$values) + st_ie
    # creating inverse square root of v matrix
    vinvroot <- 1 / sqrt(v)
    # creating inverse of w * sqrt(v)
    t_w_vinvroot <- as.vector(vinvroot) * t(w)
    # creating the transpose
    w_vinvroot <- t(t_w_vinvroot)
    
    
    # storing the output we will need for the iterative smw
    c_t <- chol(make_sigma(de = t_de, r_mx = r_t_small, ie = t_ie))
    c_s <- chol(make_sigma(de = s_de, r_mx = r_s_small, ie = s_ie))
    
    ist_zt <- w_vinvroot %*% multiply_z(t_w_vinvroot, "temporal", n_s, n_t, "right")
    # st x st %*% (st x st * st x t) = st x t
    tr_ist_zt <- t(ist_zt)
    c_mt <- chol(chol2inv(c_t) + multiply_z(ist_zt, "temporal", n_s, n_t, "p_left"))
    # t(multiply_z(mx = t(ist_zt), z_type = "temporal", n_s = n_s)))
    # (t x t + tr(tr(st x t) * st x t) = t x t
    ic_mt <- chol2inv(c_mt)
    # t x t
    
    istpt_zs <- w_vinvroot %*% multiply_z(mx = t_w_vinvroot, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right") -
      ist_zt %*% (ic_mt %*% multiply_z(mx = tr_ist_zt, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right"))
    # st x st * (st x st * st x s) - st x t (t x t * (tr(st x t) * st x s)) =
    # st x s - st x t * t x s = st x s
    tr_istpt_zs <- t(istpt_zs)
    c_ms <- chol(chol2inv(c_s) + multiply_z(mx = istpt_zs, z_type = "spatial", n_s = n_s, n_t = n_t, side = "p_left"))
    # t(multiply_z(mx = tr_istpt_zs, z_type = "spatial", n_s = n_s)))
    # s x s + tr(tr(st x s) * st x s) = s x s
    ic_ms <- chol2inv(c_ms)
    # s x s
    
    # now to implement algorithm
    
    if (dense) {
      siginv_o <- w_vinvroot %*% (t_w_vinvroot %*% xyc_o) -
        # st x st * (st x st * st x p) = st x p
        ist_zt %*% (ic_mt %*% (tr_ist_zt %*% xyc_o)) -
        # st x t * (t x t * (t x st * st x p)) = st x p
        istpt_zs %*% (ic_ms %*% (tr_istpt_zs %*% xyc_o))
      # st x s * (s x s * (s x st * st x p)) = st x p
    } else {
      d_oo <- w_vinvroot[o_index, o_index, drop = FALSE] %*% (t_w_vinvroot[o_index, o_index, drop = FALSE] %*% xyc_o) +
        w_vinvroot[o_index, m_index, drop = FALSE] %*% (t_w_vinvroot[m_index, o_index, drop = FALSE] %*% xyc_o) -
        ist_zt[o_index, , drop = FALSE] %*%
        (ic_mt %*% (tr_ist_zt[, o_index, drop = FALSE] %*% xyc_o)) -
        istpt_zs[o_index, , drop = FALSE] %*%
        (ic_ms %*% (tr_istpt_zs[, o_index, drop = FALSE] %*% xyc_o))
      
      d_om <- w_vinvroot[o_index, o_index, drop = FALSE] %*% t_w_vinvroot[o_index, m_index, drop = FALSE] +
        w_vinvroot[o_index, m_index, drop = FALSE] %*% t_w_vinvroot[m_index, m_index, drop = FALSE] -
        ist_zt[o_index, , drop = FALSE] %*%
        (ic_mt %*% tr_ist_zt[, m_index, drop = FALSE]) -
        istpt_zs[o_index, , drop = FALSE] %*%
        (ic_ms %*% tr_istpt_zs[, m_index, drop = FALSE])
      
      d_mm <- w_vinvroot[m_index, o_index, drop = FALSE] %*% t_w_vinvroot[o_index, m_index, drop = FALSE] +
        w_vinvroot[m_index, m_index, drop = FALSE] %*% t_w_vinvroot[m_index, m_index, drop = FALSE] -
        ist_zt[m_index, , drop = FALSE] %*%
        (ic_mt %*% tr_ist_zt[, m_index, drop = FALSE]) -
        istpt_zs[m_index, , drop = FALSE] %*%
        (ic_ms %*% tr_istpt_zs[, m_index, drop = FALSE])
      
      # return the correct object
      c_mm <- chol(d_mm)
      siginv_o <- d_oo - d_om %*% (chol2inv(c_mm) %*% (t(d_om) %*% xyc_o))
    }
    
    if (logdet) {
      logdet <- sum(log(v)) +
        2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
        2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))
      
      if (!dense) {
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
      }
    } else {
      logdet <- NULL
    }
  }
  
  output <- list(sigmainv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

#' @name invert
#'
#' @method invert sum_with_error
#' @noRd

invert.sum_with_error <- function(invert_object) {
  if (invert_object$chol) {
    
    # make the covariance matrix
    sigma <- make_stcovariance.sum_with_error(
      covparam_object = invert_object$covparams,
      h_s_large = invert_object$h_s_large,
      h_t_large = invert_object$h_t_large,
      s_cor = invert_object$s_cor,
      t_cor = invert_object$t_cor
    )
    
    # adding condition number stability
    diag(sigma) <- diag(sigma) + invert_object$condition
    
    # finding the Cholesky decomposition
    chol_sigma <- chol(sigma)
    
    # computing the inverse
    siginv <- chol2inv(chol_sigma)
    
    # multiplying this inverse on the right
    siginv_o <- siginv %*% invert_object$xyc_o
    
    # computing the log determinant if requested
    if (invert_object$logdet) {
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {
    
    # saving the required objects to invert the matrix - this needs to be cleaned up
    # at some point
    r_s_small <- make_r(
      h = invert_object$h_s_small,
      range = invert_object$covparams[["s_range"]],
      structure = invert_object$s_cor
    )
    r_t_small <- make_r(
      h = invert_object$h_t_small,
      range = invert_object$covparams[["t_range"]],
      structure = invert_object$t_cor
    )
    s_de <- invert_object$covparams[["s_de"]]
    s_ie <- invert_object$covparams[["s_ie"]]
    t_de <- invert_object$covparams[["t_de"]]
    t_ie <- invert_object$covparams[["t_ie"]]
    st_ie <- invert_object$covparams[["st_ie"]]
    xyc_o <- invert_object$xyc_o
    condition <- invert_object$condition
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index
    n_s <- invert_object$n_s
    n_t <- invert_object$n_t
    n_st <- n_s * n_t
    
    
    dense <- length(o_index) == n_st
    
    
    
    
    ## adding diagonal tolerances for invertibility stability - for the correlation matrices, they are
    ## also rescaled so the diagonal is 1
    r_s_small <- r_s_small / (1 + condition)
    diag(r_s_small) <- 1
    r_t_small <- r_t_small / (1 + condition)
    diag(r_t_small) <- 1
    s_ie <- s_ie + condition
    t_ie <- t_ie + condition
    st_ie <- st_ie + condition
    
    # storing the output we will need for the iterative smw
    c_t <- chol(make_sigma(de = t_de, r_mx = r_t_small, ie = t_ie))
    c_s <- chol(make_sigma(de = s_de, r_mx = r_s_small, ie = s_ie))
    
    c_mt <- chol(chol2inv(c_t) + multiply_z(z_type = "temporal", n_s = n_s, n_t = n_t, side = "pz_z") / st_ie)
    ic_mt <- chol2inv(c_mt)
    istpt <- -multiply_z(multiply_z(ic_mt, z_type = "temporal", n_s = n_s, n_t = n_t, side = "p_right"),
                         z_type = "temporal", n_s = n_s, n_t = n_t, side = "left"
    ) / (st_ie^2)
    diag(istpt) <- diag(istpt) + 1 / st_ie
    
    
    istpt_zs <- multiply_z(mx = istpt, z_type = "spatial", n_s = n_s, n_t = n_t, side = "right")
    c_ms <- chol(chol2inv(c_s) + multiply_z(mx = istpt_zs, z_type = "spatial", n_s = n_s, n_t = n_t, side = "p_left"))
    ic_ms <- chol2inv(c_ms)
    
    
    siginv <- istpt - istpt_zs %*% (ic_ms %*% t(istpt_zs))
    
    if (dense) {
      siginv_o <- siginv %*% xyc_o
    } else {
      c_mm <- chol(siginv[m_index, m_index])
      siginv_o <- (siginv[o_index, o_index] - siginv[o_index, m_index] %*% (chol2inv(c_mm) %*% siginv[m_index, o_index])) %*% xyc_o
    }
    
    
    if (logdet) {
      logdet <- n_st * (log(st_ie)) +
        2 * sum(log(diag(c_t))) + 2 * sum(log(diag(c_mt))) +
        2 * sum(log(diag(c_s))) + 2 * sum(log(diag(c_ms)))
      
      if (!dense) {
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
      }
    } else {
      logdet <- NULL
    }
  }
  
  output <- list(sigmainv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}






#' @name invert
#'
#' @method invert product
#' @noRd

invert.product <- function(invert_object) {
  if (invert_object$chol) {
    
    # make the covariance matrix
    sigma <- make_stcovariance.product(
      covparam_object = invert_object$covparams,
      h_s_large = invert_object$h_s_large,
      h_t_large = invert_object$h_t_large,
      s_cor = invert_object$s_cor,
      t_cor = invert_object$t_cor
    )
    
    # adding condition number stability
    diag(sigma) <- diag(sigma) + invert_object$condition
    
    # finding the Cholesky decomposition
    chol_sigma <- chol(sigma)
    
    # computing the inverse
    siginv <- chol2inv(chol_sigma)
    
    # multiplying this inverse on the right
    siginv_o <- siginv %*% invert_object$xyc_o
    
    # computing the log determinant if requested
    if (invert_object$logdet) {
      logdet <- 2 * sum(log(diag(chol_sigma)))
    } else {
      logdet <- NULL
    }
  } else {
    
    # saving the required objects to invert the matrix - this needs to be cleaned up
    # at some point
    r_s_small <- make_r(
      h = invert_object$h_s_small,
      range = invert_object$covparams[["s_range"]],
      structure = invert_object$s_cor
    )
    r_t_small <- make_r(
      h = invert_object$h_t_small,
      range = invert_object$covparams[["t_range"]],
      structure = invert_object$t_cor
    )
    st_de <- invert_object$covparams[["st_de"]]
    v_s <- invert_object$covparams[["v_s"]]
    v_t <- invert_object$covparams[["v_t"]]
    xyc_o <- invert_object$xyc_o
    condition <- invert_object$condition
    logdet <- invert_object$logdet
    o_index <- invert_object$o_index
    m_index <- invert_object$m_index
    n_s <- invert_object$n_s
    n_t <- invert_object$n_t
    n_st <- n_s * n_t
    
    
    dense <- length(o_index) == n_st
    
    r_s_small <- r_s_small / (1 + condition)
    diag(r_s_small) <- 1
    r_t_small <- r_t_small / (1 + condition)
    diag(r_t_small) <- 1
    
    scale_r_s_small <- make_sigma(r_mx = r_s_small, v_ie = v_s, e = 1, scale = TRUE)
    c_scale_r_s_small <- chol(scale_r_s_small)
    scale_r_t_small <- make_sigma(r_mx = r_t_small, v_ie = v_t, e = 1, scale = TRUE)
    c_scale_r_t_small <- chol(scale_r_t_small)
    
    siginv <- kronecker(chol2inv(c_scale_r_t_small), chol2inv(c_scale_r_s_small)) / st_de
    
    if (dense) {
      siginv_o <- siginv %*% xyc_o
    } else {
      c_mm <- chol(siginv[m_index, m_index])
      siginv_o <- (siginv[o_index, o_index] - siginv[o_index, m_index] %*% (chol2inv(c_mm) %*% siginv[m_index, o_index])) %*% xyc_o
    }
    
    
    if (logdet) {
      logdet <- n_st * log(st_de) +
        n_s * 2 * sum(log(diag(c_scale_r_t_small))) +
        n_t * 2 * sum(log(diag(c_scale_r_s_small)))
      if (!dense) {
        logdet <- logdet + 2 * sum(log(diag(c_mm)))
      }
    } else {
      logdet <- NULL
    }
  }
  output <- list(sigmainv_o = siginv_o, logdet = logdet)
  output_non_null <- output[!unlist(lapply(output, is.null))]
  return(output_non_null)
}

make_covest_object <- function(initial = NULL,
                               estmethod,
                               stcov,
                               s_cor,
                               t_cor,
                               data_object,
                               weights = NULL,
                               chol = NULL,
                               condition = NULL,
                               max_options = NULL,
                               optim_options = NULL,
                               stempsv_options = NULL) {
  
  # setting default optim options
  if (is.null(optim_options)) {
    if (estmethod == "reml") {
      # lower tolerance for REML
      control <- list(reltol = 1e-4, maxit = 2000)
    } else {
      # higher tolerance for CWLS
      control <- list(reltol = 1e-8, maxit = 10000)
    }
    # Nelder-Mead as a default
    optim_options <- list(method = "Nelder-Mead", control = control)
  }
  if (is.null(max_options)) {
    # setting null values for the maxes - these NULLS are given to
    # appropriate defaults later
    max_options <- list(max_v = NULL, max_s_range = NULL, max_t_range = NULL)
  }
  
  if (is.null(stempsv_options)) {
    if (estmethod == "svwls") {
      # defaults for the semivariogram options - the NULLS are given to
      # appropriate defaults later
      stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
    } else {
      stempsv_options <- NULL
    }
  }
  
  # find the linear model residuals
  # beta ols
  beta_ols <- chol2inv(chol(t(data_object$ordered_xo) %*% data_object$ordered_xo)) %*%
    (t(data_object$ordered_xo) %*% data_object$ordered_yo)
  # computing the ols residuals
  lmod_r <- as.vector(data_object$ordered_yo - data_object$ordered_xo %*% beta_ols)
  
  # computing the ols sample variance
  lmod_s2 <- sum(lmod_r^2) / (nrow(data_object$ordered_xo) - ncol(data_object$ordered_xo))
  
  # a check to make sure that the large distance matrices are provided if required
  if (is.null(data_object$h_s_large) || is.null(data_object$h_t_large)) {
    stop("h_s_large and h_t_large must be non NULL in data_object: if using stlmm, set h_options$h_large = TRUE")
  }
  
  # setting initial values if none are provided (sample variance evenly
  # divided among all variance components)
  if (is.null(initial)) {
    initial <- initial(
      s_de = 1,
      s_ie = 1,
      t_de = 1,
      t_ie = 1,
      st_de = 1,
      st_ie = 1,
      v_s = 0.5,
      v_t = 0.5,
      s_range = max(data_object$h_s_small) / 2, # 2 chosen so that it is half the max observed distance
      t_range = max(data_object$h_t_small) / 2, # 2 chosen so that it is half the max observed distance
      estmethod = estmethod,
      stcov = stcov
    )
    
    # storing the parameter names
    vparm_names <- c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie")
    
    # summing the relevant parameters for each st covariance type
    numparams <- sum(vparm_names %in% names(initial))
    
    # scaling initial value variance parameters
    initial[names(initial) %in% vparm_names] <- lmod_s2 / numparams
  }
  
  # provide default value for the maximum possible variance
  if (is.null(max_options$max_v)) {
    # find the sample variance - this is needed if the initial values
    # are supplied by the user
    s2 <- sum(initial[!(names(initial) %in% c("s_range", "t_range"))])
    
    # setting a max variance equal to 4 times the sample variance
    max_options$max_v <- 4 * s2
  }
  
  # provide default value for the maximum possible spatial range
  if (is.null(max_options$max_s_range)) {
    max_options$max_s_range <- 4 * max(data_object$h_s_small)
  }
  
  # provide default value for the maximum possible temporal range
  if (is.null(max_options$max_t_range)) {
    max_options$max_t_range <- 4 * max(data_object$h_t_small)
  }
  
  # giving a warning message
  if (!identical(class(initial), c(estmethod, stcov))) {
    stop("class of initial value parameter vector must match estmethod and stcov used in stlmm")
  }
  
  # storing profiled log odds initial variance parameters for
  # unconstrained optimization
  initial_plo <- r2plo(covparam_object = initial, max_options = max_options)
  
  # making the semivariogram if required
  if (estmethod == "svwls") {
    
    # storing the squared difference of residuals
    # timing
    hresp_start <- Sys.time()
    h_response <- make_h(coord1 = lmod_r, distmetric = "euclidean")^2
    hresp_end <- Sys.time()
    hresp_seconds <- as.numeric(hresp_end - hresp_start, units = "secs")
    
    # making the empirical semivariogram using the distance matrices
    
    # timing
    stempsv_short_start <- Sys.time()
    stempsv <- stempsv(
      h_response = h_response,
      h_s_large = data_object$h_s_large,
      h_t_large = data_object$h_t_large,
      stempsv_options = stempsv_options
    )
    stempsv_short_end <- Sys.time()
    stempsv_short_seconds <- as.numeric(stempsv_short_end - stempsv_short_start, units = "secs")
    
    # computing the final seconds
    stempsv_seconds <- data_object$hdist_seconds + hresp_seconds + stempsv_short_seconds
    
    # passing through the weights
    weights <- weights
    
    # setting relevant values to NULL for this type of
    # covariance estimation (these are only set so that
    # we can use these function inputs which are general
    # for each type of covariance estimation method)
    chol <- NULL
    logdet <- NULL
    condition <- NULL
  } else {
    
    # passing along relevant arguments for REML estimation
    stempsv <- NULL
    weights <- NULL
    chol <- chol
    logdet <- TRUE
    condition <- condition
    stempsv_seconds <- 0
  }
  
  # making the covest_object and giving it the appropriate class
  covest_object <- structure(
    list(
      chol = chol,
      condition = condition,
      initial = initial,
      initial_plo = initial_plo,
      logdet = logdet,
      max_options = max_options,
      optim_options = optim_options,
      s_cor = s_cor,
      stempsv = stempsv,
      stempsv_options = stempsv_options,
      stempsv_seconds = stempsv_seconds,
      t_cor = t_cor,
      weights = weights
    ),
    class = class(initial)
  )
  
  # returning the appropriate covest_object
  return(covest_object)
}


#' Make a Covariance Parameter Object
#'
#' @param s_de The spatial dependent variance (spatial partial sill).
#'
#' @param s_ie The spatial independent variance (spatial nugget).
#'
#' @param t_de The temporal dependent variance (temporal partial sill).
#'
#' @param t_ie The temporal independent variance (temporal nugget).
#'
#' @param st_de The spatio-temporal dependent variance (spatio-temporal partial sill).
#'
#' @param st_ie The spatio-temporal independent variance (spatio-temporal nugget).
#'
#' @param v_s The proportion of spatial dependent variance
#'   (if \code{estmethod = "product"}).
#'
#' @param v_t The proportion of temporal dependent variance
#'   (if \code{estmethod = "product"}).
#'
#' @param s_range The spatial effective range (the spatial distance at which
#'   the correlation equals 0.05 (for non-compact) or 0 (for compact correlations)
#'
#' @param t_range The spatial effective range (the spatial distance at which
#'   the correlation equals 0.05 (for non-compact) or 0 (for compact correlations)
#'
#' @param estmethod The estimation method
#'  \describe{
#'    \item{\code{reml}}{Restricted Maximum Likelihood}
#'    \item{\code{svwls}}{Semivariogram Weighted Least Squares}
#'  }
#'
#' @param stcov The spatio-temporal covariance type
#'  \describe{
#'    \item{\code{product}}{The product LMM}
#'    \item{\code{sum_with_error}}{The sum-with-error LMM}
#'    \item{\code{productsum}}{The product sum LMM}
#'  }
#'
#' @return A named vector with covariance parameters having class equal to
#' the \code{estmethod} argument (if provided) and the \code{stcov} argument.
#'
#' @seealso [initial()]
#' @noRd

make_covparam_object <- function(s_de,
                                 s_ie,
                                 t_de,
                                 t_ie,
                                 st_de,
                                 st_ie,
                                 v_s,
                                 v_t,
                                 s_range,
                                 t_range,
                                 estmethod = NULL,
                                 stcov) {
  
  # conditional call for type of spatio-temporal covariance
  covparam_object <- switch(
    stcov,
    "productsum" = covparam_object_productsum(
      s_de = s_de,
      s_ie = s_ie,
      t_de = t_de,
      t_ie = t_ie,
      st_de = st_de,
      st_ie = st_ie,
      s_range = s_range,
      t_range = t_range
    ),
    "sum_with_error" = covparam_object_sum_with_error(
      s_de = s_de,
      s_ie = s_ie,
      t_de = t_de,
      t_ie = t_ie,
      st_ie = st_ie,
      s_range = s_range,
      t_range = t_range
    ),
    "product" = covparam_object_product(
      st_de = st_de,
      v_s = v_s,
      v_t = v_t,
      s_range = s_range,
      t_range = t_range
    ),
    stop("Use a valid error structure")
  )
  
  # giving the object the appropriate estimation method and stcov class
  covparam_object <- structure(covparam_object, class = c(estmethod, stcov))
  
  # returning the object
  return(covparam_object)
}

# make product sum covariance parameter object
covparam_object_productsum <- function(s_de,
                                       s_ie,
                                       t_de,
                                       t_ie,
                                       st_de,
                                       st_ie,
                                       s_range,
                                       t_range) {
  # create the covariance parameter vector
  cov_vec <- c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie,
    s_range = s_range,
    t_range = t_range
  )
  
  # return the covariance parameter vector
  return(cov_vec)
}

# make sum with error covariance parameter object
covparam_object_sum_with_error <- function(s_de,
                                           s_ie,
                                           t_de,
                                           t_ie,
                                           st_ie,
                                           s_range,
                                           t_range) {
  # create the covariance parameter vector
  cov_vec <- c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_ie = st_ie,
    s_range = s_range,
    t_range = t_range
  )
  
  # return the covariance parameter vector
  return(cov_vec)
}

# make product covariance parameter object
covparam_object_product <- function(st_de,
                                    v_s,
                                    v_t,
                                    s_range,
                                    t_range) {
  # create the covariance parameter vector
  cov_vec <- c(
    st_de = st_de,
    v_s = v_s,
    v_t = v_t,
    s_range = s_range,
    t_range = t_range
  )
  
  # return the covariance parameter vector
  return(cov_vec)
}

#' Make the Data Object
#'
#' @param formula A formula of the form \code{y ~ x}, where \code{y} is the response variable
#'   and \code{x} are the predictor variables.
#'
#' @param xcoord A character vector specifying the column name of the x-coordinate
#'   variable in \code{data}.
#'
#' @param ycoord A character vector specifying the column name of the y-coordinate
#'   variable in \code{data}.
#'
#' @param tcoord A character vector specifying the column name of the t-coordinate (time)
#'   variable in \code{data}.
#'
#' @param data A data object containing all necessary variables.
#'
#' @param h_options A list containing options to compute distances if
#'   \code{response}, \code{xcoord}, \code{ycoord}, and \code{tcoord} are
#'   provided. Named arguments are
#'   \describe{
#'     \item{\code{h_t_distmetric}}{The temporal distance matrix (defaults to
#'     \code{"euclidean"}).}
#'     \item{\code{h_s_distmetric}}{The spatial distance matrix (defaults to
#'     \code{"euclidean"}).}
#'  }
#'
#' @return A list with relevant data ordering information.
#' @noRd

make_data_object <- function(formula, xcoord, ycoord, tcoord, data, h_options) {
  
  # setting a default for h_options
  if (is.null(h_options)) {
    h_options <- list(
      h_large = TRUE,
      h_t_distmetric = "euclidean",
      h_s_distmetric = "euclidean"
    )
  }
  
  # restoring the original data
  original_data <- data
  
  # making the original model frame
  original_stmodel_frame <- model.frame(formula, original_data, na.action = stats::na.omit)
  
  # creating the fixed design matrix
  original_xo <- model.matrix(formula, original_stmodel_frame)
  
  # creating the response column
  original_yo <- model.response(original_stmodel_frame)
  
  # order the data by space within time
  spint <- storder(
    data = original_data,
    xcoord = xcoord,
    ycoord = ycoord,
    tcoord = tcoord,
    h_options = h_options
  )
  
  # create the model frame using the provided formula
  ordered_stmodel_frame <- model.frame(formula, spint$ordered_data_o,
                                       na.action = stats::na.omit
  )
  
  # creating the fixed design matrix
  ordered_xo <- model.matrix(formula, ordered_stmodel_frame)
  
  # creating the response column
  ordered_yo <- model.response(ordered_stmodel_frame)
  
  # creating the data object
  data_object <- list(
    formula = formula,
    original_data = original_data,
    original_xo = original_xo,
    original_yo = original_yo,
    ordered_data_dense = spint$ordered_data_dense,
    ordered_data_o = spint$ordered_data_o,
    hdist_seconds = spint$hdist_seconds,
    h_s_small = spint$h_s_small,
    h_t_small = spint$h_t_small,
    n_s = spint$n_s,
    n_t = spint$n_t,
    o_index = spint$o_index,
    m_index = spint$m_index,
    h_s_large = spint$h_s_large,
    h_t_large = spint$h_t_large,
    key_s = spint$key_s,
    key_t = spint$key_t,
    ordered_xo = ordered_xo,
    ordered_yo = ordered_yo,
    h_options = h_options
  )
  
  # returning the data object
  return(data_object)
}

#' Create a Distance Matrix
#'
#' @param coord1 Coordinate 1.
#' @param coord2 Coordinate 2.
#' @param distmetric Distance metric. Defaults to \code{"euclidean"}.
#'
#' @return A distance matrix
#' @noRd

make_h <- function(coord1, coord2 = NULL, distmetric = "euclidean") {
  
  # show the available distance metrics
  distmetric <- match.arg(distmetric)
  
  # calling the appropriate distance calculation
  switch(
    distmetric,
    euclidean = eucdist(coord1, coord2),
    stop("invalid distance metric")
  )
}



# compute the euclidean distance
eucdist <- function(coord1, coord2 = NULL) {
  if (is.null(coord2)) {
    # euclidean distance if 1d
    eucdist_1d <- sqrt(outer(coord1, coord1, sqr_dif))
    
    # return the distance
    return(eucdist_1d)
  } else {
    
    # euclidean distance if 2d
    eucdist_2d <- sqrt(outer(coord1, coord1, sqr_dif) + outer(coord2, coord2, sqr_dif))
    
    # return the distance
    return(eucdist_2d)
  }
}

# a squared difference function for outer to call
sqr_dif <- function(a, b) {
  (a - b)^2
}

#' Make the Inversion Object.
#'
#' @param covparam_object A covariance parameter object from \code{make_covparam_object()}.
#'
#' @param chol Should the Cholesky decomposition be used? If \code{FALSE},
#'   efficient inversion algorithms are implemented. Defaults to \code{FALSE}.
#'
#' @param co The covariance at prediction locations (if specified)
#'
#' @param condition A small number added to the diagonals of matrices before
#'   inverting them to prevent ill-conditioning (defaults to \code{1e-4}).
#'
#' @param h_s_large A spatial distance matrix of all spatio-temporal observations (if specified)
#'
#' @param h_t_large A temporal distance matrix of all spatio-temopral observations (if specified)
#'
#' @param h_s_small A spatial distance matrix of all spatial locations (if specified)
#'
#' @param h_t_small A temporal distance matrix of all temporal locations (if specified)
#'
#' @param logdet Should the log determinant be returned? (defaults to \code{FALSE}).
#'
#' @param m_index An index of missing values (from the space time cube).
#'
#' @param o_index An index of observed values (from the space time cube).
#'
#' @param s_cor The spatial correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'   }
#'
#' @param t_cor The temporal correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'     \item{\code{tent}}{The tent (linear with sill) correlation.}
#'   }
#'
#' @param xo The fixed effects design matrix.
#'
#' @param yo A response vector.
#'
#' @return A list of relevant inversion information.
#' @noRd

make_invert_object <- function(covparam_object,
                               chol,
                               co = NULL,
                               condition,
                               h_s_large = NULL,
                               h_t_large = NULL,
                               h_s_small = NULL,
                               h_t_small = NULL,
                               logdet,
                               m_index = NULL,
                               o_index = NULL,
                               s_cor,
                               t_cor,
                               xo,
                               yo) {
  
  # making the right multiplication matrix
  xyc_o <- cbind(xo, yo, co)
  
  # error message if small distance matrices are not provided and
  # cholesky is not being used
  if (!chol & (is.null(h_s_small) | is.null(h_t_small))) {
    stop("If not using Cholesky decomposition, h_s_small and h_t_small must be provided")
  }
  
  # error message if large distance matrices are not provided and
  # cholesky is being used
  if (chol & (is.null(h_s_large) | is.null(h_t_large))) {
    stop("If using Cholesky decomposition, h_s_large and h_t_large must be provided")
  }
  # storing the spatial dimension
  n_s <- nrow(h_s_small)
  
  # storing the temporal dimension
  n_t <- nrow(h_t_small)
  
  # creating the invert_object
  invert_object <- structure(list(
    covparams = covparam_object,
    chol = chol,
    condition = condition,
    logdet = logdet,
    h_s_small = h_s_small,
    h_t_small = h_t_small,
    h_s_large = h_s_large,
    h_t_large = h_t_large,
    o_index = o_index,
    m_index = m_index,
    n_s = n_s,
    n_t = n_t,
    s_cor = s_cor,
    t_cor = t_cor,
    xo = xo,
    yo = yo,
    co = co,
    xyc_o = cbind(xo, yo, co)
  ),
  class = class(covparam_object)
  )
  return(invert_object)
}

make_newdata_h <- function(newdata_coord1, data_coord1, newdata_coord2 = NULL, data_coord2 = NULL, distmetric = "euclidean") {
  
  # show the available distance metrics
  distmetric <- match.arg(distmetric)
  
  # calling the appropriate distance calculation
  switch(
    distmetric,
    euclidean = eucdist_newdata(newdata_coord1, data_coord1, newdata_coord2, data_coord2),
    stop("invalid distance metric")
  )
}

# compute the euclidean distance
eucdist_newdata <- function(newdata_coord1, data_coord1, newdata_coord2, data_coord2) {
  if (is.null(newdata_coord2) & is.null(data_coord2)) {
    
    # euclidean distance if 1d
    eucdist_1d <- sqrt(outer(newdata_coord1, data_coord1, sqr_dif))
    
    # return the distance
    return(eucdist_1d)
  } else {
    
    # euclidean distance if 2d
    eucdist_2d <- sqrt(outer(newdata_coord1, data_coord1, sqr_dif) +
                         outer(newdata_coord2, data_coord2, sqr_dif))
    
    # return the distance
    return(eucdist_2d)
  }
}

make_r <- function(h, range, structure = c("exponential", "spherical", "gaussian", "tent")) {
  
  # call the appropriate correlation function
  switch(
    structure,
    exponential = r_exp(h, range),
    spherical = r_sph(h, range),
    gaussian = r_gau(h, range),
    tent = r_tent(h, range),
    stop("Choose a valid covariance structure")
  )
}

# the correlation functions are parameterized in terms of their effective ranges

# exponential correlation
r_exp <- function(h, range) {
  
  # compute the exponential correlation
  r <- exp(-(3 * (h / range)))
  
  # return the exponential correlation
  return(r)
}

# spherical correlation
r_sph <- function(h, range) {
  
  # compute the spherical correlation
  r <- (h <= range) * (1 - (3 / 2) * (h / range) + (1 / 2) * (h / range)^3)
  
  # return the spherical correlation
  return(r)
}

# gaussian correlation
r_gau <- function(h, range) {
  
  # compute the gaussian correlation
  r <- exp(-(3 * (h / range)^2))
  
  # return the gaussian correlation
  return(r)
}

# tent correlation (only valid in 1d)
r_tent <- function(h, range) {
  
  # compute the tent correlation
  r <- (1 - h / range) * (h <= range)
  
  # return the tent correlation
  return(r)
}

make_sigma <- function(de, r_mx, ie, v_ie, e = 1, scale = FALSE) {
  
  # a warning message if the diagonal bug (being slightly greater than 1) appears
  if (any(r_mx > 1)) {
    stop("Bug - Diagonal of correlation matrix different from one")
  }
  
  # If the scale = FALSE, the standard variance parameters are used
  if (!scale) {
    
    # create the sigma matrix (equal to the weighted average of
    # dependent correlation and independent correlation)
    sigma <- de * r_mx + ie * (r_mx == 1)
    
    # return the sigma matrix
    return(sigma)
  }
  # If the scale = TRUE, the scaled variance parameters are used
  else if (scale) {
    
    # create the proportion of independent error as 1 - the proportion
    # of dependent error
    v_de <- 1 - v_ie
    
    # create the sigma matrix (equal to the weighted average of
    # dependent correlation and independent correlation)
    sigma <- e * v_de * r_mx + e * v_ie * (r_mx == 1)
    
    # return the sigma matrix
    return(sigma)
  } else {
    
    # a warning if scale is not appropriately specified
    stop("scale must be TRUE or FALSE")
  }
}


#' Make a covariance matrix
#'
#' @param covparam_object A covparam object
#'
#' @param h_s_large A spatial distance matrix of all spatio-temporal observations (if specified)
#'
#' @param h_t_large A temporal distance matrix of all spatio-temopral observations (if specified)
#'
#' @param s_cor The spatial correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'   }
#'
#' @param t_cor The temporal correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'     \item{\code{tent}}{The tent (linear with sill) correlation.}
#'   }
#'
#' @return A covariance matrix.
#' @noRd

make_stcovariance <- function(covparam_object,
                              h_s_large,
                              h_t_large,
                              s_cor,
                              t_cor) {
  
  # call the appropriate generic
  UseMethod("make_stcovariance", object = covparam_object)
}

#' @name make_stcovariance
#'
#' @method make_stcovariance productsum
#' @noRd

make_stcovariance.productsum <- function(covparam_object,
                                         h_s_large,
                                         h_t_large,
                                         s_cor,
                                         t_cor) {
  # creating the spatial correlation matrix
  r_s_large <- make_r(
    h = h_s_large,
    range = covparam_object[["s_range"]],
    structure = s_cor
  )
  
  # creating the temporal correlation matrix
  r_t_large <- make_r(
    h = h_t_large,
    range = covparam_object[["t_range"]],
    structure = t_cor
  )
  
  # creating the spatio-temporal correlation matrix
  r_st <- r_s_large * r_t_large
  
  # make the spatial sigma matrix
  sigma_s <- make_sigma(
    de = covparam_object[["s_de"]],
    r_mx = r_s_large,
    ie = covparam_object[["s_ie"]]
  )
  
  # make the temporal sigma matrix
  sigma_t <- make_sigma(
    de = covparam_object[["t_de"]],
    r_mx = r_t_large,
    ie = covparam_object[["t_ie"]]
  )
  
  # make the spatio-temporal sigma matrix
  sigma_st <- make_sigma(
    de = covparam_object[["st_de"]],
    r_mx = r_st,
    ie = covparam_object[["st_ie"]]
  )
  
  # making the overall covariance matrix
  sigma <- sigma_s + sigma_t + sigma_st
  
  # returning the overall covariance matrix
  return(sigma)
}

#' @name make_stcovariance
#'
#' @method make_stcovariance sum_with_error
#' @noRd

make_stcovariance.sum_with_error <- function(covparam_object,
                                             h_s_large,
                                             h_t_large,
                                             s_cor,
                                             t_cor) {
  # creating the spatial correlation matrix
  r_s_large <- make_r(
    h = h_s_large,
    range = covparam_object[["s_range"]],
    structure = s_cor
  )
  
  # creating the temporal correlation matrix
  r_t_large <- make_r(
    h = h_t_large,
    range = covparam_object[["t_range"]],
    structure = t_cor
  )
  
  # creating the spatio-temporal correlation matrix
  r_st <- r_s_large * r_t_large
  
  # make the spatial sigma matrix
  sigma_s <- make_sigma(
    de = covparam_object[["s_de"]],
    r_mx = r_s_large,
    ie = covparam_object[["s_ie"]]
  )
  
  # make the temporal sigma matrix
  sigma_t <- make_sigma(
    de = covparam_object[["t_de"]],
    r_mx = r_t_large,
    ie = covparam_object[["t_ie"]]
  )
  
  # making the spatio-temporal matrix
  sigma_st <- make_sigma(
    de = 0,
    r_mx = r_st,
    ie = covparam_object[["st_ie"]]
  )
  
  # make the overall covariance matrix
  sigma <- sigma_s + sigma_t + sigma_st
  
  # return the overall covariance matrix
  return(sigma)
}


#' @name make_stcovariance
#'
#' @method make_stcovariance product
#' @noRd

make_stcovariance.product <- function(covparam_object,
                                      h_s_large,
                                      h_t_large,
                                      s_cor,
                                      t_cor) {
  # make the spatial scaled correlation matrix
  # (e = 1 because the variance must be 1 as these are
  # correlation matrices)
  r_s <- make_sigma(
    r_mx = make_r(
      h = h_s_large,
      range = covparam_object[["s_range"]],
      structure = s_cor
    ),
    v_ie = covparam_object[["v_s"]],
    e = 1,
    scale = TRUE
  )
  
  # make the temporal scaled correlation matrix
  # (e = 1 because the variance must be 1 as these are
  # correlation matrices)
  r_t <- make_sigma(
    r_mx = make_r(
      h = h_t_large,
      range = covparam_object[["t_range"]],
      structure = t_cor
    ),
    v_ie = covparam_object[["v_t"]],
    e = 1,
    scale = TRUE
  )
  
  # make the spatio-temporal scaled correlation matrix
  r_st <- r_s * r_t
  
  # make overall covariance matrix (the independent
  # error must be zero for the product covariance)
  sigma_st <- make_sigma(
    de = covparam_object[["st_de"]],
    r_mx = r_st,
    ie = 0
  )
  
  # making it explicit that sigma is the same as sigma_st
  sigma <- sigma_st
  
  # returning the overall covariance
  return(sigma)
}

#' Make a Semivariogram Matrix
#'
#' @inheritParams make_stcovariance
#'
#' @return A semivariogram matrix
#' @noRd

make_stsemivariogram <- function(covparam_object,
                                 h_s_large,
                                 h_t_large,
                                 s_cor,
                                 t_cor) {
  
  # call the appropriate generic
  UseMethod("make_stsemivariogram", object = covparam_object)
}

#' @name make_stsemivariogram
#'
#' @method make_stsemivariogram productsum
#' @noRd

make_stsemivariogram.productsum <- function(covparam_object,
                                            h_s_large,
                                            h_t_large,
                                            s_cor,
                                            t_cor) {
  
  # make the productsum semivariogram
  # taking the variance parameters from the covparam_object
  variances <- c(covparam_object[c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie")])
  
  # computing the semivariogram
  gamma <- sum(variances) -
    make_stcovariance.productsum(
      covparam_object = covparam_object,
      h_s_large = h_s_large,
      h_t_large = h_t_large,
      s_cor = s_cor,
      t_cor = t_cor
    )
  
  # returning the semivariogram
  return(gamma)
}

#' @name make_stsemivariogram
#'
#' @method make_stsemivariogram sum_with_error
#' @noRd

make_stsemivariogram.sum_with_error <- function(covparam_object,
                                                h_s_large,
                                                h_t_large,
                                                s_cor,
                                                t_cor) {
  
  # taking the variance parameters from the covparam_object
  variances <- c(covparam_object[c("s_de", "s_ie", "t_de", "t_ie", "st_ie")])
  
  # computing the semivariogram
  gamma <- sum(variances) -
    make_stcovariance.sum_with_error(
      covparam_object = covparam_object,
      h_s_large = h_s_large,
      h_t_large = h_t_large,
      s_cor = s_cor,
      t_cor = t_cor
    )
  
  # returning the semivariogram
  return(gamma)
}

#' @name make_stsemivariogram
#'
#' @method make_stsemivariogram product
#' @noRd

make_stsemivariogram.product <- function(covparam_object, h_s_large, h_t_large,
                                         s_cor, t_cor) {
  
  # taking the variance parameters from the covparam_object
  variances <- c(covparam_object[c("st_de")])
  
  # computing the semivariogram
  gamma <- sum(variances) -
    make_stcovariance.product(
      covparam_object = covparam_object,
      h_s_large = h_s_large, h_t_large = h_t_large,
      s_cor = s_cor, t_cor = t_cor
    )
  
  # returning the semivariogram
  return(gamma)
}

make_z_factor <- function(factor_index) {
  
  # if factor_index is not a character vector, store it as one
  if (!is.character(factor_index)) {
    factor_index <- as.character(factor_index)
  }
  
  # make the fixed effect spatial design matrix
  z <- stats::model.matrix(~ factor_index - 1)
  
  # return the matrix
  return(z)
}

minus2loglik <- function(invert_object, invert_output) {
  
  # appropriate generic for likelihood based estimation
  UseMethod("minus2loglik", object = invert_object)
}

# compute the restricted log likelihood
minus2loglik.reml <- function(invert_object, invert_output) {
  
  # storing the sample size
  n <- nrow(invert_output$sigmainv_o)
  
  # computing the current beta estimate
  betaest_output <- betaest(
    xo = invert_object$xo,
    sigmainv_xyo = invert_output$sigmainv_o,
    condition = invert_object$condition,
    return_estlist = TRUE
  )
  
  # computing sigma inverse times the residual vector
  siginv_r <- betaest_output$estlist$sigmainv_yo -
    betaest_output$estlist$sigmainv_xo %*% betaest_output$betahat
  
  # computing the quadratic residual function
  r_siginv_r <- t(invert_object$yo - invert_object$xo %*% betaest_output$betahat) %*%
    siginv_r
  
  # storing n - p
  nminusp <- n - betaest_output$estlist$p
  
  # computing the restricted likelihood
  m2ll <- invert_output$logdet +
    (nminusp) * (log(r_siginv_r) + 1 + log(2 * pi / nminusp)) +
    betaest_output$estlist$ldet_cicb
  
  # returning the restricted likelihood
  return(m2ll)
}

multiply_z <- function(mx, z_type, n_s, n_t, side = c("right", "left", "p_right", "p_left", "pz_z", "z_pz")) {
  
  # show the appropriate multiplication side arguments
  side <- match.arg(side)
  
  # performing the appropriate multiplication
  switch(side,
         right = multiply_z_r(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         left = multiply_z_l(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         p_right = multiply_zp_r(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         p_left = multiply_zp_l(mx = mx, z_type = z_type, n_s = n_s, n_t = n_t),
         pz_z = multiply_zp_z(z_type = z_type, n_s = n_s, n_t = n_t),
         z_pz = multiply_z_zp(z_type = z_type, n_s = n_s, n_t = n_t),
         stop("Please choose an approporiate multiplication function for Z")
  )
}

# multiplying by z on the right (AZ)
multiply_z_r <- function(mx, z_type, n_s, n_t) {
  
  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)) {
    mx <- matrix(mx, ncol = 1)
  }
  
  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t
  
  # find the row dimension of the matrix A
  n_row <- nrow(mx)
  
  # multiplying by Z_s
  if (z_type == "spatial") {
    
    # iterating through the appropriate row sums
    a_zs <- vapply(
      1:n_s,
      function(a) rowSums(mx[, seq(a, n_st, by = n_s)]),
      double(n_row)
    )
    
    # returning the appropriate product
    return(a_zs)
  } else if (z_type == "temporal") { # multiplying by Z_t
    
    # iterating through the appropriate row sums
    a_zt <- vapply(
      1:n_t,
      function(a) rowSums(mx[, seq(n_s * (a - 1) + 1, n_s * a, by = 1), drop = FALSE]),
      double(n_row)
    )
    
    # returning the appropriate product
    return(a_zt)
  } else {
    # returning an error message if no proper z product is chosen
    stop("inappropriate z type")
  }
}

# multiplying by z transpose on the left (Z'A)
# if A is symmetric this equals AZ
multiply_zp_l <- function(mx, z_type, n_s, n_t) {
  
  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)) {
    mx <- matrix(mx, ncol = 1)
  }
  
  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t
  
  # find the column dimension of the matrix A
  n_col <- ncol(mx)
  
  # multiplying by Z_s
  if (z_type == "spatial") {
    
    # iterating through the appropriate column sums
    zps_a <- t(vapply(
      1:n_s,
      function(a) colSums(mx[seq(a, n_st, by = n_s), , drop = FALSE]),
      double(n_col)
    ))
    
    # returning the appropriate product
    return(zps_a)
  } else if (z_type == "temporal") { #   # multiplying by Z_s
    
    # iterating through the appropriate column sums
    zpt_a <- t(vapply(
      1:n_t,
      function(a) colSums(mx[seq(n_s * (a - 1) + 1, n_s * a, by = 1), , drop = FALSE]),
      double(n_col)
    ))
    
    # returning the appropriate product
    return(zpt_a)
  } else {
    # returning an error message if no proper z product is chosen
    stop("innapropriate Z type")
  }
}

# multiplying by z on the left (ZA)
multiply_z_l <- function(mx, z_type, n_s, n_t) {
  
  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)) {
    mx <- matrix(mx, ncol = 1)
  }
  
  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t
  
  # multiplying by Z_s
  if (z_type == "spatial") {
    
    # computing the appropriate sums
    zs_a <- mx[rep(seq(1, n_s), times = n_t), , drop = FALSE]
    
    # returning the appropriate product
    return(zs_a)
  } else if (z_type == "temporal") {
    
    # computing the appropriate sums
    zt_a <- mx[rep(seq(1, n_t), each = n_s), , drop = FALSE]
    
    # returning the appropriate product
    return(zt_a)
  } else {
    stop("innapropriate Z type")
  }
}

# multiplying by z transpose on the right (AZ')
# if A is symmetric this equals ZA
multiply_zp_r <- function(mx, z_type, n_s, n_t) {
  
  # is mx is a 1d vector, turn it into a matrix
  if (is.vector(mx)) {
    mx <- matrix(mx, ncol = 1)
  }
  
  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t
  
  # multiplying by Z_s
  if (z_type == "spatial") {
    
    # computing the appropriate sums
    a_zps <- mx[, rep(seq(1, n_s), times = n_t), drop = FALSE]
    
    # returning the appropriate product
    return(a_zps)
  } else if (z_type == "temporal") { # multiplying by Z_t
    
    # computing the appropriate sums
    a_zpt <- mx[, rep(seq(1, n_t), each = n_s), drop = FALSE]
    
    # returning the appropriate product
    return(a_zpt)
  } else {
    stop("innapropriate Z type")
  }
}

# multiplying z transpose and z (of the same type)
multiply_zp_z <- function(z_type, n_s, n_t) {
  
  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t
  
  # multiplying Z_s
  if (z_type == "spatial") {
    
    # computing the appropriate sums
    zps_zs <- diag(n_t, nrow = n_s, ncol = n_s)
    
    # returning the appropriate product
    return(zps_zs)
  } else if (z_type == "temporal") { # multiplying Z_t
    
    # computing the appropriate sums
    zpt_zt <- diag(n_s, nrow = n_t, ncol = n_t)
    
    # returning the appropriate product
    return(zpt_zt)
  } else {
    stop("inappropriate z type")
  }
}

# multiplying z and z transpose (of the same type)
multiply_z_zp <- function(z_type, n_s, n_t) {
  
  # compute the sample size of a dense data rectangle
  n_st <- n_s * n_t
  
  # multiplying by Z_s
  if (z_type == "spatial") {
    
    # computing the appropriate sums
    zs_zps <- kronecker(matrix(1, nrow = n_t, ncol = n_t), diag(1, nrow = n_s, ncol = n_s))
    
    # returning the appropriate product
    return(zs_zps)
  } else if (z_type == "temporal") {
    
    # computing the appropriate sums
    zt_zpt <- kronecker(diag(1, nrow = n_t, ncol = n_t), matrix(1, nrow = n_s, ncol = n_s))
    
    # returning the appropriate product
    return(zt_zpt)
  } else {
    stop("inappropriate z type")
  }
}

#' Predict
#'
#' Compute best linear unbiased predictions (Kriging).
#'
#' @param object A model object of class \code{stlmm}.
#'
#' @param newdata A data frame containing columns whose names match the names
#'   of the x-coordinate, y-coordinate, t-coordinate, and predictor variables
#'   in \code{object}.
#'
#' @param interval The interval type. \code{"none"} implies point estimates,
#'   \code{"confidence"} implies point estimates whose standard errors are related
#'   to the mean. \code{"predction"} implies point estimates whose standard errors
#'   are related to a prediction.
#'
#' @param se.fit Should the standard error of the point estimate be returned?
#'   Defaults to \code{TRUE}.
#'
#' @param predcov Should the appropriate full covariance matrix of predictions
#' be returned? Deafults to \code{FALSE}.
#'
#' @param ... Additional arguments
#'
#' @return A list containing the point estimates, standard errors (if requested), and
#'   prediction covariance matrix (if requested).
#' @noRd

predict.stlmm <- function(object,
                          newdata,
                          interval = c("none", "confidence", "prediction"),
                          se.fit = TRUE,
                          predcov = FALSE,
                          ...) {
  
  # show the interval types in documentation
  interval <- match.arg(interval)
  
  # calling the appropriate interval
  pred_output <- switch(
    interval,
    "none" = predict.stlmm_none(
      object = object,
      newdata = newdata,
      ...
    ),
    "confidence" = predict.stlmm_confidence(
      object = object,
      newdata = newdata,
      se.fit = se.fit,
      predcov = predcov,
      ...
    ),
    "prediction" = predict.stlmm_prediction(
      object = object,
      newdata = newdata,
      se.fit = se.fit,
      predcov = predcov,
      ...
    ),
    stop("Must choose confidence or prediction interval")
  )
  
  # returning the interval argument in the prediction output
  pred_output$interval <- interval
  
  # returning the prediction output
  return(pred_output)
}


predict.stlmm_none <- function(object,
                               newdata,
                               ...) {
  
  # creating the model frame
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
                                       na.action = stats::na.omit
  )
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)
  
  # computing the estimate
  fit <- newdata_xo %*% object$Coefficients
  
  # storing the output as a list
  pred_output <- list(fit = fit)
  
  # returning the output
  return(pred_output)
}

predict.stlmm_confidence <- function(object,
                                     newdata,
                                     se.fit,
                                     predcov,
                                     ...) {
  # creating the model frame
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
                                       na.action = stats::na.omit
  )
  
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)
  
  # computing the estimate
  fit <- newdata_xo %*% object$Coefficients
  
  
  if (predcov) {
    # computing the full covariance matrix of the estimates
    predcov <- newdata_xo %*% object$CovCoefficients %*% t(newdata_xo)
    if (se.fit) {
      # computing the standard errors
      se.fit <- sqrt(diag(predcov))
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  } else {
    # setting the full covariance matrix NULL if requested
    predcov <- NULL
    if (se.fit) {
      # the number of predictions
      npred <- nrow(newdata_xo)
      
      # computing the standard errors
      var.fit <- vapply(
        1:npred,
        function(x) newdata_xo[x, , drop = FALSE] %*% object$CovCoefficients %*% t(newdata_xo[x, , drop = FALSE]),
        double(1)
      )
      se.fit <- sqrt(var.fit)
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  }
  
  # storing output in a list
  pred_output <- list(fit = fit, se.fit = se.fit, predcov = predcov)
  
  # removing the non-NULL elements
  pred_output <- pred_output[!vapply(pred_output, function(x) is.null(x), logical(1))]
  
  # returnign the output
  return(pred_output)
}

predict.stlmm_prediction <- function(object,
                                     newdata,
                                     se.fit,
                                     predcov,
                                     ...) {
  
  # creating the model frame
  newdata_stmodel_frame <- model.frame(object$formula[-2], newdata,
                                       na.action = stats::na.omit
  )
  # creating the fixed design matrix
  # -2 is to remove response from formula
  newdata_xo <- model.matrix(object$formula[-2], newdata_stmodel_frame)
  
  # make spatial distance matrices in 1d
  if (is.null(object$coordnames$ycoord)) {
    
    # distances between new observations and data
    newdata_data_h_s_large <- make_newdata_h(
      newdata_coord1 = newdata[[object$coordnames$xcoord]],
      data_coord1 = object$coords$xcoord,
      distmetric = object$h_options$h_s_distmetric
    )
    
    # distances between new obsevations and new observations
    newdata_h_s_large <- make_h(
      coord1 = newdata[[object$coordnames$xcoord]],
      distmetric = object$h_options$h_s_distmetric
    )
  } else { # make spatial distance matrices in 2d
    
    # distances between new observations and data
    newdata_data_h_s_large <- make_newdata_h(
      newdata_coord1 = newdata[[object$coordnames$xcoord]],
      data_coord1 = object$coords$xcoord,
      newdata_coord2 = newdata[[object$coordnames$ycoord]],
      data_coord2 = object$coords$ycoord,
      distmetric = object$h_options$h_s_distmetric
    )
    # distances between new observations and new observations
    newdata_h_s_large <- make_h(
      coord1 = newdata[[object$coordnames$xcoord]],
      coord2 = newdata[[object$coordnames$ycoord]],
      distmetric = object$h_options$h_s_distmetric
    )
  }
  
  # make temporal distance matrices
  # distances between new observations and data
  newdata_data_h_t_large <- make_newdata_h(
    newdata_coord1 = newdata[[object$coordnames$tcoord]],
    data_coord1 = object$coords$tcoord,
    distmetric = object$h_options$h_t_distmetric
  )
  
  # distances between new observations and new observations
  newdata_h_t_large <- make_h(
    coord1 = newdata[[object$coordnames$tcoord]],
    distmetric = object$h_options$h_t_distmetric
  )
  
  
  # make covariance matrix between new observations and data
  newdata_data_stcovariance <- make_stcovariance(
    covparam_object = object$CovarianceParameters,
    h_s_large = newdata_data_h_s_large,
    h_t_large = newdata_data_h_t_large,
    s_cor = object$CovarianceForms[["s_cor"]],
    t_cor = object$CovarianceForms[["t_cor"]]
  )
  
  
  # make the object to invert and multiply by covariance on the right
  invert_object <- make_invert_object(
    covparam_object = object$CovarianceParameters,
    chol = object$chol,
    condition = object$condition,
    co = t(newdata_data_stcovariance),
    h_s_small = object$data_object$h_s_small,
    h_t_small = object$data_object$h_t_small,
    h_s_large = object$data_object$h_s_large,
    h_t_large = object$data_object$h_t_large,
    logdet = FALSE,
    m_index = object$data_object$m_index,
    o_index = object$data_object$o_index,
    s_cor = object$CovarianceForms[["s_cor"]],
    t_cor = object$CovarianceForms[["t_cor"]],
    xo = NULL,
    yo = NULL
  )
  
  # compute the inverse
  invert_output <- invert(invert_object)
  
  # store sigmainverse %*% X
  invxo <- object$invert_output$sigmainv_o[, 1:(ncol(object$invert_output$sigmainv_o) - 1), drop = FALSE]
  
  # store cov(new, data) %*% sigmainverse %*% X
  newdata_invxo <- newdata_data_stcovariance %*% invxo
  
  # computing the BLUP estimates
  fit <- newdata_xo %*% object$Coefficients + # beta hat estimates
    newdata_data_stcovariance %*% object$invert_output$sigmainv_o[, ncol(object$invert_output$sigmainv_o), drop = FALSE] -
    newdata_invxo %*% object$Coefficients # blup residuals
  
  
  
  if (predcov) {
    # computing the covariance of new and new
    newdata_stcovariance <- make_stcovariance(
      covparam_object = object$CovarianceParameters,
      h_s_large = newdata_h_s_large,
      h_t_large = newdata_h_t_large,
      s_cor = object$CovarianceForms[["s_cor"]],
      t_cor = object$CovarianceForms[["t_cor"]]
    )
    
    # storing Xo - cov(new, data) %*% sigmainverse %*% X
    H <- newdata_xo - newdata_invxo
    
    # computing the prediction covariance matrix
    predcov <- newdata_stcovariance - newdata_data_stcovariance %*% invert_output$sigmainv_o +
      H %*% object$CovCoefficients %*% t(H)
    
    if (se.fit) {
      # computing the standard errors
      se.fit <- sqrt(diag(predcov))
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  } else {
    # setting the predcov matrix to NULL if not requested
    predcov <- NULL
    if (se.fit) {
      # storing the parameter names
      vparm_names <- c(
        "s_de", "s_ie", "t_de", "t_ie",
        "st_de", "st_ie"
      )
      
      # computing the overall variance by taking only the variance parameters in the
      # covariance parameter object
      varsum <- sum(object$CovarianceParameters[names(object$CovarianceParameters) %in% vparm_names])
      
      # computing the standard error
      se.fit <- vapply(
        1:nrow(newdata_xo),
        function(x) {
          H <- newdata_xo[x, , drop = FALSE] - newdata_invxo[x, , drop = FALSE]
          predvar <- varsum - newdata_data_stcovariance[x, , drop = FALSE] %*% invert_output$sigmainv_o[, x, drop = FALSE] +
            H %*% object$CovCoefficients %*% t(H)
          se.fit <- sqrt(predvar)
        },
        double(1)
      )
    } else {
      # setting the standard errors to NULL if not requested
      se.fit <- NULL
    }
  }
  
  # storing output in a list
  pred_output <- list(fit = fit, se.fit = se.fit, predcov = predcov)
  
  # returnign the non-NULL elements
  pred_output <- pred_output[!vapply(pred_output, function(x) is.null(x), logical(1))]
  
  # returning the output
  return(pred_output)
}

residuals.stlmm <- function(object, type = "raw", ...) {
  
  # calling the function for the raw residuals
  if (type == "raw") {
    
    # computing the raw residuals
    residuals <- object$model$Response - object$model$FixedDesignMatrix %*% object$Coefficients
  }
  
  # returning the type of the residuals as an attribute
  attr(residuals, "type") <- type
  
  # returning the residuals
  return(residuals)
}

#' Compute the Empirical Spatio-Temporal Semivariogram
#'
#' @param response A vector of response variables (not required if
#'   \code{h_response} is provided).
#'
#' @param xcoord A vector of x-coordinates (not required if
#'   \code{h_s_large} is provided).
#'
#' @param ycoord A vector of y-coordinates (not required if
#'   \code{h_s_large} is provided).
#'
#' @param tcoord A vector of t-coordinates (not required if
#'   \code{h_t_large} is provided).
#'
#' @param h_options A list containing options to compute distances if
#'   \code{response}, \code{xcoord}, \code{ycoord}, and \code{tcoord} are
#'   provided. Named arguments are
#'   \describe{
#'     \item{\code{h_t_distmetric}}{The temporal distance matrix (defaults to
#'     \code{"euclidean"}).}.
#'     \item{\code{h_s_distmetric}}{The spatial distance matrix (defaults to
#'     \code{"euclidean"}).}
#'  }
#'
#' @param h_response A distance matrix of paired differences of the response
#'   variable.
#'
#' @param h_s_large A distance matrix of paired differences of the spatial
#'   locations.
#'
#' @param h_t_large A distance matrix of paired differences of the temporal
#'   locations.
#'
#' @param stempsv_options A list containing additional options. Named arguments
#'   are
#'   \describe{
#'     \item{\code{n_s_lag}}{The number of spatial distance classes (defaults to 16).}
#'     \item{\code{n_t_lag}}{The number of temporal distance classes (defaults to 16).}
#'     \item{\code{h_s_max}}{The maximum spatial distance. Deafaults to half the
#'       maximum distance in the spatial domain.}
#'     \item{\code{h_t_max}}{The maximum temporal distance. Deafaults to half the
#'       maximum distance in the temporal domain.}
#'   }
#'
#' @return A data frame whose columns are
#'   \describe{
#'     \item{\code{gammahat}}{The estimated empirical semivariogram value for
#'       for the spatio-temporal distance class.}
#'     \item{\code{n}}{The number of pairs for the spatio-temporal distance class}
#'     \item{\code{h_s_avg}}{The average spatial distance in the spatio-temporal
#'       distance class.}
#'     \item{\code{h_t_avg}}{The average temporal distance in the spatio-temporal
#'       distance class.}
#'   }
#' @noRd

stempsv <- function(response,
                    xcoord,
                    ycoord = NULL,
                    tcoord,
                    h_options = NULL,
                    h_response = NULL,
                    h_s_large = NULL,
                    h_t_large = NULL,
                    stempsv_options = NULL) {
  
  # setting default options if none are given
  if (is.null(stempsv_options)) {
    stempsv_options <- list(n_s_lag = 16, n_t_lag = 16, h_s_max = NULL, h_t_max = NULL)
  }
  
  # setting default h options if none are given
  if (is.null(h_options)) {
    h_options <- list(h_t_distmetric = "euclidean", h_s_distmetric = "euclidean")
  }
  
  # creating a large spatial distance matrix if not provided
  if (is.null(h_s_large)) {
    h_s_large <- make_h(coord1 = xcoord, coord2 = ycoord, distmetric = h_options$h_s_distmetric)
  }
  
  # creating a large temporal distance matrix if not provided
  if (is.null(h_t_large)) {
    h_t_large <- make_h(coord1 = tcoord, distmetric = h_options$h_t_distmetric)
  }
  
  # creating a large squared response matrix if not provided
  if (is.null(h_response)) {
    h_response <- make_h(coord1 = response, distmetric = "euclidean")^2
  }
  
  # only storing the upper triangular portions of the spatial, temporal,
  # and squared distance matrices
  h_s_large <- h_s_large[upper.tri(h_s_large, diag = F)]
  h_t_large <- h_t_large[upper.tri(h_t_large, diag = F)]
  h_response <- h_response[upper.tri(h_response, diag = F)]
  
  # setting a max for the spatial disance if none is provided
  # (as the max distance over 2)
  if (is.null(stempsv_options$h_s_max)) {
    stempsv_options$h_s_max <- max(h_s_large) / 2
  }
  
  # setting a max for the temporal disance if none is provided
  # (as the max distance over 2)
  if (is.null(stempsv_options$h_t_max)) {
    stempsv_options$h_t_max <- max(h_t_large) / 2
  }
  
  hmax_index <- (h_s_large <= stempsv_options$h_s_max) & (h_t_large <= stempsv_options$h_t_max)
  h_s_large <- h_s_large[hmax_index]
  h_t_large <- h_t_large[hmax_index]
  h_response <- h_response[hmax_index]
  
  s_lags <- c(-0.1, seq(0, stempsv_options$h_s_max, length.out = stempsv_options$n_s_lag))
  h_s_index <- cut(h_s_large, s_lags, include.lowest = T)
  t_lags <- c(-0.1, seq(0, stempsv_options$h_t_max, length.out = stempsv_options$n_t_lag))
  h_t_index <- cut(h_t_large, t_lags, include.lowest = T)
  st_index <- interaction(h_s_index, h_t_index)
  
  
  h_s_avg <- tapply(h_s_large, st_index, mean)
  h_t_avg <- tapply(h_t_large, st_index, mean)
  gammahat <- tapply(h_response, st_index, FUN = function(x) mean(x) / 2)
  n <- tapply(h_response, st_index, length)
  
  return_output <- na.omit(data.frame(gammahat, n, h_s_avg, h_t_avg))
  attr(return_output, "na.action") <- NULL
  row.names(return_output) <- NULL
  return(return_output)
}

storder <- function(data, xcoord, ycoord = NULL, tcoord, h_options) {
  
  # find unique temporal coordinates
  key_t <- unique(data[, tcoord, drop = FALSE])
  
  # order the unique temporal coodrinates
  key_t <- key_t[order(key_t[[tcoord]]), , drop = FALSE]
  
  # compute the ordered small temporal distance matrix
  h_t_small <- make_h(
    coord1 = key_t[[tcoord]],
    distmetric = h_options$h_t_distmetric
  )
  
  # record the number of unique temporal observations
  n_t <- nrow(key_t)
  
  # create an index for the unique temporal locations
  key_t$tindex <- seq.int(1, n_t)
  
  # find unique spatial coordinates
  key_s <- unique(data[, c(xcoord, ycoord), drop = FALSE])
  
  # compute the small distance matrix in 1d
  if (is.null(ycoord)) {
    
    # order the unique spatial coodrinates
    key_s <- key_s[order(key_s[[xcoord]]), , drop = FALSE]
    
    # compute the ordered small spatial distance matrix
    h_s_small <- make_h(
      coord1 = key_s[[xcoord]],
      distmetric = h_options$h_s_distmetric
    )
  } else { # compute the small distance matrix in 1d
    
    # order the unique spatial coodrinates
    key_s <- key_s[order(key_s[[ycoord]], key_s[[xcoord]]), , drop = FALSE]
    
    # compute the ordered small spatial distance matrix
    h_s_small <- make_h(
      coord1 = key_s[[xcoord]],
      coord2 = key_s[[ycoord]],
      distmetric = h_options$h_s_distmetric
    )
  }
  
  # record the number of unique spatial observations
  n_s <- nrow(key_s)
  
  # create an index for the unique temporal locations
  key_s$sindex <- seq.int(1, n_s)
  
  # merge the data and the spatial data key
  # and then merge that data and the temporal data key
  data <- merge(merge(data, key_s), key_t)
  
  # create a grid containing every combination of spatial and temporal indices
  full_grid <- expand.grid(sindex = key_s$sindex, tindex = key_t$tindex)
  
  # create an overall index
  full_grid$index <- seq.int(1, n_t * n_s)
  
  # merge the grid and data
  data <- merge(full_grid, data, all = TRUE)
  
  # order the data by the ordered index (by space within time)
  data <- data[order(data$index), , drop = FALSE]
  
  # create a new variable, observed, which is logical indicating whether
  # the index value was observed in the original data
  data$observed <- !(is.na(data[[tcoord]]) & is.na(data[[xcoord]]))
  
  # find indicies of the observed data and "missing" data - which
  # are index values having no observation
  o_index <- data$index[data$observed]
  m_index <- data$index[!data$observed]
  
  # create a subsetted data frame of only the observed values
  ordered_data_o <- data[o_index, , drop = FALSE]
  
  # compute large distance matrices if needed
  if (h_options$h_large) {
    
    # compute the large distance matrices in 1d
    if (is.null(ycoord)) {
      
      # save the times
      hdist_start <- Sys.time()
      # compute the large spatial distance matrix
      h_s_large <- make_h(
        coord1 = ordered_data_o[[xcoord]],
        distmetric = h_options$h_s_distmetric
      )
      
      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = ordered_data_o[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
      hdist_end <- Sys.time()
      hdist_seconds <- as.numeric(hdist_end - hdist_start, units = "secs")
    } else { # compute the large distance matrices in 2d
      
      # compute the large spatial distance matrix
      # save the times
      hdist_start <- Sys.time()
      h_s_large <- make_h(
        coord1 = ordered_data_o[[xcoord]],
        coord2 = ordered_data_o[[ycoord]],
        distmetric = h_options$h_s_distmetric
      )
      
      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = ordered_data_o[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
      hdist_end <- Sys.time()
      hdist_seconds <- as.numeric(hdist_end - hdist_start, units = "secs")
    }
  } else { # set the large distance matrices equal to NULL if not requested
    h_s_large <- NULL
    h_t_large <- NULL
  }
  
  # return the relevant output
  return(list(
    ordered_data_dense = data,
    ordered_data_o = ordered_data_o,
    h_s_small = h_s_small,
    h_t_small = h_t_small,
    n_s = n_s,
    n_t = n_t,
    o_index = o_index,
    m_index = m_index,
    hdist_seconds = hdist_seconds,
    h_s_large = h_s_large,
    h_t_large = h_t_large,
    key_s = key_s,
    key_t = key_t
  ))
}

#' Simulate a Spatio-Temporal Random Variable
#'
#' @param object An covariance matrix or a \code{"covparam"} object.
#'
#' @param mu A mean vector with length equal to the number of rows in \code{data}
#'   or a scalar.
#'
#' @param size The number of independent simulations
#'
#' @param condition A small number added to the diagonals of matrices before
#'   inverting them to prevent ill-conditioning (defaults to \code{1e-4}).
#'
#' @param error The random error type
#' \describe{
#'   \item{\code{normal}}{All random effects are Gaussian and mutually independent.}
#'   \item{\code{component_squared}}{Gaussian, mutually independent random effects
#'     are simulated, squared, and rescaled to match the sample variance on the
#'     Gaussian scale.}
#'   \item{\code{sum_squared}}{Gaussian, mutually independent random effects
#'     are simulated, summed, squared, and rescaled to match the sample variance on the
#'     Gaussian scale.}
#' }
#'
#' @param xcoord A character vector specifying the column name of the x-coordinate
#'   variable in \code{data}.
#'
#' @param ycoord A character vector specifying the column name of the y-coordinate
#'   variable in \code{data}.
#'
#' @param tcoord A character vector specifying the column name of the t-coordinate (time)
#'   variable in \code{data}.
#'
#' @param data A data object containing all necessary variables.
#'
#' @param s_cor The spatial correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'   }
#'
#' @param t_cor The temporal correlation
#'   \describe{
#'     \item{\code{exponential}}{The exponential correlation (the default).}
#'     \item{\code{spherical}}{The spherical correlation.}
#'     \item{\code{gaussian}}{The Gaussian correlation.}
#'     \item{\code{tent}}{The tent (linear with sill) correlation.}
#'   }
#'
#' @param chol Should the Cholesky decomposition be used? If \code{FALSE},
#'   efficient inversion algorithms are implemented. Defaults to \code{FALSE}.
#'
#' @param h_options A list containing options to compute distances if
#'   \code{response}, \code{xcoord}, \code{ycoord}, and \code{tcoord} are
#'   provided. Named arguments are
#'   \describe{
#'     \item{\code{h_t_distmetric}}{The temporal distance matrix (defaults to
#'     \code{"euclidean"}).}
#'     \item{\code{h_s_distmetric}}{The spatial distance matrix (defaults to
#'     \code{"euclidean"}).}
#'  }
#'
#' @param ... Additional arguments.
#'
#'
#' @return A vector of random variables (if \code{size = 1}) or a matrix of
#'   random variables (if \code{size > 1}) whose columns indicate seprate
#'   simulations. The row order corresponds to the rows of the covariance
#'   matrix (if \code{object} is a matrix) or the rows of \code{data} (if
#'   \code{object} is a \code{covparam} object.
#' @noRd

strnorm <- function(object, mu, size, condition = 1e-4, error, ...) {
  
  # dispatching the appropriate generic
  UseMethod("strnorm", object = object)
}

#' @name strnorm
#' @method strnorm matrix
#' @noRd

strnorm.matrix <- function(object, mu, size, condition = 1e-4, error, ...) {
  # simulating a random variable if a covariance matrix is provided
  # record the sample size
  n_st <- nrow(object)
  
  # return an error if the length of mu does not equal the sample size
  if (length(mu) != 1 & length(mu) != n_st) {
    stop("choose mu as a vector having size nst or a single scalar")
  }
  
  # taking the lower triangular cholesky decomposition
  # (chol() returns the upper triangular)
  chol_siginv <- t(chol(object))
  
  # returning the simulated random vector as many times as "size" indicates
  return(
    vapply(
      1:size,
      FUN = function(x) mu + chol_siginv %*% rnorm(n_st),
      double(n_st)
    )
  )
}


#' @name strnorm
#'
#' @method strnorm default
#' @noRd

strnorm.default <- function(object,
                            mu,
                            size,
                            condition = 1e-4,
                            error,
                            xcoord,
                            ycoord = NULL,
                            tcoord,
                            data,
                            s_cor,
                            t_cor,
                            chol = FALSE,
                            h_options = NULL,
                            ...) {
  # the deafult simulation method
  # the user does not have to provide a covariance matrix
  # setting default h options if none are provided
  if (is.null(h_options)) {
    h_options <- list(
      h_large = TRUE,
      h_t_distmetric = "euclidean",
      h_s_distmetric = "euclidean"
    )
  }
  
  # simulation the random vector by cholesky decompositing
  # the full covariance matrix
  if (chol) {
    
    # compute the large spatial distance matrices in 1d
    if (is.null(ycoord)) {
      
      # compute the large spatial distance matrix
      h_s_large <- make_h(
        coord1 = data[[xcoord]],
        distmetric = h_options$h_s_distmetric
      )
      
      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = data[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
    } else { # compute the large spatial distance matrices in 2d
      
      # compute the large spatial distance matrix
      h_s_large <- make_h(
        coord1 = data[[xcoord]],
        coord2 = data[[ycoord]],
        distmetric = h_options$h_s_distmetric
      )
      
      # compute the large temporal distance matrix
      h_t_large <- make_h(
        coord1 = data[[tcoord]],
        distmetric = h_options$h_t_distmetric
      )
    }
    
    # make the spatio-temporal covariance
    st_covariance <- make_stcovariance(
      covparam_object = object,
      h_s_large = h_s_large,
      h_t_large = h_t_large,
      s_cor = s_cor,
      t_cor = t_cor
    )
    
    # simulate the random vector by calling strnorm.matrix
    strnorm_sim <- strnorm.matrix(
      object = st_covariance,
      mu = mu,
      size = size,
      condition = condition
    )
  } else {
    # rename the object as covparam_object (should provide a check for this later)
    covparam_object <- object
    
    # storing an original key of provided spatio-temporal locations
    data$original_key <- seq.int(1, nrow(data))
    
    # order the data by space within time
    spint <- storder(
      data = data,
      xcoord = xcoord,
      ycoord = ycoord,
      tcoord = tcoord,
      h_options = h_options
    )
    
    # compute the small spatial correlation matrix
    r_s_small <- make_r(
      h = spint$h_s_small,
      range = covparam_object[["s_range"]],
      structure = s_cor
    )
    
    # adding a small diagonal constant for inversion stability
    diag(r_s_small) <- diag(r_s_small) + condition
    
    # compute the small temporal correlation matrix
    r_t_small <- make_r(
      h = spint$h_t_small,
      range = covparam_object[["t_range"]],
      structure = t_cor
    )
    
    # adding a small diagonal constant for inversion stability
    diag(r_t_small) <- diag(r_t_small) + condition
    
    
    # simulating the random vector
    strnorm_sim <- strnorm_small(
      covparam_object = covparam_object,
      mu = mu,
      size = size,
      r_s_small = r_s_small,
      r_t_small = r_t_small,
      error = error
    )
    
    # removing the spatio-temporal observations not provided
    strnorm_sim <- strnorm_sim[spint$ordered_data_dense$observed, , drop = FALSE]
    
    # ordering by the original data
    strnorm_sim <- strnorm_sim[order(spint$ordered_data_o$original_key), , drop = FALSE]
    
    # adding the mean
    strnorm_sim <- mu + strnorm_sim
    
    # return the random vector
    return(strnorm_sim)
  }
}

# the function that actually simulates the random vector in strnorm.deafult
# and object is a covariance parameter object
strnorm_small <- function(covparam_object, mu, size, r_s_small, r_t_small, error) {
  
  # calling the appropriate generic
  UseMethod("strnorm_small", object = covparam_object)
}


strnorm_small.productsum <- function(covparam_object, mu, size, r_s_small, r_t_small, error) {
  
  # computing the lower spatial cholesky
  chol_r_s_small <- t(chol(r_s_small))
  
  # computing the lower temporal cholesky
  chol_r_t_small <- t(chol(r_t_small))
  
  # storing the number of unique spatial locations
  n_s <- nrow(chol_r_s_small)
  
  # storing the number of unique temporal locations
  n_t <- nrow(chol_r_t_small)
  
  # storing the number of possible spatio-temporal locations
  n_st <- n_s * n_t
  
  # multiply the spatial cholesky by zs on the left
  zs_chol_r_s_small <- multiply_z(
    mx = chol_r_s_small,
    z_type = "spatial",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )
  
  # multiply the temporal cholesky by zt on the left
  zt_chol_r_t_small <- multiply_z(
    mx = chol_r_t_small,
    z_type = "temporal",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )
  
  # computing the lower spatio-temporal cholesky
  chol_r_st <- kronecker(chol_r_t_small, chol_r_s_small)
  
  # simulate the random error vector
  strnorm_small_sim <- vapply(
    1:size,
    function(x) {
      # spatial dependent error simulation
      s_de_sim <- zs_chol_r_s_small %*% rnorm(n_s, sd = sqrt(covparam_object[["s_de"]]))
      
      # spatial independent error simulation
      s_ie_sim <- multiply_z(
        mx = rnorm(n_s, sd = sqrt(covparam_object[["s_ie"]])),
        z_type = "spatial",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )
      
      # temporal dependent random error simulation
      t_de_sim <- zt_chol_r_t_small %*% rnorm(n_t, sd = sqrt(covparam_object[["t_de"]]))
      
      # temporal independent error simulation
      t_ie_sim <- multiply_z(
        mx = rnorm(n_t, sd = sqrt(covparam_object[["t_ie"]])),
        z_type = "temporal",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )
      
      # spatio-temporal dependent error simulation
      st_de_sim <- chol_r_st %*% rnorm(n_st, sd = sqrt(covparam_object[["st_de"]]))
      
      # spatio-temporal independent error simulation
      st_ie_sim <- rnorm(n_st, sd = sqrt(covparam_object[["st_ie"]]))
      
      # simulating the random error vector
      if (error == "normal") {
        randomerror_sim <- s_de_sim + s_ie_sim + t_de_sim + t_ie_sim + st_de_sim + st_ie_sim
        return(randomerror_sim)
      } else if (error == "component_squared") {
        normalerror_sim <- s_de_sim + s_ie_sim + t_de_sim + t_ie_sim + st_de_sim + st_ie_sim
        var_normalerror_sim <- var(as.vector(normalerror_sim))
        squarederror_sim <- sign(s_de_sim) * s_de_sim^2 +
          sign(s_ie_sim) * s_ie_sim^2 +
          sign(t_de_sim) * t_de_sim^2 +
          sign(t_ie_sim) * t_ie_sim^2 +
          sign(st_de_sim) * st_de_sim^2 +
          sign(st_ie_sim) * st_ie_sim^2
        var_squarederorr_sim <- var(as.vector(squarederror_sim))
        randomerror_sim <- squarederror_sim * sqrt(var_normalerror_sim) / sqrt(var_squarederorr_sim)
        return(randomerror_sim)
      } else if (error == "sum_squared") {
        normalerror_sim <- s_de_sim + s_ie_sim + t_de_sim + t_ie_sim + st_de_sim + st_ie_sim
        var_normalerror_sim <- var(as.vector(normalerror_sim))
        squarederror_sim <- sign(normalerror_sim) * normalerror_sim^2
        var_squarederorr_sim <- var(as.vector(squarederror_sim))
        randomerror_sim <- squarederror_sim * sqrt(var_normalerror_sim) / sqrt(var_squarederorr_sim)
        return(randomerror_sim)
      } else {
        stop("choose a valid error type")
      }
    },
    double(n_st)
  )
  
  # return the random error vector
  return(strnorm_small_sim)
}

strnorm_small.sum_with_error <- function(covparam_object, mu, size, r_s_small, r_t_small) {
  
  # storing the lower spatial cholesky
  chol_r_s_small <- t(chol(r_s_small))
  
  # storing the lower temporal cholesky
  chol_r_t_small <- t(chol(r_t_small))
  
  # storing the number of unique spatial locations
  n_s <- nrow(chol_r_s_small)
  
  # storing the number of unique temporal locations
  n_t <- nrow(chol_r_t_small)
  
  # storing the number of possible spatio-temporal locations
  n_st <- n_s * n_t
  
  # multiply the spatial cholesky by zs on the left
  zs_chol_r_s_small <- multiply_z(
    mx = chol_r_s_small,
    z_type = "spatial",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )
  
  # multiply the spatial cholesky by zs on the left
  zt_chol_r_t_small <- multiply_z(
    mx = chol_r_t_small,
    z_type = "temporal",
    n_s = n_s,
    n_t = n_t,
    side = "left"
  )
  
  # simulate the random error vector
  strnorm_small_sim <- vapply(
    1:size,
    function(x) {
      
      # spatial dependent error simulation
      s_de_sim <- zs_chol_r_s_small %*% rnorm(n_s, sd = sqrt(covparam_object[["s_de"]]))
      
      # spatial independent error simulation
      s_ie_sim <- multiply_z(
        mx = rnorm(n_s, sd = sqrt(covparam_object[["s_ie"]])),
        z_type = "spatial",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )
      
      # temporal dependent error simulation
      t_de_sim <- zt_chol_r_t_small %*% rnorm(n_t, sd = sqrt(covparam_object[["t_de"]]))
      
      # temporal independent error simulation
      t_ie_sim <- multiply_z(
        mx = rnorm(n_t, sd = sqrt(covparam_object[["t_ie"]])),
        z_type = "temporal",
        n_s = n_s,
        n_t = n_t,
        side = "left"
      )
      
      # spatio-temporal independent error simulation
      st_ie_sim <- rnorm(n_st, sd = sqrt(covparam_object[["st_ie"]]))
      
      # simulating the random error vector
      randomerror_sim <- s_de_sim + s_ie_sim + t_de_sim + t_ie_sim + st_ie_sim
      return(randomerror_sim)
    },
    double(n_st)
  )
  
  # return the random error vector
  return(strnorm_small_sim)
}


strnorm_small.product <- function(covparam_object, mu, size, r_s_small, r_t_small) {
  
  # store the lower spatial cholesky
  chol_r_s_small <- t(chol(r_s_small))
  
  # store the lower temporal cholesky
  chol_r_t_small <- t(chol(r_t_small))
  
  # storing the number of unique spatial locations
  n_s <- nrow(chol_r_s_small)
  
  # storing the number of unique temporal locations
  n_t <- nrow(chol_r_t_small)
  
  # storing the number of possible spatio-temporal locations
  n_st <- n_s * n_t
  
  # store the lower spatio-temporal cholesky
  chol_r_st <- kronecker(chol_r_t_small, chol_r_s_small)
  
  # simulate the random error vector
  strnorm_small_sim <- vapply(
    1:size,
    function(x) {
      
      # spatio-temporal dependent error simulation
      st_de_sim <- chol_r_st %*% rnorm(n_st, sd = sqrt(covparam_object[["st_de"]]))
      
      # simulating the random error vector
      randomerror_sim <- st_de_sim
      return(randomerror_sim)
    },
    double(n_st)
  )
  
  # return the random error vector
  return(strnorm_small_sim)
}

#' Summarize a Spatio-Temporal Linear Mixed Model
#'
#' @param object A \code{stlmm} object.
#' @param ... Additional arguments.
#'
#' @name summary
#'
#' @method summary stlmm
#'
#' @return A list containing several objects
#'   \describe{
#'     \item{\code{Call}}{The original function call.}
#'     \item{\code{FixedEffects}}{Fixed effects estimates and standard errors.}
#'     \item{\code{CovarianceParameters}}{Covariance parameter estimates.}
#'     \item{\code{CovarianceForms}}{The spatial, temopral, and spatio-temporal correlation forms.}
#'     \item{\code{Residuals}}{Raw residuals.}
#'     \item{\code{ObjectiveFn}}{Objective function value.}
#'   }
#' @noRd

summary.stlmm <- function(object, ...) {
  
  # store the formula
  call <- object$formula
  
  # store the coefficients
  regcoefs <- as.vector(object$Coefficients)
  
  # store the covariance of beta hat
  regvar <- object$CovCoefficients
  
  # store the number of fixed effect predictors
  p <- ncol(object$model$FixedDesignMatrix)
  
  # store the number of observations
  n <- nrow(object$model$FixedDesignMatrix)
  
  # store the standard errors
  sereg <- sqrt(diag(as.matrix(regvar)))
  
  # store a vector of t statistics
  tvec <- regcoefs / sereg
  
  # store a vector of p-values appropriately rounded
  pvec <- round(100000 * (1 - pt(abs(regcoefs / sereg), df = n - p)) * 2) / 100000
  
  # save the fixed effect data frame
  fixed.effect.estimates <- data.frame(
    Estimate = regcoefs,
    Std.Error = sereg,
    t.value = tvec,
    prob.t = pvec
  )
  
  # assign row names
  rownames(fixed.effect.estimates) <- object$NamesCoefficients
  
  # save the covariance parameters as a list and then data frame
  covmodels <- as.list(object$CovarianceParameters)
  covmodelout <- data.frame(covmodels, stringsAsFactors = FALSE)
  
  # save the covariance information as a list and then data frame
  covinfo <- as.list(object$CovarianceForms)
  covinfoout <- data.frame(covinfo, stringsAsFactors = FALSE)
  
  # save the residuals
  resid_vec <- object$Residuals
  
  # save the objective functions
  objective <- object$Objective
  
  # store the final summary output
  output <- structure(
    list(
      Call = call,
      FixedEffects = fixed.effect.estimates,
      CovarianceParameters = covmodelout,
      CovarianceForms = covinfoout,
      Residuals = resid_vec,
      ObjectiveFn = objective
    ),
    class = "summary.stlmm"
  )
  
  # return the summary output
  return(output)
}

varest <- function(invert_object, invert_output) {
  
  # calling the appropriate generic
  UseMethod("varest", object = invert_object)
}

# compute the overall variance for reml estimation
varest.reml <- function(invert_object, invert_output) {
  
  # store the sample size
  n <- nrow(invert_output$sigmainv_o)
  
  # compute the beta estimate and return estmation information (via return_estlist)
  betaest_output <- betaest(
    xo = invert_object$xo,
    sigmainv_xyo = invert_output$sigmainv_o,
    condition = invert_object$condition,
    return_estlist = TRUE
  )
  
  # compute sigma inverse times the residual vector
  siginv_r <- betaest_output$estlist$sigmainv_yo - betaest_output$estlist$sigmainv_xo %*% betaest_output$betahat
  
  # compute the quadratic residual function
  r_siginv_r <- t(invert_object$yo - invert_object$xo %*% betaest_output$betahat) %*% siginv_r
  
  # compute the overall variance
  ws_l2 <- as.vector(r_siginv_r / (n - betaest_output$estlist$p))
  
  # return the overall variance
  return(ws_l2) # thanks Russ and John!
}

r2plo <- function(covparam_object, ...) {
  # calling the appropriate estimation method / stcov type generic (class/sub)
  UseMethod("r2plo", object = covparam_object)
}

# the regular to profiled log odds for semivariogram-weighted least squares
r2plo.svwls <- function(covparam_object, ...) {
  UseMethod("r2plo.svwls", object = covparam_object)
}

# the svwls productsum r2plo
r2plo.svwls.productsum <- function(covparam_object, max_options) {
  
  # giving covparam_object a new name
  params <- covparam_object
  
  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  
  # computing the overall spatial variance
  svar <- sum(s_vc)
  
  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  
  # computing the overall temporal variance
  tvar <- sum(t_vc)
  
  # storing the spatio-temporal variance components
  st_vc <- c(params[["st_de"]], params[["st_ie"]])
  
  # computing the overall spatio-temporal variance
  stvar <- sum(st_vc)
  
  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)
  
  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  
  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)
  
  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar
  
  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar
  
  # computing the proportion of spatial-temporal that is independent error
  n_st <- params[["st_ie"]] / stvar
  
  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)
  
  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    var_prop = pmin(1, (svar + tvar + stvar) / max_options$max_v),
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )
  
  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))
  
  # returning the logit of the profiled parameters
  return(pparm)
}

# the svwls sum with error r2plo
r2plo.svwls.sum_with_error <- function(covparam_object, max_options) {
  
  # giving covparam_object a new name
  params <- covparam_object
  
  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  
  # computing the overall spatial variance
  svar <- sum(s_vc)
  
  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  
  # computing the overall temporal variance
  tvar <- sum(t_vc)
  
  # computing the spatio-temporal independent error
  st_vc <- params[["st_ie"]]
  
  # restoring it to make it clear this is the only parameter
  stvar <- sum(st_vc)
  
  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)
  
  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  
  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)
  
  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar
  
  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar
  
  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)
  
  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    var_prop = pmin(1, (svar + tvar + stvar) / max_options$max_v),
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )
  
  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))
  
  # returning the logit of the profiled parameters
  return(pparm)
}

# the svwls product r2plo
r2plo.svwls.product <- function(covparam_object, max_options) {
  
  # giving covparam_object a new name
  params <- covparam_object
  
  # storing the profiled parameters
  pparm <- c(v_s = params[["v_s"]], v_t = params[["v_t"]])
  
  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    var_prop = pmin(1, params[["st_de"]] / max_options$max_v),
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )
  
  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))
  
  # returning the logit of the profiled parameters
  return(pparm)
}

# the regular to profiled log odds for semivariogram-weighted least squares
r2plo.reml <- function(covparam_object, ...) {
  UseMethod("r2plo.reml", object = covparam_object)
}

# the reml productsum r2plo
r2plo.reml.productsum <- function(covparam_object, max_options) {
  # giving covparam_object a new name
  params <- covparam_object
  
  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  
  # computing the overall spatial variance
  svar <- sum(s_vc)
  
  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  
  # computing the overall temporal variance
  tvar <- sum(t_vc)
  
  # storing the spatio-temporal variance components
  st_vc <- c(params[["st_de"]], params[["st_ie"]])
  
  # computing the overall spatio-temporal variance
  stvar <- sum(st_vc)
  
  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)
  
  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  
  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)
  
  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar
  
  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar
  
  # computing the proportion of spatial-temporal that is independent error
  n_st <- params[["st_ie"]] / stvar
  
  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t, n_st = n_st)
  
  # storing the proportion of spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )
  
  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))
  
  # returning the logit of the profiled parameters
  return(pparm)
}

# the reml sum with error r2plo
r2plo.reml.sum_with_error <- function(covparam_object, max_options) {
  
  # giving covparam_object a new name
  params <- covparam_object
  
  # storing the spatial variance components
  s_vc <- c(params[["s_de"]], params[["s_ie"]])
  
  # computing the overall spatial variance
  svar <- sum(s_vc)
  
  # storing the temporal variance components
  t_vc <- c(params[["t_de"]], params[["t_ie"]])
  
  # computing the overall temporal variance
  tvar <- sum(t_vc)
  
  # computing the spatio-temporal independent error
  st_vc <- params[["st_ie"]]
  
  # restoring it to make it clear this is the only parameter
  stvar <- sum(st_vc)
  
  # storing the overall variances
  vc <- c(s_vc, t_vc, st_vc)
  
  # computing the proportion of main effect variance
  lambda <- (svar + tvar) / (svar + tvar + stvar)
  
  # computing the proportion of main effect variance that is spatial
  alpha <- svar / (svar + tvar)
  
  # computing the proportion of spatial variance that is independent error
  n_s <- params[["s_ie"]] / svar
  
  # computing the proportion of temporal variance that is independent error
  n_t <- params[["t_ie"]] / tvar
  
  # storing the profiled parameters
  pparm <- c(lambda = lambda, alpha = alpha, n_s = n_s, n_t = n_t)
  
  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )
  
  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))
  
  # returning the logit of the profiled parameters
  return(pparm)
}

# the reml product r2plo
r2plo.reml.product <- function(covparam_object, max_options) {
  
  # giving covparam_object a new name
  params <- covparam_object
  
  # storing the profiled parameters
  pparm <- c(v_s = params[["v_s"]], v_t = params[["v_t"]])
  
  # storing the proportion of variance/spatial range/temporal range relative
  # to their maximums
  pparm <- c(
    pparm,
    srange_prop = params[["s_range"]] / max_options$max_s_range,
    trange_prop = params[["t_range"]] / max_options$max_t_range
  )
  
  # storing the logit of the profiled parameters
  pparm <- log(pparm / (1 - pparm))
  
  # returning the logit of the profiled parameters
  return(pparm)
}

plo2r <- function(par, covest_object, ...) {
  # calling the appropriate estimation method / stcov type generic (class/sub)
  UseMethod("plo2r", object = covest_object)
}

# the profiled log odds to regular for semivariogram-weighted least squares
plo2r.svwls <- function(par, covest_object) {
  
  # calling the appropriate estimation method generic
  UseMethod("plo2r.svwls", object = covest_object)
}

# the svwls productsum plo2r
plo2r.svwls.productsum <- function(par, covest_object) {
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] <- 7
  par[par < -7] <- -7
  
  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))
  
  # store lambda, alpha, n_s, n_t, n_st
  lambda <- invlogit[["lambda"]] # proportion of main effect variance
  alpha <- invlogit[["alpha"]] # proportion of spatial main effect variance
  n_s <- invlogit[["n_s"]] # proportion of spatial main efffect ind error
  n_t <- invlogit[["n_t"]] # proportion of temporal main effect ind error
  n_st <- invlogit[["n_st"]] # proportion of interaction ind error
  
  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_de <- (1 - lambda) * (1 - n_st)
  st_ie <- (1 - lambda) * n_st
  # overall variance capped by max_v
  ov_var <- covest_object$max_options$max_v * invlogit[["var_prop"]]
  
  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie
  ),
  # range capped by max_s_range
  s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
  # range capped by max_t_range
  t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  
  # returning the parameters
  return(rparm)
}


plo2r.svwls.sum_with_error <- function(par, covest_object) {
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] <- 7
  par[par < -7] <- -7
  
  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))
  
  # store lambda, alpha, n_s, n_t
  lambda <- invlogit[["lambda"]]
  alpha <- invlogit[["alpha"]]
  n_s <- invlogit[["n_s"]]
  n_t <- invlogit[["n_t"]]
  # n_st is automatically 1
  
  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_ie <- (1 - lambda)
  
  # overall variance capped by max_v
  ov_var <- covest_object$max_options$max_v * invlogit[["var_prop"]]
  
  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_ie = st_ie
  ),
  # range parameter capped by max_s_range
  s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
  # range parameter capepd by max_t_range
  t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  # returning the parameters
  return(rparm)
}

plo2r.svwls.product <- function(par, covest_object) {
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] <- 7
  par[par < -7] <- -7
  
  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))
  
  # store v_s and v_t
  v_s <- invlogit[["v_s"]]
  v_t <- invlogit[["v_t"]]
  
  # overall variance capped by max_v
  ov_var <- covest_object$max_options$max_v * invlogit[["var_prop"]]
  
  # storing the parameter vector
  rparm <- c(
    v_s = v_s,
    v_t = v_t,
    st_de = ov_var,
    # range parameter capped by max_s_range
    s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
    # range parameter capped by max_t_range
    t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  
  # return the variance parameters
  return(rparm)
}


# profiled log odds generic for reml
plo2r.reml <- function(par, covest_object, ...) {
  UseMethod("plo2r.reml", object = covest_object)
}

plo2r.reml.productsum <- function(par, covest_object, ov_var) {
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] <- 7
  par[par < -7] <- -7
  
  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))
  
  # store lambda, alpha, n_s, n_t, n_st
  lambda <- invlogit[["lambda"]] # proportion of main effect variance
  alpha <- invlogit[["alpha"]] # proportion of spatial main effect variance
  n_s <- invlogit[["n_s"]] # proportion of spatial main efffect ind error
  n_t <- invlogit[["n_t"]] # proportion of temporal main effect ind error
  n_st <- invlogit[["n_st"]] # proportion of interaction ind error
  
  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_de <- (1 - lambda) * (1 - n_st)
  st_ie <- (1 - lambda) * n_st
  
  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = st_de,
    st_ie = st_ie
  ),
  # range capped by max_s_range
  s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
  # range capped by max_t_range
  t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  
  # returning the parameters
  return(rparm)
}


plo2r.reml.sum_with_error <- function(par, covest_object, ov_var) {
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] <- 7
  par[par < -7] <- -7
  
  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))
  
  # store lambda, alpha, n_s, n_t
  lambda <- invlogit[["lambda"]]
  alpha <- invlogit[["alpha"]]
  n_s <- invlogit[["n_s"]]
  n_t <- invlogit[["n_t"]]
  # n_st is automatically 1
  
  # transform the plo parameters to regular parameters
  s_de <- lambda * alpha * (1 - n_s)
  s_ie <- lambda * alpha * n_s
  t_de <- lambda * (1 - alpha) * (1 - n_t)
  t_ie <- lambda * (1 - alpha) * n_t
  st_ie <- (1 - lambda)
  
  # storing the parameter vector
  rparm <- c(ov_var * c(
    s_de = s_de,
    s_ie = s_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_ie = st_ie
  ),
  # range parameter capped by max_s_range
  s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
  # range parameter capepd by max_t_range
  t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  # returning the parameters
  return(rparm)
}

plo2r.reml.product <- function(par, covest_object, ov_var) {
  
  # set reasonable bounds on what the exponentiated plo parameters can be - this was a frustrating bug to find
  par[par > 7] <- 7
  par[par < -7] <- -7
  
  # inverse logit the parameter vector
  invlogit <- exp(par) / (1 + exp(par))
  
  # store v_s and v_t
  v_s <- invlogit[["v_s"]]
  v_t <- invlogit[["v_t"]]
  
  # storing the parameter vector
  rparm <- c(
    v_s = v_s,
    v_t = v_t,
    st_de = ov_var,
    # range parameter capped by max_s_range
    s_range = covest_object$max_options$max_s_range * invlogit[["srange_prop"]],
    # range parameter capped by max_t_range
    t_range = covest_object$max_options$max_t_range * invlogit[["trange_prop"]]
  )
  
  # return the variance parameters
  return(rparm)
}