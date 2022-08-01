#' Spatial Plot
#' 
#' Produces a spatial plot of a response variable faceted by time.
#' @param xcoords a vector of x-coordinates
#' @param ycoords a vector of y-coordinates
#' @param times a vector of times
#' @param response, a vector of the response variable
#' 
#' @return a plot.
#' @examples 
#' sim_data <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'        sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'        sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'        sigma_nugget_spacetime = 0.5)
#'  df <- sim_data$out_df
#'  plot_sp(df$xcoords, df$ycoords, df$times, df$response)
#'  
#' @import ggplot2
#' @export plot_sp 

plot_sp <- function(xcoords, ycoords, times, response) {
  df <- dplyr::tibble(xcoords = xcoords, ycoords = ycoords, times = times,
                response = response)
  sp_plot <- ggplot2::ggplot(data = df,
                             ggplot2::aes(x = xcoords, y = ycoords,
                                            colour = response)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ times) +
    ggplot2::scale_colour_viridis_c()
  return(sp_plot)
}

#' Temporal Plot
#' 
#' Produces a temporal plot of a response variable grouped by location
#' @param xcoords a vector of x-coordinates
#' @param ycoords a vector of y-coordinates
#' @param times a vector of times
#' @param response, a vector of the response variable
#' @param alpha a line transparency parameter
#' 
#' @return a plot.
#' @examples 
#' sim_data <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'        sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'        sigma_parsil_time = 1, rho = 0.99, sigma_nugget_time = 0,
#'        sigma_nugget_spacetime = 0.5)
#'  df <- sim_data$out_df
#'  plot_t(df$xcoords, df$ycoords, df$times, df$response)
#' @import ggplot2
#' @export plot_t

plot_t <- function(xcoords, ycoords, times, response, alpha = 0.4) {
  df <- dplyr::tibble(xcoords = xcoords, ycoords = ycoords, times = times,
                      response = response)
  ggplot2::ggplot(data = df, ggplot2::aes(x = times, y = response)) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(xcoords, ycoords)),
              alpha = alpha)
}


#' Covariance Plot
#' 
#' Produces a plot of the covariance with time distance on the x-axis for various
#' degrees of spatial distance.
#' 
#' @param stlmfit_obj an object fit with \code{stlmfit()}
#' @param sp_epstol a range of values used for colours for various spatial distances
#' @param t_epstol a range of values used for colours for various temporal distances
#' @param xaxis_var is \code{"spatial"} if the spatial distances should go on the 
#' x-axis and \code{"temporal"} if the temporal distances should go on the x-axis
#' 
#' @return a plot with covariance on the y-axis and temporal distance on the x-axis
#' coloured by various spatial distances in \code{sp_epstol}.
#' @examples 
#' set.seed(07262022)
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sp_de = 0.5, sp_range = 4, sp_ie = 0.5,
#'       t_de = 0.5, t_range = 0.7, t_ie = 0.5,
#'       spt_ie = 0.5)
#'       
#' samp_obj <- sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' samp_data <- samp_obj$df_full
#' samp_data <- samp_data |>
#'  dplyr::mutate(predwts = dplyr::if_else(times == max(times),
#'   true = 1, false = 0))
#' samp_data$x <- rnorm(nrow(samp_data), 0, 1)
#' samp_data <- samp_data
#' stlmfit_obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#'  
#' plot_cov(stlmfit_obj = stlmfit_obj, sp_epstol = c(0.01, 0.1, 0.5, Inf),
#' t_epstol = c(0.2, 1, 4))

plot_cov <- function(stlmfit_obj, sp_epstol = c(0.2, 4, 20, Inf),
                         t_epstol = c(0.2, 2, 6),
                         xaxis_var = "temporal") {
  
  cov_parms <- stlmfit_obj$cov_parms
  
  sp_de <- cov_parms[1]
  sp_ie <- cov_parms[2]
  t_de <- cov_parms[4]
  t_ie <- cov_parms[5]
  spt_de <- cov_parms[7]
  spt_ie <- cov_parms[8]
  total_var <- sum(sp_de, sp_ie, t_de, t_ie, spt_de, spt_ie)
  sp_range <- cov_parms[3]
  t_range <- cov_parms[6]
  
  rangetol <- 1 / 4
  sp_rangetol <- sp_range + sp_range * rangetol
  # t_rangetol <- t_range + t_range * rangetol
  t_rangetol <- max(t_epstol) ## for the six years
  
  h_sp_seq <- seq(0, sp_rangetol, length.out = 200)
  h_sp_seq <- append(h_sp_seq, sp_epstol, sp_range / 2)
  h_t_seq <- seq(0, t_rangetol, length.out = 200)
  h_t_seq <- append(h_t_seq, t_range / 2)
  h_sp <- rep(h_sp_seq, times = length(h_t_seq))
  h_t <- rep(h_t_seq, each = length(h_sp_seq))
  
  # make covparam object
  covparams <- DumelleEtAl2021STLMM::make_covparam_object(
    s_de = sp_de,
    s_ie = sp_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = spt_de,
    st_ie = spt_ie,
    s_range = sp_range * 3, ## dumelle uses effective range
    t_range = t_range * 3, ## dumelle uses effective range
    stcov = "productsum"
  )
  
  names(covparams) <- c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie",
                        "s_range", "t_range")
  
  # make stcovariance
  sigma <- DumelleEtAl2021STLMM::make_stcovariance(covparams, h_sp, h_t,
                                                   "exponential", "exponential")
  
  gamma <- DumelleEtAl2021STLMM::make_stsemivariogram(covparams, h_sp, h_t,
                                                      "exponential",
                                                      "exponential")
  
  data <- tibble::tibble(
    h_s = h_sp,
    h_t = h_t,
    sigma = sigma,
    gamma = gamma
  )
  
  t_plot_pos <- data |>
    dplyr::filter((h_sp %in% c(0, sp_epstol)) & (h_t > 0.2))
  t_plot_zero <- data |>
    dplyr::filter((h_sp %in% c(0, sp_epstol)) & (h_t == 0))
  
  sp_plot_pos <- data |> 
    dplyr::filter(h_s > 0 & h_t %in% c(0, min(h_t[h_t > 0]), max(h_t)))
  sp_plot_zero <- data |>
    dplyr::filter(h_s == 0 & h_t %in% c(0, min(h_t[h_t > 0]), max(h_t)))
  
  if (xaxis_var == "temporal") {
    
    ggplot(data = t_plot_pos) +
      geom_line(mapping = aes(x = h_t, y = sigma, colour = as.factor(h_s)), size = 1.5) +
      geom_point(data = t_plot_zero, mapping = aes(x = h_t, y = sigma, colour = as.factor(h_s)), size = 4) +
      scale_colour_viridis_d(name = "Spatial Distance", begin = 0, end = 0.9) +
      theme_minimal() +
      labs(x = "Temporal Distance (years)",
           y = "Estimated Covariance")
    
  } else if (xaxis_var == "spatial") {
    
    ggplot(data = sp_plot_pos) +
      geom_line(mapping = aes(x = h_s, y = sigma, colour = as.factor(h_t)), size = 1.5) +
      geom_point(data = sp_plot_zero, mapping = aes(x = h_s, y = sigma, colour = as.factor(h_t)), size = 4) +
      scale_colour_viridis_d(end = 0.9) +
      theme_minimal()  +
      labs(x = "Spatial Distance",
           y = "Estimated Covariance")
    
  } else {
    
    stop("xaxis_var must be either 'spatial' or 'temporal'")
    
  }
  
  # plot_df_obj <- list(t_plot_pos, t_plot_zero, sp_plot_pos, sp_plot_zero)
  # class(plot_df_obj) <- "make_plot_df"
  # 
  # return(plot_df_obj)
}








