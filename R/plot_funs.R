#' Covariance Plot
#' 
#' Produces a plot of the covariance with time distance on the x-axis for various
#' degrees of spatial distance.
#' 
#' @param stlmfit_obj an object fit with \code{stlmfit()}
#' @param sp_epstol a vector of the values used for colours for various spatial distances. The default is the a vector of \code{0}, the minimum spatial distance between two sites, the maximum spatial distance between two sites, and \code{Inf}.
#' @param t_max the maximum value for temporal distance to be plotted on the x-axis. The default is the maximum temporal distance between two time points in the data.
#' @param ... extra options to be passed to `\code{ggplot}`
#' 
#' @return a plot with covariance on the y-axis and temporal distance on the x-axis
#' coloured by various spatial distances in \code{sp_epstol}.
#' @examples 
#' obj <- stlmfit(formula = response_na ~ x, data = samp_data,
#'  xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
#' plot_cov(obj)
#' plot_cov(stlmfit_obj = obj, sp_epstol = c(0.2, 0.4, 1.2, 2, Inf),
#' t_max = 1.5)
#' @import ggplot2
#' @export plot_cov

plot_cov <- function(stlmfit_obj, sp_epstol = NULL,
                         t_max = NULL, ... ) {
  
  if (is.null(sp_epstol) == TRUE) {
    sp_epstol <- c(stlmfit_obj$minimax_vec["min_dist_sp"],
                  stlmfit_obj$minimax_vec["max_dist_sp"])
  }
  
  if(is.null(t_max) == TRUE) {
    t_max <- stlmfit_obj$minimax_vec["max_dist_t"]
  }

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

  
  h_sp_seq <- seq(0, sp_rangetol, length.out = 200)
  h_sp_seq <- append(h_sp_seq, sp_epstol, sp_range / 2)
  h_t_seq <- seq(0, t_max, length.out = 200)
  h_t_seq <- append(h_t_seq, t_range / 2)
  h_sp <- rep(h_sp_seq, times = length(h_t_seq))
  h_t <- rep(h_t_seq, each = length(h_sp_seq))
  
  # make covparam object
  covparams <- make_covparam_object(
    s_de = sp_de,
    s_ie = sp_ie,
    t_de = t_de,
    t_ie = t_ie,
    st_de = spt_de,
    st_ie = spt_ie,
    s_range = sp_range, 
    t_range = t_range,
    stcov = "productsum"
  )
  
  names(covparams) <- c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie",
                        "s_range", "t_range")

  
  sigma <- make_stcovariance(covparams, h_sp, h_t,
                      stlmfit_obj$summary_stlmm$CovarianceForms[["s_cor"]],
                          stlmfit_obj$summary_stlmm$CovarianceForms[["t_cor"]])
  
  
  data <- tibble::tibble(
    h_s = h_sp,
    h_t = h_t,
    sigma = sigma
  )
  
  t_plot_pos <- data |>
    dplyr::filter((.data$h_s %in% c(0, sp_epstol)) & (.data$h_t > 0.02))
  t_plot_zero <- data |>
    dplyr::filter((.data$h_s %in% c(0, sp_epstol)) & (.data$h_t == 0))
  
    
    ggplot(data = t_plot_pos) +
      geom_line(mapping = aes(x = h_t, y = sigma, colour = as.factor(round(.data$h_s, digits = 2))), linewidth = 1.5, ...) +
      geom_point(data = t_plot_zero, mapping = aes(x = h_t, y = sigma, colour = as.factor(round(.data$h_s, digits = 2))), size = 4) +
      scale_colour_viridis_d(name = "Spatial Distance", begin = 0, end = 0.9) +
      theme_minimal() +
      labs(x = "Temporal Distance",
           y = "Estimated Covariance")
}








