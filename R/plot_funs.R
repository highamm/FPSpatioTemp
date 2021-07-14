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
