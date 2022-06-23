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
#' @param xcoords a vector of x-coordinates
#' @param ycoords a vector of y-coordinates
#' @param times a vector of times
#' @param response, a vector of the response variable
#' @param alpha a line transparency parameter


moose_parms <- fit_obj_all$parms ## (from manu.Rmd file)

s_de <- moose_parms$sigma_parsil_spat
s_ie <- moose_parms$sigma_nugget_spat
t_de <- moose_parms$sigma_parsil_time
t_ie <- moose_parms$sigma_nugget_time
st_de <- moose_parms$sigma_parsil_spacetime
st_ie <- moose_parms$sigma_nugget_spacetime
total_var <- sum(s_de, s_ie, t_de, t_ie, st_de, st_ie)
s_range <- moose_parms$range
t_range <- moose_parms$rho
rangetol <- 1 / 4
s_rangetol <- s_range + s_range * rangetol
# t_rangetol <- t_range + t_range * rangetol
t_rangetol <- 3 * 2


brspace <- 0.075
brlen <- 0.15

# epsilon tolerance
epstol <- 1e-5
# axistext_size <- 24
# legendtext_size <- 24
# annotate_size <- 9

s_epstol <- c(0.2, 4, 20, Inf) ##  4, 10,


# setting unique distance values
h_s_seq <- seq(0, s_rangetol, length.out = 200)
h_s_seq <- append(h_s_seq, s_epstol, s_range / 2)
h_t_seq <- seq(0, t_rangetol, length.out = 200)
h_t_seq <- append(h_t_seq, t_range / 2)
h_s <- rep(h_s_seq, times = length(h_t_seq))
h_t <- rep(h_t_seq, each = length(h_s_seq))

t_epstol <- c(0.2, 2, max(h_t))

# make covparam object
covparams <- DumelleEtAl2021STLMM::make_covparam_object(
  s_de = s_de,
  s_ie = s_ie,
  t_de = t_de,
  t_ie = t_ie,
  st_de = st_de,
  st_ie = st_ie,
  s_range = s_range,
  t_range = t_range,
  stcov = "productsum"
)
names(covparams) <- c("s_de", "s_ie", "t_de", "t_ie", "st_de", "st_ie",
                      "s_range", "t_range")

# make stcovariance
sigma <- DumelleEtAl2021STLMM::make_stcovariance(covparams, h_s, h_t, "exponential", "exponential")

gamma <- DumelleEtAl2021STLMM::make_stsemivariogram(covparams, h_s, h_t,
                                                    "exponential",
                                                    "exponential")

data <- tibble::tibble(
  h_s = h_s,
  h_t = h_t,
  sigma = sigma,
  gamma = gamma
)

# temporal plotting subset
# change s_epstol to c(1, Inf)
t_plot_pos <- data |>
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t > 0.2))
t_plot_zero <- data |>
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t == 0))

s_plot_pos <- data |> 
  dplyr::filter(h_s > 0 & h_t %in% c(0, min(h_t[h_t > 0]), max(h_t)))
s_plot_zero <- data |>
  dplyr::filter(h_s == 0 & h_t %in% c(0, min(h_t[h_t > 0]), max(h_t)))


plot_s_de <- total_var - t_plot_pos$gamma[t_plot_pos$h_s == 0.2 & t_plot_pos$h_t == max(t_plot_pos$h_t)]
plot_s_ie <- t_plot_pos$gamma[t_plot_pos$h_s == 0.2 & t_plot_pos$h_t == max(t_plot_pos$h_t)] -
  t_plot_pos$gamma[t_plot_pos$h_s == 0 & t_plot_pos$h_t == max(t_plot_pos$h_t)]

ggplot(data  = t_plot_pos) +
  geom_line(mapping = aes(x = h_t, y = sigma, colour = as.factor(h_s)), size = 1.5) +
  geom_point(data = t_plot_zero, mapping = aes(x = h_t, y = sigma, colour = as.factor(h_s)), size = 4) +
  scale_colour_viridis_d(name = "Spatial Distance", begin = 0, end = 0.9) +
  theme_minimal() +
  labs(x = "Time Distance (years)",
       y = "Estimated Covariance")