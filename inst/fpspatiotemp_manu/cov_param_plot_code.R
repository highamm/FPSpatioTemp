## some code from Dumelle et al. 2021
## next steps: colour points appropriately (and just get rid of linetypes)

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
t_rangetol <- t_range * 2


brspace <- 0.075
brlen <- 0.15

# epsilon tolerance
epstol <- 1e-5
# axistext_size <- 24
# legendtext_size <- 24
# annotate_size <- 9

s_epstol <- c(0.2, Inf) ##  4, 10,


# setting unique distance values
h_s_seq <- seq(0, s_rangetol, length.out = 200)
h_s_seq <- append(h_s_seq, s_epstol, s_range / 2)
h_t_seq <- seq(0, t_rangetol, length.out = 200)
h_t_seq <- append(h_t_seq, t_range / 2)
h_s <- rep(h_s_seq, times = length(h_t_seq))
h_t <- rep(h_t_seq, each = length(h_s_seq))

t_epstol <- c(0.2, max(h_t))

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
t_plot_pos <- data %>%
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t > 0.2))
t_plot_zero <- data %>%
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t == 0))

s_plot_pos <- data %>% 
  dplyr::filter(h_s > 0 & h_t %in% c(0, min(h_t[h_t > 0]), max(h_t)))
s_plot_zero <- data %>%
  dplyr::filter(h_s == 0 & h_t %in% c(0, min(h_t[h_t > 0]), max(h_t)))


plot_s_de <- total_var - t_plot_pos$gamma[t_plot_pos$h_s == 0.2 & t_plot_pos$h_t == max(t_plot_pos$h_t)]
plot_s_ie <- t_plot_pos$gamma[t_plot_pos$h_s == 0.2 & t_plot_pos$h_t == max(t_plot_pos$h_t)] -
  t_plot_pos$gamma[t_plot_pos$h_s == 0 & t_plot_pos$h_t == max(t_plot_pos$h_t)]

ggplot(data  = t_plot_pos) +
  geom_line(mapping = aes(x = h_t, y = sigma, colour = as.factor(h_s)), size = 1.5) +
  geom_point(data = t_plot_zero, mapping = aes(x = h_t, y = sigma, colour = as.factor(h_s)), size = 4) +
  scale_colour_viridis_d(begin = 0, end = 0.9) +
  theme_minimal()


ggplot(data  = s_plot_pos) +
  geom_line(mapping = aes(x = h_s, y = sigma, colour = as.factor(h_t)), size = 1.5) +
  geom_point(data = s_plot_zero, mapping = aes(x = h_s, y = sigma, colour = as.factor(h_t)), size = 4) +
  scale_colour_viridis_d(end = 0.9) +
  theme_minimal()
