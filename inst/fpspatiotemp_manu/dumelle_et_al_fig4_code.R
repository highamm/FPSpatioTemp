## code for Dumelle et al. (2021) Figure 4
## need to modify (and try to make it shorter, or wrap some pieces
## in a plotting function)

## next step fix the max spatial distance to match the moose problem

## second next step: wrap some of this in a plotting function

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
axistext_size <- 24
legendtext_size <- 24
annotate_size <- 9

s_epstol <- c(1, 4, 10, Inf)

# setting unique distance values
h_s_seq <- seq(0, s_rangetol, length.out = 200)
h_s_seq <- append(h_s_seq, s_epstol, s_range / 2)
h_t_seq <- seq(0, t_rangetol, length.out = 200)
h_t_seq <- append(h_t_seq, t_range / 2)
h_s <- rep(h_s_seq, times = length(h_t_seq))
h_t <- rep(h_t_seq, each = length(h_s_seq))

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

str(sigma)
# make semivariogram
gamma <- DumelleEtAl2021STLMM::make_stsemivariogram(covparams, h_s, h_t,
                                                    "exponential",
                                                    "exponential")

# make data frame
data <- data.frame(
  h_s = h_s,
  h_t = h_t,
  sigma = sigma,
  gamma = gamma
)


# temporal plotting subset
t_plot_pos <- data |>
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t > 0.2))

t_plot_zero <- data |>
  dplyr::filter((h_s %in% c(0, s_epstol)) & (h_t == 0))

plot_s_de <- total_var - t_plot_pos$gamma[t_plot_pos$h_s == 10 & t_plot_pos$h_t == max(t_plot_pos$h_t)]
plot_s_ie <- t_plot_pos$gamma[t_plot_pos$h_s == 10 & t_plot_pos$h_t == max(t_plot_pos$h_t)] -
  t_plot_pos$gamma[t_plot_pos$h_s == 0 & t_plot_pos$h_t == max(t_plot_pos$h_t)]

sv_tempx_ps_reml <- ggplot(t_plot_pos) +
  geom_line(mapping = aes(x = h_t, y = gamma, linetype = as.factor(h_s)), size = 1.5) +
  geom_point(t_plot_zero, mapping = aes(x = h_t, y = gamma), size = 4) +
  labs(x = "Temporal Distance $(h_t)$ in days", y = "Semivariance") +
  scale_linetype_discrete(name = "Spatial Distance $(h_s)$ in km", labels = c(bquote(h[s] == 0), bquote(h[s] == 0^"+"), bquote(h[s] == 750), bquote(h[s] == 1500), bquote(h[s] == infinity))) +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, "$\\infty$")) +
  scale_y_continuous(breaks = c(0, 30, 60, 90, 120), labels = c(0, 30, 60, 90, 120)) +
  theme(
    axis.title.x = element_text(face = "bold", size = axistext_size, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face = "bold", size = axistext_size, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = legendtext_size),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.5, 0.5),
    legend.text = element_text(size = legendtext_size),
    legend.key.size = unit(3, "line"),
    axis.text.x = element_text(size = axistext_size, face = "bold", colour = "black"),
    axis.text.y = element_text(size = axistext_size, face = "bold", colour = "black"),
    legend.key = element_rect(fill = NA),
    legend.key.height = unit(1, "cm")
  ) +
  expand_limits(x = c(-1 * rangetol, t_rangetol + 5.5 * rangetol), y = c(0, 120)) +
  # spatial independent
  # annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = total_var - (plot_s_de + plot_s_ie), yend =  total_var - (plot_s_de + plot_s_ie), size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace, xend = t_rangetol + brspace + brlen, y = total_var - plot_s_de, yend =  total_var - plot_s_de, size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + brlen, y = total_var - (plot_s_de + plot_s_ie), yend =  total_var - plot_s_de, size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace + brlen, xend = t_rangetol + brspace + 2 * brlen, y = total_var - (plot_s_de + plot_s_ie / 2), yend =  total_var - (plot_s_de + plot_s_ie / 2), size = 1.5) +
  # annotate("segment", x = t_rangetol + brspace + 2 * brlen, xend = t_rangetol + brspace + 2 * brlen, y = total_var - (plot_s_de + plot_s_ie / 2), yend =  total_var - (plot_s_de + plot_s_ie / 2) - 40 * brlen, size = 1.5) +
  # annotate("text", x = t_rangetol + brspace + 2.75 * brlen, y = total_var - (plot_s_de + plot_s_ie / 2) - 80 * brlen, label = TeX("$\\sigma^2_{\\gamma}$"), size = annotate_size, color = "black") +
  # spatial dependent
  annotate("segment", x = t_rangetol + 4 * brspace + 2 * brlen, xend = t_rangetol + 4 * brspace + brlen + 2 * brlen, y = total_var - plot_s_de, yend = total_var - plot_s_de, size = 1.5) +
  annotate("segment", x = t_rangetol + 4 * brspace + 2 * brlen, xend = t_rangetol + 4 * brspace + brlen + 2 * brlen, y = total_var - 0, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = t_rangetol + 4 * brspace + brlen + 2 * brlen, xend = t_rangetol + 4 * brspace + brlen + 2 * brlen, y = total_var - plot_s_de, yend = total_var - 0, size = 1.5) +
  annotate("segment", x = t_rangetol + 4 * brspace + brlen + 2 * brlen, xend = t_rangetol + 4 * brspace + 2.75 * brlen + 2 * brlen, y = total_var - plot_s_de / 2, yend = total_var - plot_s_de / 2, size = 1.5) +
  annotate("text", x = t_rangetol + brspace + 5 * brlen + 5 * brlen, y = total_var - plot_s_de / 2, label = "$\\sigma^2_{\\delta}$", size = annotate_size, color = "black")

sv_tempx_ps_reml
