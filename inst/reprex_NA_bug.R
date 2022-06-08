## reprex for NA

library(DumelleEtAl2021STLMM)
library(tidyverse)

df <- tibble(year = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4),
       xcoord = c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2),
       ycoord = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2),
       resp = c(10, 9, 8, 10, 0, 4, 2, 11, NA, NA, NA, NA, 5, 6, 10, 4))

totalvar <- var(df$resp, na.rm = TRUE)

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


s_range_initial <- 2
t_range_initial <- 2

swe_scale <- total_var_initial / (total_var_initial - st_de_initial)

siminitial <-  make_covparam_object(
  s_de = s_de_initial, s_ie = s_ie_initial,
  t_ie = t_ie_initial, t_de = t_de_initial,
  st_de = st_de_initial, st_ie = st_ie_initial,
  s_range = s_range_initial, t_range = t_range_initial,
  stcov = "productsum", estmethod = "reml"
)

## errors: non-conformable arguments
stlmm(
  formula = resp ~ 1,
  data = df,
  xcoord = "xcoord",
  ycoord = "ycoord",
  tcoord = "year",
  stcov = "productsum",
  estmethod = "reml",
  s_cor = "exponential",
  t_cor = "exponential",
  initial = siminitial,
  condition = 1e-4
)

## does not error
stlmm(
  formula = resp ~ 1,
  data = df %>% dplyr::filter(!is.na(resp)),
  xcoord = "xcoord",
  ycoord = "ycoord",
  tcoord = "year",
  stcov = "productsum",
  estmethod = "reml",
  s_cor = "exponential",
  t_cor = "exponential",
  initial = siminitial,
  condition = 1e-4
)
