## dumelle et al. analysis conversion

library(DumelleEtAl2021STLMM)

library(dplyr)

# load the data
data("or_data")

# subset into training and test
or_train <- subset(or_data, TYPE == "TRAIN")
or_test <- subset(or_data, TYPE == "TEST")

s_cor <- "exponential"
t_cor <- "exponential"

s_de_initial <- 40
s_ie_initial <- 2
t_de_initial <- 2
t_ie_initial <- 2
st_de_initial <- 7
st_ie_initial <- 7
total_var_initial <- sum(
  s_de_initial, s_ie_initial, t_de_initial,
  t_ie_initial, st_de_initial, st_ie_initial
)
s_range_initial <- 1000
t_range_initial <- 4
swe_scale <- total_var_initial / (total_var_initial - st_de_initial)

siminitial <-  make_covparam_object(
    s_de = s_de_initial, s_ie = s_ie_initial,
    t_ie = t_ie_initial, t_de = t_de_initial,
    st_de = st_de_initial, st_ie = st_ie_initial,
    s_range = s_range_initial, t_range = t_range_initial,
    stcov = "productsum", estmethod = "reml"
)
  
ps_reml_mod <- stlmm(
  formula = TMAX ~ ELEVATION + TIMES + PRCP,
  data = or_train,
  xcoord = "LONGITUDE",
  ycoord = "LATITUDE",
  tcoord = "TIMES",
  stcov = "productsum",
  estmethod = "reml",
  s_cor = s_cor,
  t_cor = t_cor,
  initial = siminitial,
  condition = 0
)

data.frame(
  stcov = "productsum",
  estmethod = "reml",
  beta = c("beta0", "beta1", "beta2", "beta3"),
  est = unname(ps_reml_mod$Coefficients),
  se = sqrt(diag(ps_reml_mod$CovCoefficients))
)

data.frame(
  stcov = "productsum",
  estmethod = "reml",
  covparam = names(unclass(ps_reml_mod$CovarianceParameters)),
  value = unclass(ps_reml_mod$CovarianceParameters)
)
