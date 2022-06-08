## analysis for estimating parameters for moose data set

## UNADDRESSED ISSUE: What to do with the sites from 2019 that are outside of
## the sampling frame for all of the other years.
## for now, these sites are removed.

library(DumelleEtAl2021STLMM)
library(FPSpatioTemp)
library(dplyr)
library(tidyverse)

data(moose_14_20)
moose_14_20
moose_14_20 %>% group_by(Surveyyear) %>% count()
moose_14_20 %>% group_by(Surveyyear) %>% summarise(nsampframe = sum(samp_frame))

## drop the extra observations from 2019

moose_20 <- moose_14_20 %>% dplyr::filter(Surveyyear == 2020 & samp_frame == 1)
moose_14_20 <- semi_join(moose_14_20, moose_20, by = c("ID"))
moose_14_20 %>% group_by(Surveyyear) %>% count()


# totalvar <- var(moose_14_20$totalmoosena, na.rm = TRUE)

# s_de_initial <- totalvar / 10
# s_ie_initial <- totalvar / 10
# t_de_initial <- totalvar / 10
# t_ie_initial <- totalvar / 10
# st_de_initial <- totalvar / 10
# st_ie_initial <- totalvar / 10
# total_var_initial <- sum(
#   s_de_initial, s_ie_initial, t_de_initial,
#   t_ie_initial, st_de_initial, st_ie_initial
# )
# s_range_initial <- 50
# t_range_initial <- 3
# 
# swe_scale <- total_var_initial / (total_var_initial - st_de_initial)
# 
# siminitial <-  make_covparam_object(
#   s_de = s_de_initial, s_ie = s_ie_initial,
#   t_ie = t_ie_initial, t_de = t_de_initial,
#   st_de = st_de_initial, st_ie = st_ie_initial,
#   s_range = s_range_initial, t_range = t_range_initial,
#   stcov = "productsum", estmethod = "reml"
# )

moose_14_20 <- moose_14_20 %>% mutate(Stratname = str_to_upper(Stratname))


moose_14_20 <- moose_14_20 %>% mutate(strat2020 = if_else(Surveyyear == 2020,
                                                          true = Stratname,
                                                          false = NA_character_)) %>%
  arrange(xcoords, ycoords, Surveyyear) %>%
  fill(strat2020, .direction = "up")

moose_14_20 <- moose_14_20 %>%
  mutate(pred_ind = if_else(samp_frame == 1 & Surveyyear == 2020,
                                          true = 1, 
                                          false = 0))

moose_high <- moose_14_20 %>% dplyr::filter(strat2020 == "HIGH")
moose_high <- moose_high %>% mutate(Surveyyear2 = Surveyyear - min(Surveyyear) + 1)
moose_sampled_high <- moose_high %>% dplyr::filter(!is.na(totalmoosena))

# s_cor <- "exponential"
# t_cor <- "exponential"

# formula <- totalmoosena ~ 1
# data <- moose_high
# xcoordcol <- "xcoords"
# ycoordcol <- "ycoords"
# tcol <- "Surveyyear2"
# areacol <- NULL
# CorModel <- "exponential"
# initial <- siminitial

stlm_obj <- stlmfit(formula = totalmoosena ~ 1, data = moose_high,
                    xcoordcol = "xcoords", ycoordcol = "ycoords",
                    tcol = "Surveyyear2",
                    areacol = NULL,
                    CorModel = "exponential")

pred_obj <- predict(object = stlm_obj, wtscol = "pred_ind",
                            pred_level = 0.90)
pred_obj

# ps_reml_mod_high <- stlmm(
#   formula = totalmoosena ~ 1,
#   data = moose_sampled_high,
#   xcoord = "xcoords",
#   ycoord = "ycoords",
#   tcoord = "Surveyyear2",
#   stcov = "productsum",
#   estmethod = "reml",
#   s_cor = s_cor,
#   t_cor = t_cor,
#   initial = siminitial,
#   condition = 1e-4
# )
# 
# beta0_high <- as.vector(unname(ps_reml_mod_high$Coefficients))
# se_high <- sqrt(diag(ps_reml_mod_high$CovCoefficients))
# cov_parms_high <- unclass(ps_reml_mod_high$CovarianceParameters)


moose_low <- moose_14_20 %>% dplyr::filter(strat2020 == "LOW")
moose_low <- moose_low %>% mutate(Surveyyear2 = Surveyyear - min(Surveyyear) + 1)
moose_sampled_low <- moose_low %>% dplyr::filter(!is.na(totalmoosena))
moose_low %>% group_by(Surveyyear) %>% count()

stlm_obj_low <- stlmfit(formula = totalmoosena ~ 1, data = moose_low,
                        xcoordcol = "xcoords", ycoordcol = "ycoords",
                        tcol = "Surveyyear2",
                        areacol = NULL,
                        CorModel = "exponential")

pred_obj_low <- predict(object = stlm_obj_low, wtscol = "pred_ind",
                                pred_level = 0.90)
pred_obj_low

pred <- as.vector(pred_obj$totalpred + pred_obj_low$totalpred)
se <- as.vector(sqrt(pred_obj$predvar + pred_obj_low$predvar))
pred
pred + c(-1, 1) * 1.96 * se
## more narrow by 336 moose

# ps_reml_mod_low <- stlmm(
#   formula = totalmoosena ~ 1,
#   data = moose_sampled_low,
#   xcoord = "xcoords",
#   ycoord = "ycoords",
#   tcoord = "Surveyyear2",
#   stcov = "productsum",
#   estmethod = "reml",
#   s_cor = s_cor,
#   t_cor = t_cor,
#   initial = siminitial,
#   condition = 1e-4
# )
# 
# beta0_low <- as.vector(unname(ps_reml_mod_low$Coefficients))
# se_low <- sqrt(diag(ps_reml_mod_low$CovCoefficients))
# cov_parms_low <- unclass(ps_reml_mod_low$CovarianceParameters)
