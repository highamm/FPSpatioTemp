library(FPSpatioTemp)
library(tidyverse)
library(sptotal)

seed <- sample(1e9, size = 1)

sim_fun <- function(resp_type = "normal", 
                    nx = 10, ny = 10, ntime = 5, betavec = 0,
                    sp_de = 0.9, sp_range = sqrt(2) / 3,
                    sp_ie = 0.1, t_de = 0.7,
                    t_range = 1 / 3, t_ie = 0.3,
                    spt_de = 0.9, spt_ie = 0.4, n = 100,
                    pred_level = 0.90) {
  
  sim_start <- Sys.time()
  
  sim_obj <- sim_spatiotemp(resp_type = resp_type,
                            nx = nx, ny = ny, ntime = ntime,
                            betavec = betavec,
                            XDesign = matrix(1, nrow = nx * ny * ntime),
                            sp_de = sp_de,
                            sp_range = sp_range, 
                            sp_ie = sp_ie, 
                            t_de = t_de,
                            t_range = t_range,
                            t_ie = t_ie,
                            spt_de = spt_de,
                            spt_ie = spt_ie,
                            seed = seed)
  
  sim_end <- Sys.time()
  sim_time <- sim_end - sim_start
  
  samp_obj <- sample_spatiotemp(obj = sim_obj, n = n,
                                samp_type = "random")

  samp_df <- samp_obj$df_full
  truetotal <- samp_df |> dplyr::filter(times == max(times)) |>
    dplyr::summarise(truetotal = sum(response))
  
  samp_df <- samp_df |>
    dplyr::mutate(wts = dplyr::if_else(times == max(times),
                                       true = 1, false = 0))
  
  samp_end <- Sys.time()
  samp_time <- samp_end - sim_end
    
  fit <- stlmfit(formula = response_na ~ 1, data = samp_df,
                 xcoord = "xcoords", ycoord = "ycoords",
                 tcoord = "times")

  fit_end <- Sys.time()
  fit_time <- fit_end - samp_end
  
  pred_fit <- predict(fit, wts = "wts", pred_level = pred_level)
 
  pred_end <- Sys.time()
  pred_time <- pred_end - fit_end
  
  cov_df <- fit$cov_parms |> as.matrix() |> t() |> as.data.frame()
  fixed_df <- fit$fixed_parms |> as.matrix() |> t() |> as.data.frame()
  sim_parm_df <- data.frame(resp_type_sim = resp_type,
                            betavec_sim = betavec,
                            sp_de_sim = sp_de, sp_range_sim = sp_range,
                            sp_ie_sim = sp_ie, t_de_sim = t_de,
                            t_range_sim = t_range, t_ie_sim = t_ie,
                            spt_ie_sim = spt_ie, spt_de_sim = spt_de, n_sim = n)
  
  
  ## sptotal comparison
  oneyear_df <- samp_df |> dplyr::filter(times == max(times))
  sp_mod <- slmfit(response_na ~ 1, data = oneyear_df, xcoordcol = "xcoords",
         ycoordcol = "ycoords")
  sp_pred <- predict(sp_mod)
  
  
  srs_df <- oneyear_df |> summarise(samp_mean = mean(response_na, na.rm = TRUE),
                          samp_var = var(response_na, na.rm = TRUE),
                          N = nrow(oneyear_df),
                          n = sum(!is.na(response_na)))
  pred_srs <- srs_df$samp_mean * srs_df$N
  se_srs <- sqrt(srs_df$N ^ 2 * srs_df$samp_var / srs_df$n * 
                   (1 - srs_df$n / srs_df$N))
  lb_srs <- pred_srs - -qnorm((1 - pred_level) / 2) * se_srs
  ub_srs <- pred_srs - qnorm((1 - pred_level) / 2) * se_srs
  
  sim_output <- data.frame(pred = pred_fit$totalpred,
                           se = sqrt(pred_fit$predvar),
                           lb = pred_fit$lb, ub = pred_fit$ub,
                           truetotal = as.numeric(truetotal),
                           cov_df,
                           fixed_df,
                           sim_parm_df,
                           sim_time = sim_time,
                           samp_time = samp_time,
                           fit_time = fit_time,
                           pred_time = pred_time,
                           pred_sptot = sp_pred$FPBK_Prediction |> as.vector(),
                           se_sptot = sp_pred$PredVar |> sqrt() |> as.vector(),
                           lb_sptot = sp_pred$conf_bounds[1],
                           ub_sptot = sp_pred$conf_bounds[2],
                           pred_srs = pred_srs,
                           se_srs = se_srs,
                           lb_srs = lb_srs,
                           ub_srs = ub_srs)
  
  return(sim_output)
}      

## do 1 simulation
## old sim parms: nx = 20, ny = 20, ntime = 10, n = 500
# 
# sim_output_df <- sim_fun(resp_type = "lognormal",
#                          nx = 10, ny = 10, ntime = 10, betavec = 10,
#                          sp_de = 0, sp_ie = 0, t_de = 0, t_ie = 0,
#                          spt_de = 1.5, spt_ie = 0.5,
#                          sp_range = sqrt(2) / 3, t_range = 1e-10, n = 100,
#                          pred_level = 0.90)


const <- 2.886
sim_parm_df <- tibble::tribble(~sp_de, ~sp_ie, ~t_de, ~t_ie, ~spt_de, ~spt_ie,
                ~sp_range, ~t_range, ~n, ~resp_type,
                0.5, 0.1675, 0.5, 0.1675, 0.5, 0.1675,
                sqrt(2) / 3, 1 / 3, 250, "normal",
                0.5, 0.1675, 0.5, 0.1675, 0.5, 0.1675,
                sqrt(2) / 3, 1 / 3, 500, "normal",
                0, 0, 0, 0, 0, 2,
                1e-10, 1e-10, 250, "normal",
                0, 0, 0, 0, 0, 2,
                1e-10, 1e-10, 500, "normal",
                0, 0, 0, 0, 1.5, 0.5,
                sqrt(2) / 3, 1e-10, 250, "normal",
                0, 0, 0, 0, 1.5, 0.5,
                sqrt(2) / 3, 1e-10, 500, "normal",
                0, 0, 0, 1.5, 0.25, 0.25,
                sqrt(2) / 3, 1e-10, 250, "normal",
                0, 0, 0, 1.5, 0.25, 0.25,
                sqrt(2) / 3, 1e-10, 500, "normal",
                
                ## lognormal
                0.5 / const, 0.1675 / const, 0.5 / const, 0.1675 / const,
                0.5 / const, 0.1675 / const, 
                sqrt(2) / 3, 1 / 3, 250, "lognormal",
                0.5 / const, 0.1675 / const, 0.5 / const, 0.1675 / const,
                0.5 / const, 0.1675 / const,
                sqrt(2) / 3, 1 / 3, 500, "lognormal",
                0, 0, 0, 0, 0, 2 / const,
                1e-10, 1e-10, 250, "lognormal",
                0, 0, 0, 0, 0, 2 / const,
                1e-10, 1e-10, 500, "lognormal",
                0, 0, 0, 0, 1.5 / const, 0.5 / const,
                sqrt(2) / 3, 1e-10, 250, "lognormal",
                0, 0, 0, 0, 1.5 / const, 0.5 / const,
                sqrt(2) / 3, 1e-10, 500, "lognormal",
                0, 0, 0, 1.5 / const, 0.25 / const, 0.25 / const,
                sqrt(2) / 3, 1e-10, 250, "lognormal",
                0, 0, 0, 1.5 / const, 0.25 / const, 0.25 / const,
                sqrt(2) / 3, 1e-10, 500, "lognormal")
# 
# sim_parm_df <- tibble::tribble(~sp_de, ~sp_ie, ~t_de, ~t_ie, ~spt_de, ~spt_ie,
#                                ~sp_range, ~t_range, ~n, ~resp_type,
#                                0.5, 0.1675, 0.5, 0.1675, 0.5, 0.1675,
#                                sqrt(2) / 3, 1 / 3, 100, "normal",
#                                0.5, 0.1675, 0.5, 0.1675, 0.5, 0.1675,
#                                sqrt(2) / 3, 1 / 3, 100, "lognormal")
#    

sim_output_df <- pmap_dfr(sim_parm_df, sim_fun,
                          nx = 10, ny = 10, ntime = 10,
                          betavec = 0, pred_level = 0.90) 



sim_output_df <- sim_output_df |>
  mutate(conf_ind = if_else(truetotal >= lb & truetotal <= ub,
                            true = 1, false = 0),
         conf_ind_sptot = if_else(truetotal >= lb_sptot & truetotal <= ub_sptot,
                                  true = 1, false = 0),
         conf_ind_srs = if_else(truetotal >= lb_srs & truetotal <= ub_srs,
                 true = 1, false = 0),
         seed = seed)

# sim_output_df |> summarise(coverage = mean(conf_ind),
#                             rmspe = sqrt(sum((pred - truetotal) ^ 2)),
#                             medci = median(ub - lb),
#                             meanse = mean(se))

## sim settings: 
library(tidyverse)
write_csv(sim_output_df, file = paste("sims/sim_output", seed, ".csv",
                                      sep = ""),
          col_names = TRUE, 
          append = TRUE)



# sims_df <- read_csv("inst/simulations/sim_output.csv")
# 
# sims_df |> summarise(coverage = mean(conf_ind),
#                      rmspe = sqrt(sum((pred - truetotal) ^ 2)),
#                      medci = median(ub - lb),
#                      meanse = mean(se))
