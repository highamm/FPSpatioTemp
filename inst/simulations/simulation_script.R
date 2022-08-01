library(FPSpatioTemp)
library(tidyverse)


seed <- sample(1e9, size = 1)

sim_fun <- function(nx = 10, ny = 10, ntime = 5, betavec = 10,
                    sp_de = 0.9, sp_range = 5,
                    sp_ie = 0.1, t_de = 0.7,
                    t_range = 3, t_ie = 0.3,
                    spt_ie = 0.4, spt_de = 0.9, n = 100,
                    pred_level = 0.90) {
  
  sim_start <- Sys.time()
  
  sim_obj <- sim_spatiotemp(nx = nx, ny = ny, ntime = ntime,
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
  
  pred_fit <- predict(fit, wtscol = "wts", pred_level = pred_level)
 
  pred_end <- Sys.time()
  pred_time <- pred_end - fit_end
  
  cov_df <- fit$cov_parms |> as.matrix() |> t() |> as.data.frame()
  fixed_df <- fit$fixed_parms |> as.matrix() |> t() |> as.data.frame()
  
  sim_output <- data.frame(pred = pred_fit$totalpred,
                           se = sqrt(pred_fit$predvar),
                           lb = pred_fit$lb, ub = pred_fit$ub,
                           truetotal = as.numeric(truetotal),
                           cov_df,
                           fixed_df,
                           sim_time = sim_time,
                           samp_time = samp_time,
                           fit_time = fit_time,
                           pred_time = pred_time)
  
  return(sim_output)
}      

## do 1 simulation
## old sim parms: nx = 20, ny = 20, ntime = 10, n = 500
sim_output_df <- sim_fun(nx = 6, ny = 6, ntime = 3, n = 30)

sim_output_df <- sim_output_df |>
  mutate(conf_ind = if_else(truetotal >= lb & truetotal <= ub,
                            true = 1, false = 0),
         seed = seed)

# sim_output_df |> summarise(coverage = mean(conf_ind),
#                             rmspe = sqrt(sum((pred - truetotal) ^ 2)),
#                             medci = median(ub - lb),
#                             meanse = mean(se))

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
