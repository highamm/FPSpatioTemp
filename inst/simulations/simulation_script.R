library(FPSpatioTemp)
library(tidyverse)

seed <- sample(1e9, size = 1)

sim_fun <- function(nx = 10, ny = 10, ntime = 5, betavec = 10,
                    sigma_parsil_spat = 0.9, range = 5,
                    sigma_nugget_spat = 0.1, sigma_parsil_time = 0.7,
                    rho = 3, sigma_nugget_time = 0.3,
                    sigma_nugget_spacetime = 0.4, n = 100,
                    pred_level = 0.90) {
  
  nx <- 20; ny <- 20; ntime <- 10
  n <- 500
  
  sim_start <- Sys.time()
  
  sim_obj <- sim_spatiotemp(nx = nx, ny = ny, ntime = ntime,
                            betavec = betavec,
                            XDesign = matrix(1, nrow = nx * ny * ntime),
                            sigma_parsil_spat = sigma_parsil_spat,
                            range = range, 
                            sigma_nugget_spat = sigma_nugget_spat, 
                            sigma_parsil_time = sigma_parsil_time,
                            rho = rho,
                            sigma_nugget_time = sigma_nugget_time,
                            sigma_nugget_spacetime = sigma_nugget_spacetime,
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
                 xcoordcol = "xcoords", ycoordcol = "ycoords",
                 tcol = "times")

  fit_end <- Sys.time()
  fit_time <- fit_end - samp_end
  
  pred_fit <- predict(fit, wtscol = "wts", pred_level = pred_level)
 
  pred_end <- Sys.time
  pred_time <- pred_end - fit_end
  
  sim_output <- data.frame(pred = as.vector(pred_fit$totalpred),
                           se = as.vector(sqrt(pred_fit$predvar)),
                           lb = pred_fit$lb, ub = pred_fit$ub,
                           truetotal = as.numeric(truetotal),
                           parms = fit$parms,
                           sim_time = sim_time,
                           samp_time = samp_time,
                           fit_time = fit_time,
                           pred_time = pred_time)
  
  return(sim_output)
}          

## do 1 simulation
sim_1 <- sim_fun()


## don't really need this anymore
nrep <- 1
sim_output <- replicate(nrep, sim_fun(), simplify = "matrix")                 
sim_output_df <- t(sim_output) |> as_tibble()  

sim_output_df <- tidyr::unnest(sim_output_df,
                               cols = c(pred, se, lb, ub, truetotal,
                                        parms.yo,
                                        parms.sigma_parsil_spat,
                                        parms.range,
                                        parms.sigma_nugget_spat,
                                        parms.sigma_parsil_time,
                                        parms.rho,
                                        parms.sigma_nugget_time,
                                        parms.sigma_nugget_spacetime,
                                        parms.sigma_parsil_spacetime))

sim_output_df <- sim_output_df |>
  mutate(conf_ind = if_else(truetotal >= lb & truetotal <= ub,
                            true = 1, false = 0),
         seed = seed)

sim_output_df |> summarise(coverage = mean(conf_ind),
                            rmspe = sqrt(sum((pred - truetotal) ^ 2)),
                            medci = median(ub - lb),
                            meanse = mean(se))

library(tidyverse)
# write_csv(sim_output_df, file = paste("sim_output", a, ".csv", sep = ""))
write_csv(sim_output_df, file = paste("sims/sim_output", seed, ".csv",
                                      sep = ""),
          col_names = TRUE, 
          append = TRUE)
# }



# sims_df <- read_csv("inst/simulations/sim_output.csv")
# 
# sims_df |> summarise(coverage = mean(conf_ind),
#                      rmspe = sqrt(sum((pred - truetotal) ^ 2)),
#                      medci = median(ub - lb),
#                      meanse = mean(se))
