library(tidyverse)
library(parallel)
library(FPSpatioTemp)

sim_fun <- function(nx = 10, ny = 10, ntime = 5, betavec = 10,
                    sigma_parsil_spat = 0.9, range = 5,
                    sigma_nugget_spat = 0.1, sigma_parsil_time = 0.7,
                    rho = 0.8, sigma_nugget_time = 0.3,
                    sigma_nugget_spacetime = 0.4, n = 100,
                    pred_level = 0.90) {
  
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
                            seed = round(runif(1, min = 1, max = 10000000)))
  
  samp_obj <- sample_spatiotemp(obj = sim_obj, n = n,
                                samp_type = "random")
  
  samp_df <- samp_obj$df_full
  truetotal <- samp_df |> dplyr::filter(times == max(times)) |>
    dplyr::summarise(truetotal = sum(response))
  
  samp_df <- samp_df |>
    dplyr::mutate(wts = dplyr::if_else(times == max(times),
                                       true = 1, false = 0))
  
  fit <- stlmfit(formula = response_na ~ 1, data = samp_df,
                 xcoordcol = "xcoords", ycoordcol = "ycoords",
                 tcol = "times")
  pred_fit <- predict(fit, wtscol = "wts", pred_level = pred_level)
  
  
  
  sim_output <- data.frame(pred = as.vector(pred_fit$totalpred),
                           se = as.vector(sqrt(pred_fit$predvar)),
                           lb = pred_fit$lb, ub = pred_fit$ub,
                           truetotal = as.numeric(truetotal),
                           parms = fit$parms)
  
  return(sim_output)
  
}

## do 1 simulation
sim_1 <- sim_fun()

nrep <- 2

seeds <- sample(1e9, size = nrep)

n_cluster <- detectCores() # find cores (48 on mine)
cluster <- makeCluster(n_cluster) # make cluster

clusterEvalQ(cluster, library(FPSpatioTemp)) # export DvMsp to cluster
clusterEvalQ(cluster, library(tidyverse))

parLapply(
  cluster, # the cluster
  seeds,
  sim_fun
)

stopCluster(cluster)
