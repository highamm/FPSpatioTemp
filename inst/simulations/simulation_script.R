library(FPSpatioTemp)
library(dplyr)

## nrows and ncols for spatial grid
nx <- 10; ny <- 10;

## number of time points
ntime <- 5

## overall mean
betavec <- 10

## set covariance parameters
sigma_parsil_spat <- 0.9
range <- 5
sigma_nugget_spat <- 0.1
sigma_parsil_time <- 0.7
rho <- 0.8
sigma_nugget_time <- 0.3
sigma_nugget_spacetime <- 0.4

samp_size <- 100

sim_fun <- function() {
  sim_obj <- sim_spatiotemp(nx = nx, ny = ny, ntime = ntime, betavec = betavec,
                            XDesign = matrix(1, nrow = nx * ny * ntime),
                            sigma_parsil_spat = sigma_parsil_spat,
                            range = range, 
                            sigma_nugget_spat = sigma_nugget_spat, 
                            sigma_parsil_time = sigma_parsil_time,
                            rho = rho,
                            sigma_nugget_time = sigma_nugget_time,
                            sigma_nugget_spacetime = sigma_nugget_spacetime,
                            seed = round(runif(1, min = 1, max = 10000000)))
  
  samp_obj <- sample_spatiotemp(sim_obj$out_df, n = samp_size,
                                samp_type = "random", seed = sim_obj$seed)
  
  samp_df <- samp_obj$df_full
  truetotal <- samp_df %>% filter(times == max(times)) %>%
    summarise(truetotal = sum(response))
  
  samp_df <- samp_df %>% mutate(wts = if_else(times == max(times),
                                              true = 1, false = 0))
  
  fit <- stlmfit(formula = response_na ~ 1, data = samp_df,
                 xcoordcol = "xcoords", ycoordcol = "ycoords", tcol = "times",
                 wtscol = "wts")
  
  sim_output <- data.frame(pred = fit$totalpred, se = sqrt(fit$predvar),
                           lb = fit$lb, ub = fit$ub,
                           truetotal = as.numeric(truetotal),
                           parms = fit$parms)
  
  return(sim_output)
}            

## do 1 simulation
sim_1 <- sim_fun()

sim_output <- replicate(100, sim_fun(), simplify = "matrix")                 
sim_output_df <- t(sim_output) %>% as_tibble()  
sim_output_df <- tidyr::unnest(sim_output_df, cols = c(pred, se, lb, ub, truetotal, parms.sigma_parsil_spat, parms.range, 
                                      parms.sigma_nugget_spat, parms.sigma_parsil_time, parms.rho, 
                                      parms.sigma_nugget_time, parms.sigma_nugget_spacetime))

sim_output_df <- sim_output_df %>% mutate(conf_ind = if_else(truetotal >= lb & truetotal <= ub, true = 1, false = 0))
sim_output_df %>% summarise(coverage = mean(conf_ind))
sim_output_df %>% select(lb, ub, truetotal, conf_ind) %>% print(n = Inf)
