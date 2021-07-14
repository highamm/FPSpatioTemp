## toy example

## simulate 4 spatial sites (2 x 2 grid) and 3 time points
nx <- 2; ny <- 2; ntime <- 3
betavec <- 0

## set covariance parameters
sigma_parsil_spat <- 0.9
range <- 5
sigma_nugget_spat <- 0.1
sigma_parsil_time <- 0.7
rho <- 0.8
sigma_nugget_time <- 0.3
sigma_nugget_spacetime <- 0.4

library(FPSpatioTemp)

sim_toy <- sim_spatiotemp(nx = nx, ny = ny, ntime = ntime, betavec = betavec,
                          XDesign = matrix(1, nrow = nx * ny * ntime),
                          sigma_parsil_spat = sigma_parsil_spat,
                          range = range, 
                          sigma_nugget_spat = sigma_nugget_spat, 
                          sigma_parsil_time = sigma_parsil_time,
                          rho = rho,
                          sigma_nugget_time = sigma_nugget_time,
                          sigma_nugget_spacetime = sigma_nugget_spacetime,
                          seed = round(runif(1, min = 1, max = 10000000)))

toy_df <- sim_toy$out_df
plot_sp(xcoords = toy_df$xcoords, ycoords = toy_df$ycoords,
        times = toy_df$times, response = toy_df$response)
plot_t(xcoords = toy_df$xcoords, ycoords = toy_df$ycoords,
        times = toy_df$times, response = toy_df$response)

## note that observations are simulated so that they are ordered by
## time first and then by space within time
nspat <- nx * ny
ntime <- ntime
N <- nspat * ntime

## create Zt matrix
Zt <- lapply(1:ntime, matrix, data = 0, nrow = nspat, ncol = ntime)

for (i in 1:ntime) {
  Zt[[i]][ ,i] <- 1
}
Zt <- do.call(rbind, Zt)
Zt

library(dplyr)
## create Zs matrix
onetime <- diag(1, nspat) %>% as.data.frame()
Zs <- onetime %>% slice(rep(row_number(), ntime)) %>% as.matrix()


## create Rt matrix
times <- 1:ntime
H <- abs(outer(times, times, "-")) 
Rt <- rho ^ H
Rt

## create Rs matrix
allcoords <- sim_toy$out_df[ ,c("xcoords", "ycoords")]
uniquecoords <- unique(allcoords)
distancemat <- as.matrix(dist(uniquecoords))
Rs <- exp(-distancemat / range)



## examine each variance component: delta
p1 <- sigma_parsil_spat * Zs %*% Rs %*% t(Zs)
p1

## gamma: similar to a random effect for spatial site
## matrix value is 0 if observations are at different spatial locations,
## 0.1 if at the same spatial location.
## this would allow sites to have different means, even in the 
## absence of spatial correlation from sigma_parsil_spat
p2 <- sigma_nugget_spat * Zs %*% diag(1, nrow = nrow(uniquecoords)) %*% t(Zs)
p2

## tau
p3 <- sigma_parsil_time * Zt %*% Rt %*% t(Zt)
p3

## eta: similar to a random effect for time
## matrix value is 0 if observations are at different time points
## this would allow time points to have different means, even in the 
## absence of temporal correlation from Zt

p4 <- sigma_nugget_time * Zt %*% diag(1, ntime) %*% t(Zt)
p4

## independent error
p5 <- sigma_nugget_spacetime * diag(1, N)
p5

## add them all up to get the overall covariance matrix
p1 + p2 + p3 + p4 + p5


## can mess around with the parameters too: e.g. setting sigma_parsil_spat and sigma_nugget_spat to 0 should give a covariance matrix 
## where locations don't matter:

p1 <- 0 * Zs %*% Rs %*% t(Zs)
p2 <- 0 * Zs %*% diag(1, nrow = nrow(uniquecoords)) %*% t(Zs)
p3 <- sigma_parsil_time * Zt %*% Rt %*% t(Zt)
p4 <- sigma_nugget_time * Zt %*% diag(1, ntime) %*% t(Zt)
p5 <- sigma_nugget_spacetime * diag(1, N)
p1 + p2 + p3 + p4 + p5

## setting the spatial sigmas, and sigma_parsil_time all equal to 0
## should give you a covariance matrix as if you were fitting a 
## simple random effects model with "year" as the "site-level covariate"

p1 <- 0 * Zs %*% Rs %*% t(Zs)
p2 <- 0 * Zs %*% diag(1, nrow = nrow(uniquecoords)) %*% t(Zs)
p3 <- 0 * Zt %*% Rt %*% t(Zt)
p4 <- sigma_nugget_time * Zt %*% diag(1, ntime) %*% t(Zt)
p5 <- sigma_nugget_spacetime * diag(1, N)
p1 + p2 + p3 + p4 + p5
