xcoords <- 1:20; ycoords <- 1:20
allcoords <- expand.grid(xcoords, ycoords)

ntime <- 5
times <- 1:ntime

xs <- rep(allcoords$Var1, times = ntime)
ys <- rep(allcoords$Var2, times = ntime)
ts <- rep(times, each = nrow(allcoords))
coords_times <- data.frame(xcoords = xs, ycoords = ys, ts = ts)

spatvec <- c(9, 4)
rhotime <- 0.7

betavec <- 10
XDesign <- matrix(1, nrow = nrow(coords_times), ncol = 1)

res <- sim.spat.temp(coords_times, spatvec, rhotime, betavec, XDesign)
res$df


sim_obj <- sim.spat.temp(coords_times, spatvec, rhotime, betavec, XDesign)
nsamp <- 500
method <- "Random"

give_obj <- give.nas(sim_obj, nsamp, method)

fp.predict(give_obj)

