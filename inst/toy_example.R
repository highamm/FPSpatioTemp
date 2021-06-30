set.seed(06302021)

response_year1 <- rnorm(100, 10, 5)
response_year2 <- response_year1 + rnorm(100, 0, 1)
response_year3 <- response_year2 + rnorm(100, 0, 1)
x <- 1:10; y <- 1:10
coords <- expand.grid(x, y)
allcoords <- rbind(coords, coords, coords)
times <- c(rep(1, 100), rep(2, 100), rep(3, 100))

sampindx <- sample(1:100, size = 30, replace = FALSE)
response_year1_na <- response_year1
response_year1_na[-sampindx] <- NA

sampindx2 <- sample(1:100, size = 25, replace = FALSE)
response_year2_na <- response_year2
response_year2_na[-sampindx2] <- NA

sampindx3 <- sample(1:100, size = 40, replace = FALSE)
response_year3_na <- response_year3
response_year3_na[-sampindx3] <- NA

response_all <- c(response_year1, response_year2, response_year3)
response <- c(response_year1_na, response_year2_na, response_year3_na)

df <- data.frame(response = response, response_all = response_all,
                 xcoords = allcoords[ ,1], ycoords = allcoords[ ,2],
                 times = times)

library(dplyr)
df <- df %>% mutate(predind = if_else(times == 3, true = 1, false = 0))

fit_toy <- stlmfit(formula = response ~ 1, data = df,
                   xcoordcol = "xcoords",
                   ycoordcol = "ycoords",
                   tcol = "times", wtscol = "predind")
fit_toy[[1]]
sqrt(fit_toy[[2]])
df %>% filter(times == 3) %>% summarise(total = sum(response_all))
## seems reasonable

fit_toy[[7]]
## these also seem reasonable: we see that there is a decent amount of 
## temporal correlation and little spatial correlation. Data were
## simulated to be correlated temporally, but not spatially, so these
## estimates seem okay.
