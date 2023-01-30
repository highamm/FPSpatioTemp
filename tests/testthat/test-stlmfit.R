## library(testthat)
test_that("the model runs for exponential covariance", {
  expect_error(stlmfit(formula = response_na ~ x,
                       data = samp_data, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times"), NA)
})

test_that("the model fixed effects coefficients do not change from known correct results", {
  mod <- stlmfit(formula = response_na ~ x,
                       data = samp_data, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times")
  names(mod$fixed_parms) <- NULL
  expect_equal(mod$fixed_parms, c(4.698, -0.412), tolerance = 1e-3)
})

test_that("the model does not change from the current snapshot", {
  mod <- stlmfit(formula = response_na ~ x,
                 data = samp_data, xcoord = "xcoords",
                 ycoord = "ycoords", tcoord = "times")
  expect_snapshot(mod)
})

test_that("the model covariance parameter estimates do not change from known correct results", {
  mod <- stlmfit(formula = response_na ~ x,
                 data = samp_data, xcoord = "xcoords",
                 ycoord = "ycoords", tcoord = "times")
  names(mod$cov_parms) <- NULL
  expect_equal(mod$cov_parms, c(0.0062, 0.4383,
                                0.8174, 0.0004,
                                0.0000, 1.2882,
                                1.0005, 0.0899),
               tolerance = 1e-3)
})


test_that("the model runs for gaussian covariance", {
  expect_error(stlmfit(formula = response_na ~ x,
                       data = samp_data, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times",
                       cor_model_sp = "gaussian",
                       cor_model_t = "gaussian"), NA)
})

test_that("re-ordering data frame does not change parameter estimates", {
  samp_data_reordered <- samp_data[sample(1:nrow(samp_data), replace = FALSE), ]
  mod1 <- stlmfit(formula = response_na ~ x,
                 data = samp_data, xcoord = "xcoords",
                 ycoord = "ycoords", tcoord = "times")
  mod2 <- stlmfit(formula = response_na ~ x,
                  data = samp_data_reordered, xcoord = "xcoords",
                  ycoord = "ycoords", tcoord = "times")
  expect_equal(mod1$fixed_parms, mod2$fixed_parms)
  expect_equal(mod1$cov_parms, mod2$cov_parms)
})

test_that("model is fit when some covariates are missing", {
  expect_error(stlmfit(formula = response_na ~ x_miss,
                       data = samp_data, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times"), NA)
  expect_error(stlmfit(formula = response_na ~ x_fact_miss,
                       data = samp_data, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times"), NA)
})

test_that("model is fit when some spatial or temporal coordinates are missing", {
  expect_error(stlmfit(formula = response_na ~ x,
                       data = samp_data, xcoord = "xcoords_miss",
                       ycoord = "ycoords", tcoord = "times"), NA)
  expect_error(stlmfit(formula = response_na ~ x,
                       data = samp_data, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times_miss"), NA)
})

test_that("model is fit when there is implicitly missing space-time data", {
  samp_data_miss <- samp_data[-c(1, 4, 5, 120), ]
  expect_error(stlmfit(formula = response_na ~ x,
                       data = samp_data_miss, xcoord = "xcoords",
                       ycoord = "ycoords", tcoord = "times"), NA)
})


