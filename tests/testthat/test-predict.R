test_that("prediction is equal to known correct value", {
  mod <- stlmfit(formula = response_na ~ x,
          data = samp_data, xcoord = "xcoords",
          ycoord = "ycoords", tcoord = "times")
  pred_obj <- predict(mod, "predwts")
  expect_equal(pred_obj[[1]], 147.9441, tolerance = 1e-3)
  expect_equal(pred_obj[[2]], 14.99772, tolerance = 1e-3)
  expect_equal(pred_obj[[3]], 141.5741, tolerance = 1e-3)
  expect_equal(pred_obj[[4]], 154.3141, tolerance = 1e-3)
})

test_that("default prediction is same as specifying weights of 1 to most recent year", {
  mod <- stlmfit(formula = response_na ~ x,
                 data = samp_data, xcoord = "xcoords",
                 ycoord = "ycoords", tcoord = "times")
  pred_obj <- predict(mod, "predwts")
  pred_obj2 <- predict(mod)
  expect_equal(pred_obj, pred_obj2)
})

test_that("prediction works for non-zero weights other than 0's and 1's", {
  samp_data$newweights <- samp_data$predwts * 4.2
  mod <- stlmfit(formula = response_na ~ x,
                 data = samp_data, xcoord = "xcoords",
                 ycoord = "ycoords", tcoord = "times")
  pred_obj <- predict(mod, "predwts")
  pred_obj2 <- predict(mod, "newweights")
  expect_equal(4.2 * pred_obj[[1]], pred_obj2[[1]], tolerance = 1e-3)
})

