

test_that("correlation matrices can be built from spatial coordinates", {
  dist_mat <- stats::dist(samp_data[ ,c("xcoords", "ycoords")], diag = TRUE, upper = TRUE) |> as.matrix()
  expect_error(build_r(cov_type = "exponential", range = 4, dist_mat = dist_mat), NA)
  expect_error(build_r(cov_type = "gaussian", range = 4, dist_mat = dist_mat), NA)
  expect_error(build_r(cov_type = "tent", range = 4, dist_mat = dist_mat), NA)
  expect_error(build_r(cov_type = "spherical", range = 4, dist_mat = dist_mat), NA)
})
