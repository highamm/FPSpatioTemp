test_that("z matrices are full", {
  toy_df <- tibble::tibble(`_spindex` = c(1, 2, 1),
                       `_tindex` = c(1, 1, 2),
                       .observed = c(TRUE, TRUE, TRUE))
  expect_equal(nrow(build_z_mats(toy_df)[[1]]), 4)
})

test_that("z matrices drop irrelevant rows", {
  toy_df <- tibble::tibble(`_spindex` = c(1, 2, 1, 2),
                           `_tindex` = c(1, 1, 2, 2),
                       .observed = c(TRUE, FALSE, TRUE, TRUE))
  expect_equal(nrow(build_z_mats(toy_df)[[1]]), 3)
})


