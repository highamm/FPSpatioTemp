test_that("sampler returns a data frame with the appropriate number of missing responses", {
  sim_obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4,
                            betavec = 3, sp_de = 0.5, sp_range = 4,
                            sp_ie = 0.5, t_de = 0.5, t_range = 0.7,
                            t_ie = 0.5, spt_ie = 0.5)
  samp_obj <- sample_spatiotemp(sim_obj, n = 100, samp_type = "random")
  expect_equal(sum(is.na(samp_obj$df_full$response_na)), 20)
})

test_that("sampler returns appropriate number of missing values for stratified samples", {
  sim_obj <- sim_spatiotemp(nx = 10, ny = 5, ntime = 4,
                            betavec = 3, sp_de = 0.5, sp_range = 4,
                            sp_ie = 0.5, t_de = 0.5, t_range = 0.7,
                            t_ie = 0.5, spt_ie = 0.5)
  samp_spacestrat <- sample_spatiotemp(sim_obj, n = 50, samp_type = "space_strat")
  samp_timestrat <- sample_spatiotemp(sim_obj, n = 12, samp_type = "time_strat")
  expect_equal(sum(is.na(samp_spacestrat$df_full$response_na)),
               150)
  expect_equal(sum(is.na(samp_timestrat$df_full$response_na)),
               188)
  
  n_per_space <- samp_spacestrat$df_full|>
    dplyr::filter(!is.na(response_na)) |>
    dplyr::count(xcoords, ycoords) |> dplyr::pull(n)
  expect_equal(n_per_space, rep(1, 50))
  
  n_per_time <- samp_timestrat$df_full |>
    dplyr::filter(!is.na(response_na)) |>
    dplyr::count(times) |> dplyr::pull(n)
  expect_equal(n_per_time, rep(3, 4))
})
