#' Sampling of Spatiotemporal Data

#' Randomly selects some rows to be sampled from a \code{sim_spatiotemp()}
#' data frame, and returns a new column that is the \code{response} 
#' value for any row that did get sampled and \code{NA} for any row that
#' did not get sampled.
#' 
#' @param df is a data frame from \code{sim_spatiotemp()}
#' @param n is the sample size
#' @param samp_type can be \itemize{
#'     \item \code{"random"} for randomly selected observations
#'     \item \code{"time_strat"} for the same number of observations to be selected in each time point (rounded up if N / n_time is not an integer), or
#'     \item \code{"space_strat"} for the same number of observations to be selected at each spatial location (rounded up if N / n_space is not an integer).
#'     }
#' @param seed a seed.
#' @return a list with (1) the same data frame that was input to the function with an appended column of the response `response_na` that has `NA` values for the observations that were not sampled and (2) the seed.
#' @examples 
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'       sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'       sigma_nugget_spacetime = 0.5)
#'       
#' sample_spatiotemp(df = obj$out_df, n = 70, samp_type = "random",
#'       seed = obj$seed)
#' @import dplyr
#' @export sample_spatiotemp

sample_spatiotemp <- function(df, n = 100, samp_type = "random", seed) {
  
  set.seed(seed)
  
  if (samp_type == "random") {
    df_samp <- df %>% sample_n(size = n)
  } else if (samp_type == "time_strat") {
    df_samp <- df %>% group_by(times) %>% sample_n(size = ceiling(n / length(unique(df$times))))
  } else if (samp_type == "space_strat") {
    df_samp <- df %>% group_by(xcoords, ycoords) %>% sample_n(size = ceiling(n / length(unique(cbind(df$xcoords, df$ycoords)))))
  }
  
  df_unsamp <- anti_join(df, df_samp)
  df_unsamp$response <- NA
  
  df_missing <- bind_rows(df_samp, df_unsamp)
  
  df_full <- inner_join(df, df_missing,
                        by = c("times", "xcoords", "ycoords")) %>%
    rename(response = response.x, response_na = response.y)
  ret_obj <- list(df_full = df_full, seed = seed)
  return(ret_obj)
}
