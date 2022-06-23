#' Sampling of Spatiotemporal Data

#' Randomly selects some rows to be sampled from a \code{sim_spatiotemp()}
#' data frame, and returns a new column that is the \code{response} 
#' value for any row that did get sampled and \code{NA} for any row that
#' did not get sampled.
#' 
#' @param obj is generated from \code{sim_spatiotemp}
#' @param n is the sample size
#' @param samp_type can be \itemize{
#'     \item \code{"random"} for randomly selected observations
#'     \item \code{"time_strat"} for the same number of observations to be selected in each time point (rounded up if N / n_time is not an integer), or
#'     \item \code{"space_strat"} for the same number of observations to be selected at each spatial location (rounded up if N / n_space is not an integer).
#'     }
#' @return a list with (1) the same data frame that was input to the function with an appended column of the response `response_na` that has `NA` values for the observations that were not sampled and (2) the seed.
#' @examples 
#' obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
#'       sigma_parsil_spat = 0.5, range = 4, sigma_nugget_spat = 0.5,
#'       sigma_parsil_time = 0.5, rho = 0.7, sigma_nugget_time = 0.5,
#'       sigma_nugget_spacetime = 0.5)
#'       
#' sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
#' @importFrom dplyr sample_n group_by anti_join inner_join bind_rows
#' @export sample_spatiotemp

sample_spatiotemp <- function(obj, n = 100, samp_type = "random") {
  
  if (inherits(obj, "sim_spatiotemp") == FALSE) {
    stop("obj must be a sim_spatiotemp class object.")
  }
  
  
  df <- obj$out_df
  
  ## need to come up with a better fix for this:
  ## if not there, gives a warning about no binding
  times <- obj$out_df$times
  xcoords <- obj$out_df$xcoords
  ycoords <- obj$out_df$ycoords
  response <- obj$out_df$response
  
  df <- tibble(times, xcoords, ycoords, response)
  
  set.seed(obj$seed)
  
  if (samp_type == "random") {
    df_samp <- df |> dplyr::sample_n(size = n)
  } else if (samp_type == "time_strat") {
    df_samp <- df |> group_by(times) |> dplyr::sample_n(size = ceiling(n / length(unique(df$times))))
  } else if (samp_type == "space_strat") {
    df_samp <- df |> dplyr::group_by(xcoords, ycoords) |> dplyr::sample_n(size = ceiling(n / length(unique(cbind(df$xcoords, df$ycoords)))))
  }
  
  df_unsamp <- dplyr::anti_join(df, df_samp)
  df_unsamp$response <- NA
  
  df_missing <- dplyr::bind_rows(df_samp, df_unsamp)
  
  df_full <- dplyr::inner_join(df, df_missing,
                        by = c("times", "xcoords", "ycoords"))
  
  names(df_full) <- c("times", "xcoords", "ycoords", "response", "response_na")
  
  ret_obj <- list(df_full = df_full, seed = obj$seed)
  return(ret_obj)
}
