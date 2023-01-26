#' Order Spatiotemporal Data
#' 
#' Order a spatiotemporal data set by space within time.
#' 
#' @param data is the name of a data.frame or tibble.
#' @param xcoord is the name of the column with the spatial x-coordinates.
#' @param ycoord is the name of the column with the spatial y-coordinates.
#' @param tcoord is the name of the column with the temporal indeces.
#' 
#' @return a list with \itemize{
#'   \item the ordered data set, with an added logical variable, .observed, that
#'   is TRUE if the row was observed in the original data set.
#'   \item the large spatial distance matrix
#'   \item the small spatial distance matrix
#'   \item the large temporal distance matrix
#'   \item the small temporal distance matrix
#'   }
#' @import stats
#' @export order_spt


order_spt <- function(data, xcoord, ycoord, tcoord) {

  # browser()
  
  # tcoord <- substitute(tcoord)
  # xcoord <- substitute(xcoord)
  # ycoord <- substitute(ycoord)
  
  data_ord <- data |> dplyr::filter(!is.na(.data[[xcoord]]) & !is.na(.data[[ycoord]]) & !is.na(.data[[tcoord]])) |>
    dplyr::arrange(.data[[tcoord]], .data[[xcoord]], .data[[ycoord]])
  
  key_t <- data_ord |> dplyr::distinct(.data[[tcoord]]) |>
    dplyr::arrange(.data[[tcoord]])
  key_t <- key_t |>
    dplyr::mutate(`_tindex` = seq.int(from = 1, to = nrow(key_t)))
  
  h_t_small <- stats::dist(key_t[[tcoord]], diag = TRUE, upper = TRUE) |>
    as.matrix()
  
  h_t_large <- stats::dist(data_ord |> dplyr::select(dplyr::any_of(tcoord)),
                           diag = TRUE, upper = TRUE) |>
    as.matrix()
  
  # record the number of unique temporal observations
  n_t <- nrow(key_t)
  
  key_sp <- data_ord |> dplyr::distinct(.data[[xcoord]], .data[[ycoord]]) |>
    dplyr::arrange(.data[[xcoord]], .data[[ycoord]])
  key_sp <- key_sp |>
    dplyr::mutate(`_spindex` = seq.int(from = 1, to = nrow(key_sp)))
  
  h_sp_small <- stats::dist(key_sp |> dplyr::select(-"_spindex"),
                            diag = TRUE, upper = TRUE) |>
    as.matrix()
  
  h_sp_large <- stats::dist(data_ord |> dplyr::select(dplyr::any_of(c(xcoord, ycoord))),
                            diag = TRUE, upper = TRUE) |>
    as.matrix()
  
  n_sp <- nrow(key_sp)
  data <- dplyr::full_join(data_ord, key_sp) |> dplyr::full_join(key_t)
  
  
  full_grid <- tidyr::expand_grid(`_spindex` = key_sp$`_spindex`,
                                  `_tindex` = key_t$`_tindex`) |>
    dplyr::mutate(`_index` = seq.int(from = 1, to = n_sp * n_t))
  
  full_data <- dplyr::left_join(full_grid, data) |>
    dplyr::arrange(.data$`_index`) |>
    dplyr::mutate(.observed = !is.na(.data[[tcoord]]) & !is.na(.data[[xcoord]]) & !is.na(.data[[ycoord]])) |>
    dplyr::arrange(.data[[tcoord]], .data[[xcoord]], .data[[ycoord]])
  
  spt_obj <- list(full_data = full_data, h_sp_large = h_sp_large,
                  h_sp_small = h_sp_small,
                  h_t_large = h_t_large, h_t_small = h_t_small)
  return(spt_obj)
  
}
