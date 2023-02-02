#' TOC Moose Survey Data - Small
#'
#' A data set that can be used with the \code{FPSpatioTemp} package. The data set is a subset of \code{moose_complete}, to be used in the vignette.
#'
#' @format A data frame object including:
#' \describe{
#'   \item{count}{Total moose observed at the site, with NA for sites that were not sampled}
#'   \item{xcoords}{site x-coordinate, in TM}
#'   \item{ycoords}{site y-coordinate, in TM}
#'   \item{year}{Year of the survey}
#'   \item{strata}{strata (\code{LOW} or \code{HIGH}), as a factor, corresponding to the original stratification in each year}
#'   \item{area_mi}{Area of the site, in square miles}
#'   \item{elev_mean}{Average elevation in the site}
#'   \item{samp_frame}{whether or not the site was in the sampling frame for that year}
#' }
#' 
#' @examples
#' data(moose_vignette)
#' names(moose_vignette)
#' summary(moose_vignette)
"moose_vignette"

#' Example Data
#' 
#' A toy data set that is small in size used for many examples and tests.
#' 
#' @format A data frame object including:
#' \describe{
#'   \item{times}{a vector of time points}
#'   \item{xcoords}{a vector of x-coordinates}
#'   \item{ycoords}{a vector of y-coordinates}
#'   \item{response}{a response variable with no missing values}
#'   \item{response_na}{a response variable with missing values}
#'   \item{predwts}{a vector of prediction weights}
#'   \item{x}{a covariate}
#'   \item{area}{the area of the site}
#'   \item{x_fact}{a predictor that is of type \code{<fct>}}
#'   \item{x_char}{a predictor that is of type \code{<chr>}}
#'   \item{x_miss}{x with a few NA values}
#'   \item{x_fact_miss}{x_fact with a few NA valus}
#'   \item{xcoords_miss}{xcoords with a few NA values}
#'   \item{times_miss}{times with a few NA values}
#'}
#' @examples
#' data(samp_data)
#' names(samp_data)
#' summary(samp_data)
"samp_data"

#' TOC Moose Survey Data 2014 - 2020 Subset
#'
#' A data set that can be used with the \code{FPSpatioTemp} package.
#'
#' @format A data frame object including:
#' \describe{
#'   \item{totalmoosena}{Total moose observed at the site, with NA for sites that were not sampled}
#'   \item{xcoords}{site x-coordinate, in TM}
#'   \item{ycoords}{site y-coordinate, in TM}
#'   \item{Surveyyear}{Year of the survey}
#'   \item{newstrat}{Stratification variable with the most recent stratification in 2020 for every year}
#'   \item{strat_origin}{Stratification variable with the original stratification for each year}
#'   \item{samp_frame}{whether or not the site was in the sampling frame for that year}
#'   \item{yearind}{whether or not the Surveyyear variable is equal to 2020}
#' }
#' 
#' #' @examples
#' data(moose_14_20)
#' names(moose_14_20)
#' summary(moose_14_20)
"moose_14_20"

