% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_spatiotemp.R
\name{sample_spatiotemp}
\alias{sample_spatiotemp}
\title{Sampling of Spatiotemporal Data
Randomly selects some rows to be sampled from a \code{sim_spatiotemp()}
data frame, and returns a new column that is the \code{response} 
value for any row that did get sampled and \code{NA} for any row that
did not get sampled.}
\usage{
sample_spatiotemp(obj, n = 100, samp_type = "random")
}
\arguments{
\item{obj}{is generated from \code{sim_spatiotemp}}

\item{n}{is the sample size}

\item{samp_type}{can be \itemize{
\item \code{"random"} for randomly selected observations
\item \code{"time_strat"} for the same number of observations to be selected in each time point (rounded up if N / n_time is not an integer), or
\item \code{"space_strat"} for the same number of observations to be selected at each spatial location (rounded up if N / n_space is not an integer).
}}
}
\value{
a list with (1) the same data frame that was input to the function with an appended column of the response `response_na` that has `NA` values for the observations that were not sampled and (2) the seed.
}
\description{
Sampling of Spatiotemporal Data
Randomly selects some rows to be sampled from a \code{sim_spatiotemp()}
data frame, and returns a new column that is the \code{response} 
value for any row that did get sampled and \code{NA} for any row that
did not get sampled.
}
\examples{
obj <- sim_spatiotemp(nx = 6, ny = 5, ntime = 4, betavec = 3,
      sp_de = 0.5, sp_range = 4, sp_ie = 0.5,
      t_de = 0.5, t_range = 0.7, t_ie = 0.5,
      spt_ie = 0.5)
      
sample_spatiotemp(obj = obj, n = 70, samp_type = "random")
sample_spatiotemp(obj = obj, n = 40, samp_type = "time_strat")
sample_spatiotemp(obj = obj, n = 60, samp_type = "space_strat")

}
