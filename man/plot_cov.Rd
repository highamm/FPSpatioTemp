% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_funs.R
\name{plot_cov}
\alias{plot_cov}
\title{Covariance Plot}
\usage{
plot_cov(stlmfit_obj, sp_epstol = NULL, t_max = NULL, ...)
}
\arguments{
\item{stlmfit_obj}{an object fit with \code{stlmfit()}}

\item{sp_epstol}{a vector of the values used for colours for various spatial distances. The default is the a vector of \code{0}, the minimum spatial distance between two sites, the maximum spatial distance between two sites, and \code{Inf}.}

\item{t_max}{the maximum value for temporal distance to be plotted on the x-axis. The default is the maximum temporal distance between two time points in the data.}

\item{...}{extra options to be passed to `\code{ggplot}`}
}
\value{
a plot with covariance on the y-axis and temporal distance on the x-axis
coloured by various spatial distances in \code{sp_epstol}.
}
\description{
Produces a plot of the covariance with time distance on the x-axis for various
degrees of spatial distance.
}
\examples{
obj <- stlmfit(formula = response_na ~ x, data = samp_data,
 xcoord = "xcoords", ycoord = "ycoords", tcoord = "times") 
plot_cov(obj)
plot_cov(stlmfit_obj = obj, sp_epstol = c(0.2, 0.4, 1.2, 2, Inf),
t_max = 1.5)
}
