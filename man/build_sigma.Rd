% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_sigma.R
\name{build_sigma}
\alias{build_sigma}
\title{Build Spatio-temporal Covariance Matrix}
\usage{
build_sigma(
  sp_de = 0.9,
  sp_ie = 0.1,
  t_de = 0.7,
  t_ie = 0.3,
  spt_de = 0.4,
  spt_ie = 0.5,
  R_sp,
  R_t,
  Z_sp,
  Z_t,
  model_type = "product_sum"
)
}
\arguments{
\item{sp_de}{is the spatial dependent error variance (spatial partial sill).}

\item{sp_ie}{is the spatial independent error variance (spatial nugget).}

\item{t_de}{is the temporal dependent error variance (temporal partial sill).}

\item{t_ie}{is the temporal independent error variance (temporal nugget).}

\item{spt_de}{is the spatio-temporal dependent error variance (spatio-temporal partial sill).}

\item{spt_ie}{is the spatio-temporal independent error variance (spatio-temporal nugget).}

\item{R_sp}{is the spatial correlation matrix}

\item{R_t}{is the temporal correlation matrix}

\item{Z_sp}{is the spatial random effects matrix}

\item{Z_t}{is the temporal random effects matrix}

\item{model_type}{is either \code{"product_sum"} (by default) or \code{"sum_with_error"}}
}
\value{
a spatio-temporal covariance matrix.
}
\description{
Construct a spatio-temporal covariance matrix.
}
