% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gr_pois_log_loglik}
\alias{gr_pois_log_loglik}
\title{gr_pois_log_loglik cpp}
\usage{
gr_pois_log_loglik(x, beta, y, time, time_zero, h, type)
}
\arguments{
\item{x}{design matrix}

\item{beta}{vcm vector}

\item{y}{response vector}

\item{time}{vector}

\item{time_zero}{grid point double}

\item{h}{bandwidth double}

\item{type}{kernel function int}
}
\description{
Gradient function for vcm
}
