% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulateModels.R
\name{GG1}
\alias{GG1}
\title{Simulate a G/G/1 queue}
\usage{
GG1(
  arrivalDistribution = Exp(3),
  serviceDistribution = Exp(6),
  staClients = 100L,
  nClients = 1000L,
  historic = FALSE,
  nsim = 10L,
  nproc = 1L
)
}
\arguments{
\item{arrivalDistribution}{arrival distribution (object from package
\emph{distr}).}

\item{serviceDistribution}{service‑time distribution (object from
\emph{distr}).}

\item{staClients}{integer, number of customers discarded as burn‑in
(stabilisation stage).}

\item{nClients}{integer, number of customers collected for statistics.}

\item{historic}{logical, record evolution of the statistics.}

\item{nsim}{integer, number of independent replications.}

\item{nproc}{integer, CPU cores to use.  If \code{nproc = 1} the function
runs sequentially.}
}
\value{
If \code{nsim == 1} a single \code{GG1Sim} object; otherwise an object
of the same class containing aggregated statistics (\code{mean}, \code{sd},
etc.) produced by \code{combineSimulations()}.
}
\description{
Convenience wrapper around \code{GG1Sim$new()} plus parallel replication.
}
