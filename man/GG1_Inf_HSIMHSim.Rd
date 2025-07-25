% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulateModels.R
\name{GG1_Inf_HSIMHSim}
\alias{GG1_Inf_HSIMHSim}
\title{GG1_Inf_HSIMHSim – simulated G/G/1/Inf/H queue (R6)}
\arguments{
\item{arrivalDistribution, serviceDistribution}{Objects from package
\strong{distr} with the inter-arrival and service-time laws.}

\item{H}{Integer >= 1, size of the customer population.}

\item{staClients}{Warm-up customers discarded from statistics.}

\item{nClients}{Customers collected for statistics.}

\item{historic}{Logical, store the whole trajectory?}
}
\description{
Discrete-event simulation of a single-server queue fed
by a \emph{finite} population of \code{H} sources.  At any moment each source
is either \strong{in} the system (being served or waiting) or \strong{outside}
and generating its own inter-arrival time.  The total population is
constant and no arrivals are lost.
}
