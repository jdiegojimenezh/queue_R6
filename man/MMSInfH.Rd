% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComplexModels.R
\name{MMSInfH}
\alias{MMSInfH}
\title{MMSInfH-class: M/M/s/Inf/H finite-population queue}
\arguments{
\item{lambda}{Mean arrival rate (> 0).}

\item{mu}{Mean service rate  (> 0).}

\item{s}{Servers (integer >= 1).}

\item{h}{Population size (>= s).}
}
\description{
MMSInfH-class: M/M/s/Inf/H finite-population queue
}
\section{Public methods}{

\describe{
\item{\code{$initialize()}}{Constructor.}
\item{\code{$lambda()} / \code{$mu()}}{Return arrival and service
rates.}
\item{\code{$Pn(n)}}{Steady‑state probability
\eqn{\Pr\{N = n\}}.}
\item{\code{$FW(x)}}{CDF of the time \emph{in the system},
\eqn{F_W(x)}.}
\item{\code{$FWq(x)}}{CDF of the \emph{waiting} time in queue,
\eqn{F_{W_q}(x)}.}
\item{\code{$Qn(n)}}{Probability that the queue length equals
\eqn{n}.}
\item{\code{$maxCustomers()}}{Practical upper bound for the number
of customers the model can hold (may be \eqn{\infty}).}
\item{\code{$print()}}{Pretty printer for the console.}
}
}

