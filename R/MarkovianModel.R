#########################################################################
########################### MARKOVIAN MODEL ############################
#########################################################################

# Class MarkovianModel --------------------------------------------------
#' @encoding UTF-8
#' @title `MarkovianModel` (abstract R6 class)
#' @description Abstract R6 base class for queueing systems whose
#'   inter-arrival *and* service times follow exponential (Markovian)
#'   distributions.
#'
#' @section Public fields:
#' \describe{
#'   \item{`arrivalDistribution`}{An object of class
#'     \code{distr::Exp} that represents inter-arrival times.}
#'   \item{`serviceDistribution`}{An object of class
#'     \code{distr::Exp} that represents service times.}
#'   \item{`servers`}{\code{integer}. Number of parallel servers.}
#'   \item{`out`}{A \code{list} of performance metrics (\eqn{\rho},
#'     \eqn{L}, \eqn{W}, etc.) produced by subclasses.}
#' }
#'
#' @section Public methods:
#' \describe{
#'   \item{\code{$initialize()}}{Constructor.}
#'   \item{\code{$lambda()} / \code{$mu()}}{Return arrival and service
#'     rates.}
#'   \item{\code{$Pn(n)}}{Steady‑state probability
#'     \eqn{\Pr\{N = n\}}.}
#'   \item{\code{$FW(x)}}{CDF of the time *in the system*,
#'     \eqn{F_W(x)}.}
#'   \item{\code{$FWq(x)}}{CDF of the *waiting* time in queue,
#'     \eqn{F_{W_q}(x)}.}
#'   \item{\code{$Qn(n)}}{Probability that the queue length equals
#'     \eqn{n}.}
#'   \item{\code{$maxCustomers()}}{Practical upper bound for the number
#'     of customers the model can hold (may be \eqn{\infty}).}
#'   \item{\code{$print()}}{Pretty printer for the console.}
#' }
#'
#' @param arrivalDistribution Object created with
#'   \code{distr::Exp(rate)}.
#' @param serviceDistribution Object created with
#'   \code{distr::Exp(rate)}.
#' @param n,x,... Arguments passed on to the corresponding methods (see
#'   details above).
#'
#' @seealso \code{\link{MM1}}, \code{\link{Pn}}, \code{\link{FW}}
#'
#' @examples
#' m <- MM1$new(4, 6)
#' m$Pn(0:3)
#' FWq(m, 1)
#'
#' @usage NULL
#' @format An \code{R6ClassGenerator} object.
#' @docType class
#' @import R6
#' @importFrom distr Exp rate
#' @export
MarkovianModel <- R6::R6Class(
  classname = "MarkovianModel",
  public = list(
    # --------------------------- Fields --------------------------------#
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    servers = NA_integer_,
    out = list(),

    # ------------------------ Constructor ------------------------------#
    initialize = function(arrivalDistribution = distr::Exp(1),
                          serviceDistribution = distr::Exp(1)) {
      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
    },

    # ---------------------------- Rates --------------------------------#
    lambda = function() distr::rate(self$arrivalDistribution),
    mu     = function() distr::rate(self$serviceDistribution),

    # ------------------------ Abstract API -----------------------------#
    Pn  = function(n) stop("Pn not implemented."),
    FW  = function(x) stop("FW not implemented."),
    FWq = function(x) stop("FWq not implemented."),
    Qn  = function(n) stop("Qn not implemented."),
    maxCustomers = function() Inf,

    # --------------------------- Printer -------------------------------#
    print = function(...) {
      cat("Model:", class(self)[1], "\n")
      if (length(self$out))
        cat("rho =", round(self$out$rho, 4),
            "\tL =", round(self$out$l, 4),
            "\tW =", round(self$out$w, 4), "\n")
    }
  )
)

#########################################################################
## Convenience wrapper functions (for S3 back‑compatibility) -----------#
#########################################################################

#' @encoding UTF-8
#' @title Steady-state probability \eqn{\Pr\{N = n\}}
#' @description S3 frontend that delegates to \code{qm$Pn(n)}.
#' @param qm An object that inherits from \code{MarkovianModel}.
#' @param n  Non‑negative integer vector.
#' @export
Pn <- function(qm, n) qm$Pn(n)

#' @encoding UTF-8
#' @title CDF of the time in system \eqn{F_W(x)}
#' @inheritParams Pn
#' @param x Non‑negative numeric vector.
#' @export
FW <- function(qm, x) qm$FW(x)

#' @encoding UTF-8
#' @title CDF of the waiting time in queue \eqn{F_{W_q}(x)}
#' @inheritParams FW
#' @export
FWq <- function(qm, x) qm$FWq(x)

#' @encoding UTF-8
#' @title Probability that the queue length equals \eqn{n}
#' @inheritParams Pn
#' @export
Qn <- function(qm, n) qm$Qn(n)

#' @encoding UTF-8
#' @title Upper bound on the number of customers supported
#' @inheritParams Pn
#' @export
maxCustomers <- function(qm) qm$maxCustomers()
