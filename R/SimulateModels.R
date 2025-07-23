#' GG1Sim – simulated G/G/1 queue (R6)
#'
#' @description R6 class that simulates a single‑server queue with a
#'   general (i.i.d.) arrival distribution and a general service
#'   distribution – commonly referred to as a *G/G/1* system.  The
#'   algorithm is exactly the same as the legacy S3 implementation
#'   but the results are stored inside the object under field `out`.
#'
#' @section Stored metrics (slot `out`):
#' * `pn`   – empirical steady‑state probability vector `P{N = n}`.
#' * `l`    – mean number of customers in the system *L*.
#' * `lq`   – mean number of customers in the queue *Lq*.
#' * `w`    – mean waiting time in the system *W*.
#' * `wq`   – mean waiting time in the queue *Wq*.
#' * `eff`  – empirical efficiency `W / (W - Wq)`.
#' * `rho`  – empirical traffic intensity `L - Lq`.
#' * `historic` (optional) – matrix with the evolution of the variables
#'   during the run when `historic = TRUE`.
#'
#' @param arrivalDistribution arrival distribution (object from package
#'   *distr*).
#' @param serviceDistribution service‑time distribution (object from
#'   *distr*).
#' @param staClients integer, number of customers discarded as burn‑in
#'   (stabilisation stage).
#' @param nClients integer, number of customers collected for statistics.
#' @param historic logical, record evolution of the statistics.
#' @name GG1Sim
#' @export
GG1Sim <- R6::R6Class(
  "GG1Sim",
  public = list(
    #•• object fields --------------------------------------------------
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,   # list filled at the end

    #•• constructor ----------------------------------------------------
    initialize = function(arrivalDistribution, serviceDistribution,
                          staClients = 100L, nClients = 1000L,
                          historic = FALSE) {
      # basic checks ---------------------------------------------------
      if (!belong(arrivalDistribution, "distr"))
        stop("'arrivalDistribution' must be a 'distr' object", call. = FALSE)
      if (!belong(serviceDistribution, "distr"))
        stop("'serviceDistribution' must be a 'distr' object", call. = FALSE)
      if (!is.numeric(staClients) || staClients < 0)
        stop("'staClients' must be >= 0", call. = FALSE)
      if (!is.numeric(nClients)  || nClients  <= 0)
        stop("'nClients' must be > 0", call. = FALSE)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    # print method -----------------------------------------------------
    print = function(...) {
      cat("Model: GG1Sim\n")
      if (is.list(self$out$l)) {
        df <- data.frame(
          mean = c(self$out$l$mean,  self$out$lq$mean,
                   self$out$w$mean, self$out$wq$mean,
                   self$out$rho$mean, self$out$eff$mean),
          sd   = c(self$out$l$sd,  self$out$lq$sd,
                   self$out$w$sd, self$out$wq$sd,
                   self$out$rho$sd, self$out$eff$sd),
          row.names = c("L", "Lq", "W", "Wq", "Intensity", "Efficiency")
        )
        print(df)
      } else {
        cat(sprintf("\nL =\t%.5g\tW =\t%.5g\tIntensity =\t%.5g\n",
                    self$out$l, self$out$w, self$out$rho))
        cat(sprintf("Lq =\t%.5g\tWq =\t%.5g\tEfficiency =\t%.5g\n\n",
                    self$out$lq, self$out$wq, self$out$eff))
      }
      invisible(self)
    }
  ),

  private = list(
    #•• core discrete‑event simulation --------------------------------
    simulate_core = function() {
      n_sta <- self$staClients
      n_sim <- self$nClients

      # pre‑generate inter‑arrival and service times -------------------
      tArr  <- r(self$arrivalDistribution)(n_sta + n_sim)
      tServ <- r(self$serviceDistribution)(n_sta + n_sim)
      if (any(tArr  < 0)) stop("Arrival distribution produced negatives")
      if (any(tServ < 0)) stop("Service distribution produced negatives")

      # state variables ------------------------------------------------
      iArr <- iServ <- 1L
      sys  <- sim  <- 0L
      a    <- 0       # next arrival
      b    <- -1      # next departure (‑1 means server idle)
      cron <- c_sum <- d_sum <- 0

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = n_sta + n_sim, ncol = 7,
                       dimnames = list(NULL, c("L", "Lq", "W", "Wq",
                                               "Clients", "Intensity",
                                               "tClient")))
      }

      # -------- stabilisation stage ----------------------------------
      while (sim < n_sta) {
        tMin <- if (sys > 0) min(a, b) else a
        cron <- cron + tMin

        if (tMin == a) {                # arrival event
          sim  <- sim + 1L
          if (sys == 0) {
            b <- tServ[iServ]; iServ <- iServ + 1L
          } else {
            c_sum <- c_sum + tMin * sys
            d_sum <- d_sum + tMin * (sys - 1)
            b     <- b - tMin
          }
          a   <- tArr[iArr]; iArr <- iArr + 1L
          sys <- sys + 1L
        } else {                        # departure event
          c_sum <- c_sum + tMin * sys
          d_sum <- d_sum + tMin * (sys - 1)
          sys   <- sys - 1L
          if (sys == 0) {
            b <- -1
          } else {
            b <- tServ[iServ]; iServ <- iServ + 1L
          }
          a <- a - tMin
        }

        if (self$historic) {
          l  <- c_sum / cron
          lq <- d_sum / cron
          w  <- c_sum / sim
          wq <- d_sum / sim
          hist[sim, ] <- c(l, lq, w, wq, sys, l - lq, cron)
        }
      }

      # -------- main simulation stage --------------------------------
      acum_sta <- cron
      sim  <- 0L
      cron <- c_sum <- d_sum <- 0
      tnClients <- numeric(n_sim)

      while (sim < n_sim) {
        tMin <- if (sys > 0) min(a, b) else a
        cron <- cron + tMin
        tnClients[sys + 1L] <- tnClients[sys + 1L] + tMin

        if (tMin == a) {                # arrival
          sim <- sim + 1L
          if (sys == 0) {
            b <- tServ[iServ]; iServ <- iServ + 1L
          } else {
            c_sum <- c_sum + tMin * sys
            d_sum <- d_sum + tMin * (sys - 1)
            b     <- b - tMin
          }
          sys <- sys + 1L
          a   <- tArr[iArr]; iArr <- iArr + 1L
        } else {                        # departure
          c_sum <- c_sum + tMin * sys
          d_sum <- d_sum + tMin * (sys - 1)
          sys   <- sys - 1L
          if (sys == 0) {
            b <- -1
          } else {
            b <- tServ[iServ]; iServ <- iServ + 1L
          }
          a <- a - tMin
        }

        if (self$historic && sim > 0) {
          l  <- c_sum / cron
          lq <- d_sum / cron
          w  <- c_sum / sim
          wq <- d_sum / sim
          hist[sim + n_sta, ] <- c(l, lq, w, wq, sys, l - lq,
                                   acum_sta + cron)
        }
      }

      # -------- statistics -------------------------------------------
      l  <- c_sum / cron
      lq <- d_sum / cron
      w  <- c_sum / n_sim
      wq <- d_sum / n_sim
      rho <- l - lq
      eff <- w / (w - wq)

      min_prob <- which(tnClients > (.Machine$double.eps ^ 0.5))
      if (length(min_prob) > 0) {
        pn <- tnClients[seq_len(max(min_prob))] / cron
      } else {
        pn <- c(0, 0)
      }

      if (self$historic) {
        self$out <- list(historic = hist, l = l, lq = lq, w = w, wq = wq,
                         pn = pn, rho = rho, eff = eff)
      } else {
        self$out <- list(l = l, lq = lq, w = w, wq = wq,
                         pn = pn, rho = rho, eff = eff)
      }
    }
  )
)

# ---------------------------------------------------------------------
# Functional constructor – user‑friendly wrapper ----------------------
# ---------------------------------------------------------------------
#' Simulate a G/G/1 queue
#'
#' Convenience wrapper around `GG1Sim$new()` plus parallel replication.
#'
#' @inheritParams GG1Sim
#' @param nsim  integer, number of independent replications.
#' @param nproc integer, CPU cores to use.  If `nproc = 1` the function
#'   runs sequentially.
#' @return If `nsim == 1` a single `GG1Sim` object; otherwise an object
#'   of the same class containing aggregated statistics (`mean`, `sd`,
#'   etc.) produced by `combineSimulations()`.
#' @export
GG1 <- function(arrivalDistribution = Exp(3),
                serviceDistribution = Exp(6),
                staClients = 100L,
                nClients   = 1000L,
                historic   = FALSE,
                nsim       = 10L,
                nproc      = 1L) {

  if (!is.numeric(nsim)  || nsim  <= 0) stop("'nsim' must be > 0")
  if (!is.numeric(nproc) || nproc <= 0) stop("'nproc' must be > 0")

  sim_fun <- function(arrivalDistribution, serviceDistribution,
                      staClients, nClients, historic) {
    GG1Sim$new(arrivalDistribution, serviceDistribution,
               staClients, nClients, historic)
  }

  params <- list(arrivalDistribution, serviceDistribution,
                 staClients, nClients, historic)

  ParallelizeSimulations(sim_fun, params, nsim, nproc) |>
    (
      function(res) {
        if (inherits(res, "GG1Sim")) return(res) else combineSimulations(res)
      }
    )()
}

########################################################################
## GGSSim  –  simulated G/G/s queue (R6)                              ##
########################################################################

#' GGSSim – simulated multiserver G/G/s queue (R6)
#'
#' @description R6 class that simulates a queue with *s* identical
#'   servers, general i.i.d. inter–arrival and service–time
#'   distributions (G/G/s).  It is a direct R6 port of the original
#'   S3 function `G_G_S()`.  Results are stored in field `out`.
#'
#' @section Stored metrics (slot `out`):
#' * `pn`, `l`, `lq`, `w`, `wq`, `rho`, `eff` – as in `GG1Sim`.
#' * `historic` matrix is present when `historic = TRUE`.
#'
#' @param arrivalDistribution,serviceDistribution  *distr* objects.
#' @param servers     integer, number of servers *s* (>= 1).
#' @param staClients  integer, warm-up customers.
#' @param nClients    integer, customers collected for stats.
#' @param historic    logical, record full trajectory?
#' @name GGSSim
#' @export
GGSSim <- R6::R6Class(
  "GGSSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    servers             = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          servers     = 2L,
                          staClients  = 100L,
                          nClients    = 1000L,
                          historic    = FALSE) {

      # argument checks ------------------------------------------------
      if (!belong(arrivalDistribution, "distr"))
        stop("'arrivalDistribution' must be a 'distr' object", call. = FALSE)
      if (!belong(serviceDistribution, "distr"))
        stop("'serviceDistribution' must be a 'distr' object", call. = FALSE)
      if (!is.numeric(servers)  || servers  < 1)
        stop("'servers' must be >= 1", call. = FALSE)
      if (!is.numeric(staClients) || staClients < 0)
        stop("'staClients' must be >= 0", call. = FALSE)
      if (!is.numeric(nClients)  || nClients  <= 0)
        stop("'nClients' must be > 0", call. = FALSE)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$servers             <- as.integer(servers)
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GGSSim  (G/G/", self$servers, ")\n", sep = "")
      with(self$out, {
        cat(sprintf("L   = %.5f | Lq  = %.5f\n", l,  lq))
        cat(sprintf("W   = %.5f | Wq  = %.5f\n", w,  wq))
        cat(sprintf("rho = %.5f | Eff = %.5f\n", rho, eff))
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      s   <- self$servers
      n_w <- self$staClients
      n_m <- self$nClients
      n_total <- n_w + n_m

      # pre-generate times --------------------------------------------
      tArr  <- r(self$arrivalDistribution)(n_total)
      tServ <- r(self$serviceDistribution)(n_total)
      if (any(tArr < 0) || any(tServ < 0))
        stop("distributions produced negative times")

      # state variables ----------------------------------------------
      iArr <- 1L; iServ <- 0L
      busy <- integer()            # indices of customers in service
      sys  <- 0L                   # customers in system
      simC <- 0L                   # counted customers
      t    <- 0                    # clock

      cCum <- dCum <- 0            # integral of L and Lq

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = n_total, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      # ---------------- warm-up phase -------------------------------
      while (simC < n_w) {
        nextArr  <- tArr[iArr]
        nextServ <- if (length(busy)) min(tServ[busy]) else Inf
        dt       <- min(nextArr, nextServ)

        # advance clock and areas
        t    <- t + dt
        cCum <- cCum + dt * sys
        if (sys > s) dCum <- dCum + dt * (sys - s)

        # decrement residual times
        if (length(busy)) tServ[busy] <- tServ[busy] - dt
        tArr[iArr] <- tArr[iArr] - dt

        if (nextArr <= nextServ) {            # arrival
          simC <- simC + 1L
          sys  <- sys + 1L
          iArr <- iArr + 1L
          if (sys <= s) {                     # start service
            iServ <- iServ + 1L; busy <- c(busy, iServ)
          }
        } else {                              # departure
          idx  <- busy[ which.min(tServ[busy]) ]
          busy <- setdiff(busy, idx)
          sys  <- sys - 1L
          if (sys >= s) {                     # pull next waiting job
            iServ <- iServ + 1L; busy <- c(busy, iServ)
          }
        }

        if (self$historic) {
          L  <- cCum / t
          Lq <- dCum / t
          W  <- if (simC) cCum / simC else NA
          Wq <- if (simC) dCum / simC else NA
          hist[simC, ] <- c(L,Lq,W,Wq,sys,(L-Lq)/s,t)
        }
      }

      tWarm <- t
      cCum <- dCum <- 0
      simC <- 0L
      pnVec <- numeric(n_m)

      # ---------------- main phase ----------------------------------
      while (simC < n_m) {
        nextArr  <- tArr[iArr]
        nextServ <- if (length(busy)) min(tServ[busy]) else Inf
        dt       <- min(nextArr, nextServ)

        # advance clock and areas
        t    <- t + dt
        cCum <- cCum + dt * sys
        if (sys > s) dCum <- dCum + dt * (sys - s)
        pnVec[min(sys, n_m) + 1L] <- pnVec[min(sys, n_m) + 1L] + dt

        # decrement residual times
        if (length(busy)) tServ[busy] <- tServ[busy] - dt
        tArr[iArr] <- tArr[iArr] - dt

        if (nextArr <= nextServ) {            # arrival
          simC <- simC + 1L
          sys  <- sys + 1L
          iArr <- iArr + 1L
          if (sys <= s) {                     # start service
            iServ <- iServ + 1L; busy <- c(busy, iServ)
          }
        } else {                              # departure
          idx  <- busy[ which.min(tServ[busy]) ]
          busy <- setdiff(busy, idx)
          sys  <- sys - 1L
          if (sys >= s) {
            iServ <- iServ + 1L; busy <- c(busy, iServ)
          }
        }

        if (self$historic && simC > 0) {
          L  <- cCum / (t - tWarm)
          Lq <- dCum / (t - tWarm)
          W  <- cCum / simC
          Wq <- dCum / simC
          hist[simC + n_w, ] <- c(L,Lq,W,Wq,sys,(L-Lq)/s,t)
        }
      }

      # statistics ----------------------------------------------------
      L  <- cCum / (t - tWarm)
      Lq <- dCum / (t - tWarm)
      W  <- cCum / n_m
      Wq <- dCum / n_m
      rho <- (L - Lq) / s
      eff <- W / (W - Wq)

      out <- list(l = L, lq = Lq, w = W, wq = Wq,
                  rho = rho, eff = eff,
                  pn = pnVec / (t - tWarm))
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Functional constructor for replication / parallel                  ##
########################################################################

#' Simulate a G/G/s queue (wrapper)
#'
#' @inheritParams GGSSim
#' @param nsim  integer, number of replications.
#' @param nproc integer, CPU cores (1 = sequential).
#' @return `GGSSim` object or aggregated result when `nsim > 1`.
#' @export
GGS <- function(arrivalDistribution = Exp(3),
                serviceDistribution = Exp(6),
                servers      = 2L,
                staClients   = 100L,
                nClients     = 1000L,
                historic     = FALSE,
                nsim         = 10L,
                nproc        = 1L) {

  if (nsim < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  sim_fun <- function(arrivalDistribution,
                      serviceDistribution,
                      servers,
                      staClients,
                      nClients,
                      historic) {
    GGSSim$new(arrivalDistribution,
               serviceDistribution,
               servers,
               staClients,
               nClients,
               historic)
  }

  params <- list(arrivalDistribution,
                 serviceDistribution,
                 servers,
                 staClients,
                 nClients,
                 historic)

  res <- ParallelizeSimulations(sim_fun, params, nsim, nproc)
  if (inherits(res, "GGSSim")) res else combineSimulations(res)
}

########################################################################
## GG1KSim – simulated finite-capacity G/G/1/K queue (R6)             ##
########################################################################

#' GG1KSim – simulated G/G/1/K queue (R6)
#'
#' @description Discrete-event simulation of a single-server queue with
#'   finite capacity `K` (system size = `K + 1`, including the job in
#'   service).  Arrivals finding the system full are **lost**.  This is
#'   a direct R6 refactor of the legacy function `G_G_1_K()`.
#'
#' @param arrivalDistribution,serviceDistribution  objects from package
#'   **distr** defining inter-arrival and service-time laws.
#' @param K           integer, maximum queue size (>= 1).  The system
#'   thus holds at most `K + 1` customers (1 in service, *K* waiting).
#' @param staClients  warm-up customers to discard.
#' @param nClients    accepted customers on which statistics are based.
#' @param historic    logical, store full trajectory?
#' @export
GG1KSim <- R6::R6Class(
  "GG1KSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    K                   = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          K          = 2L,
                          staClients = 100L,
                          nClients   = 1000L,
                          historic   = FALSE) {

      if (!belong(arrivalDistribution, "distr"))
        stop("'arrivalDistribution' must be a distr object", call. = FALSE)
      if (!belong(serviceDistribution, "distr"))
        stop("'serviceDistribution' must be a distr object", call. = FALSE)
      if (!is.numeric(K) || K < 1)
        stop("'K' must be >= 1", call. = FALSE)
      if (!is.numeric(staClients) || staClients < 0)
        stop("'staClients' must be >= 0", call. = FALSE)
      if (!is.numeric(nClients)  || nClients  <= 0)
        stop("'nClients' must be > 0", call. = FALSE)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$K                   <- as.integer(K)
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GG1KSim  (G/G/1/", self$K, ")\n", sep = "")
      with(self$out, {
        cat(sprintf("L   = %.5f | Lq  = %.5f\n", l,  lq))
        cat(sprintf("W   = %.5f | Wq  = %.5f\n", w,  wq))
        cat(sprintf("rho = %.5f | Eff = %.5f\n", rho, eff))
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      K   <- self$K
      n_w <- self$staClients
      n_m <- self$nClients
      n_total <- (n_w + n_m) * 2L          # overshoot to ensure enough times
      nServed <- 0L

      # generate inter-arrival and service times ----------------------
      tArr  <- r(self$arrivalDistribution)(n_total)
      tServ <- r(self$serviceDistribution)(n_total)
      if (any(tArr < 0) || any(tServ < 0))
        stop("distributions produced negative times")

      # state variables ----------------------------------------------
      iArr <- 1L; iServ <- 1L
      sys  <- 0L            # customers in system (0 … K+1)
      simC <- 0L            # accepted customers counted
      t    <- 0

      cCum <- dCum <- 0
      lost <- 0L            # rejected customers (for info only)

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = n_w + n_m, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      # helper: next departure time or Inf if idle
      next_dep <- function() {
        if (sys == 0) Inf else tServ[iServ]
      }

      # ---------------- warm-up phase -------------------------------
      while (simC < n_w) {
        tmin <- min(tArr[iArr], next_dep())
        # advance time and areas
        t    <- t + tmin
        cCum <- cCum + tmin * sys
        if (sys > 1) dCum <- dCum + tmin * (sys - 1)
        # update residual times
        if (sys > 0) tServ[iServ] <- tServ[iServ] - tmin
        tArr[iArr] <- tArr[iArr] - tmin

        if (tArr[iArr] <= 1e-12) {          # arrival event
          if (sys <= K) {                   # accepted
            simC <- simC + 1L
            sys  <- sys + 1L
            iArr <- iArr + 1L
          } else {                          # lost
            lost <- lost + 1L
            iArr <- iArr + 1L
          }
        } else {                            # departure event
          sys  <- sys - 1L
          if (sys > 0) {
            iServ <- iServ + 1L
          }
        }

        if (self$historic) {
          l  <- cCum / t
          lq <- dCum / t
          w  <- if (simC) cCum / simC else NA
          wq <- if (simC) dCum / simC else NA
          hist[simC, ] <- c(l,lq,w,wq,sys,l-lq,t)
        }
      }

      tWarm <- t
      cCum <- dCum <- 0
      simC <- 0L
      pnVec <- numeric(n_m)

      # ---------------- main phase ----------------------------------
      while (simC < n_m) {
        tmin <- min(tArr[iArr], next_dep())
        t    <- t + tmin
        cCum <- cCum + tmin * sys
        if (sys > 1) dCum <- dCum + tmin * (sys - 1)
        pnVec[min(sys, n_m) + 1L] <- pnVec[min(sys, n_m) + 1L] + tmin

        if (sys > 0) tServ[iServ] <- tServ[iServ] - tmin
        tArr[iArr]  <- tArr[iArr]  - tmin

        if (tArr[iArr] <= 1e-12) {          # arrival
          if (sys <= K) {                   # accepted
            simC <- simC + 1L
            sys  <- sys + 1L
            iArr <- iArr + 1L
          } else {                          # lost
            lost <- lost + 1L
            iArr <- iArr + 1L
          }
        } else {                            # departure
          sys  <- sys - 1L
          nServed <- nServed + 1L
          if (sys > 0) iServ <- iServ + 1L
        }

        if (self$historic && simC > 0) {
          l  <- cCum / (t - tWarm)
          lq <- dCum / (t - tWarm)
          w  <- cCum / simC
          wq <- dCum / simC
          hist[simC + n_w, ] <- c(l,lq,w,wq,sys,l-lq,t)
        }
      }

      # statistics ----------------------------------------------------
      lambda_eff <- nServed / (t - tWarm)
      L  <- cCum / (t - tWarm)
      Lq <- dCum / (t - tWarm)
      W  <- cCum / nServed
      Wq <- dCum / nServed
      rho <- L - Lq
      eff <- W / (W - Wq)

      out <- list(l = L, lq = Lq, w = W, wq = Wq,
                  lambda_eff = lambda_eff,
                  rho = rho, eff = eff,
                  pn = pnVec / (t - tWarm),
                  lost = lost)
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Wrapper with replication / parallel                                ##
########################################################################

#' Simulate a G/G/1/K queue (wrapper)
#'
#' @inheritParams GG1KSim
#' @param nsim  integer, number of replications.
#' @param nproc integer, CPU cores (1 = sequential).
#' @return `GG1KSim` object or aggregated result when `nsim > 1`.
#' @export
GG1K <- function(arrivalDistribution = Exp(3),
                 serviceDistribution = Exp(6),
                 K             = 2L,
                 staClients    = 100L,
                 nClients      = 1000L,
                 historic      = FALSE,
                 nsim          = 10L,
                 nproc         = 1L) {

  if (nsim < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  sim_fun <- function(arrivalDistribution,
                      serviceDistribution,
                      K,
                      staClients,
                      nClients,
                      historic) {
    GG1KSim$new(arrivalDistribution,
                serviceDistribution,
                K,
                staClients,
                nClients,
                historic)
  }

  params <- list(arrivalDistribution,
                 serviceDistribution,
                 K,
                 staClients,
                 nClients,
                 historic)

  res <- ParallelizeSimulations(sim_fun, params, nsim, nproc)
  if (inherits(res, "GG1KSim")) res else combineSimulations(res)
}

########################################################################
## GGSKSim – simulated G/G/s/K queue (R6)                             ##
########################################################################

#' GGSKSim – simulated G/G/s/K queue (R6)
#'
#' @description Discrete-event simulation of a *multi-server* queue with
#'   finite capacity **K** (maximum queue length).  The system can hold
#'   at most `s + K` customers (up to *s* in service, at most *K* in the
#'   waiting line).  Arrivals that find the system full are **lost**.
#'
#' @param arrivalDistribution,serviceDistribution objects from package
#'   **distr** defining the inter-arrival and service-time laws.
#' @param servers Integer >= 1 (number of parallel identical servers).
#' @param K       Integer >= 1, maximum queue size (capacity minus
#'   servers).
#' @param staClients Warm-up customers to discard.
#' @param nClients  Customers accepted for statistics.
#' @param historic  Logical, keep full trajectory?
#' @export
GGSKSim <- R6::R6Class(
  "GGSKSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    servers             = NULL,
    K                   = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          servers     = 2L,
                          K           = 3L,
                          staClients  = 100L,
                          nClients    = 1000L,
                          historic    = FALSE) {

      if (!belong(arrivalDistribution, "distr"))
        stop("'arrivalDistribution' must be a distr object", call. = FALSE)
      if (!belong(serviceDistribution, "distr"))
        stop("'serviceDistribution' must be a distr object", call. = FALSE)
      if (!is.numeric(servers) || servers < 1)
        stop("'servers' must be >= 1", call. = FALSE)
      if (!is.numeric(K) || K < 1)
        stop("'K' must be >= 1", call. = FALSE)
      if (!is.numeric(staClients) || staClients < 0)
        stop("'staClients' must be >= 0", call. = FALSE)
      if (!is.numeric(nClients)  || nClients  <= 0)
        stop("'nClients' must be > 0", call. = FALSE)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$servers             <- as.integer(servers)
      self$K                   <- as.integer(K)
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GGSKSim  (G/G/", self$servers, "/", self$K, ")\n", sep = "")
      with(self$out, {
        cat(sprintf("L   = %.5f | Lq  = %.5f\n", l,  lq))
        cat(sprintf("W   = %.5f | Wq  = %.5f\n", w,  wq))
        cat(sprintf("rho = %.5f | Eff = %.5f\n", rho, eff))
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      s   <- self$servers
      K   <- self$K
      n_w <- self$staClients
      n_m <- self$nClients
      n_tot <- (n_w + n_m) * 2L                   # generous cushion

      tArr  <- r(self$arrivalDistribution)(n_tot)
      tServ <- r(self$serviceDistribution)(n_tot)
      if (any(tArr < 0) || any(tServ < 0))
        stop("negative times produced by distributions")

      # --- state variables -------------------------------------------
      iArr <- 1L; iServ <- 1L
      busy <- numeric(0)          # remaining times of jobs in service
      sys  <- 0L                  # customers in system (0 … s+K)
      simC <- 0L                  # accepted customers counted
      t    <- 0

      cCum <- dCum <- 0
      lost <- 0L

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = n_w + n_m, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      next_dep <- function() if (length(busy)) min(busy) else Inf

      # ---------------- warm-up phase --------------------------------
      while (simC < n_w) {
        dt <- min(tArr[iArr], next_dep())
        t  <- t + dt
        cCum <- cCum + dt * sys
        if (sys > s) dCum <- dCum + dt * (sys - s)
        if (length(busy)) busy <- busy - dt
        tArr[iArr] <- tArr[iArr] - dt

        if (tArr[iArr] <= 1e-12) {            # arrival
          if (sys < s + K) {                  # accepted
            simC <- simC + 1L
            sys  <- sys + 1L
            iArr <- iArr + 1L
            if (length(busy) < s) {           # starts service
              busy <- c(busy, tServ[iServ]); iServ <- iServ + 1L
            }
          } else {                            # lost
            lost <- lost + 1L; iArr <- iArr + 1L
          }
        } else {                              # departure
          idx  <- which.min(busy)
          busy <- busy[-idx]
          sys  <- sys - 1L
          if (sys >= s) {                     # pull from queue
            busy <- c(busy, tServ[iServ]); iServ <- iServ + 1L
          }
        }

        if (self$historic) {
          l  <- cCum / t
          lq <- dCum / t
          w  <- if (simC) cCum / simC else NA
          wq <- if (simC) dCum / simC else NA
          hist[simC, ] <- c(l,lq,w,wq,sys,(l-lq)/s,t)
        }
      }

      tWarm <- t
      cCum <- dCum <- 0
      simC <- 0L
      pnVec <- numeric(n_m)

      # ---------------- main phase -----------------------------------
      while (simC < n_m) {
        dt <- min(tArr[iArr], next_dep())
        t  <- t + dt
        cCum <- cCum + dt * sys
        if (sys > s) dCum <- dCum + dt * (sys - s)
        pnVec[min(sys, n_m) + 1L] <- pnVec[min(sys, n_m) + 1L] + dt
        if (length(busy)) busy <- busy - dt
        tArr[iArr] <- tArr[iArr] - dt

        if (tArr[iArr] <= 1e-12) {            # arrival
          if (sys < s + K) {                  # accepted
            simC <- simC + 1L
            sys  <- sys + 1L
            iArr <- iArr + 1L
            if (length(busy) < s) {
              busy <- c(busy, tServ[iServ]); iServ <- iServ + 1L
            }
          } else {
            lost <- lost + 1L; iArr <- iArr + 1L
          }
        } else {                              # departure
          idx  <- which.min(busy)
          busy <- busy[-idx]
          sys  <- sys - 1L
          if (sys >= s) {
            busy <- c(busy, tServ[iServ]); iServ <- iServ + 1L
          }
        }

        if (self$historic && simC > 0L) {
          l  <- cCum / (t - tWarm)
          lq <- dCum / (t - tWarm)
          w  <- cCum / simC
          wq <- dCum / simC
          hist[simC + n_w, ] <- c(l,lq,w,wq,sys,(l-lq)/s,t)
        }
      }

      # ---------------- statistics -----------------------------------
      L  <- cCum / (t - tWarm)
      Lq <- dCum / (t - tWarm)
      W  <- cCum / n_m
      Wq <- dCum / n_m
      rho <- (L - Lq) / s
      eff <- W / (W - Wq)

      out <- list(l = L, lq = Lq, w = W, wq = Wq,
                  rho = rho, eff = eff,
                  pn  = pnVec / (t - tWarm))
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Convenience wrapper with replication / parallel                    ##
########################################################################

#' Simulate a G/G/s/K queue
#'
#' @inheritParams GGSKSim
#' @param nsim  integer, number of independent replications.
#' @param nproc integer, CPU workers (1 = sequential).
#' @return A `GGSKSim` object or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
GGSK <- function(arrivalDistribution = Exp(3),
                 serviceDistribution = Exp(6),
                 servers       = 2L,
                 K             = 3L,
                 staClients    = 100L,
                 nClients      = 1000L,
                 historic      = FALSE,
                 nsim          = 10L,
                 nproc         = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  sim_fun <- function(arrivalDistribution,
                      serviceDistribution,
                      servers,
                      K,
                      staClients,
                      nClients,
                      historic) {
    GGSKSim$new(arrivalDistribution,
                serviceDistribution,
                servers,
                K,
                staClients,
                nClients,
                historic)
  }

  params <- list(arrivalDistribution,
                 serviceDistribution,
                 servers,
                 K,
                 staClients,
                 nClients,
                 historic)

  res <- ParallelizeSimulations(sim_fun, params, nsim, nproc)
  if (inherits(res, "GGSKSim")) res else combineSimulations(res)
}

########################################################################
## GG1HSim – simulated G/G/1/∞/H queue (finite-source)                ##
########################################################################

#' GG1_Inf_HSIMHSim – simulated G/G/1/Inf/H queue (R6)
#'
#' @description Discrete-event simulation of a single-server queue fed
#'   by a *finite* population of `H` sources.  At any moment each source
#'   is either **in** the system (being served or waiting) or **outside**
#'   and generating its own inter-arrival time.  The total population is
#'   constant and no arrivals are lost.
#'
#' @param arrivalDistribution,serviceDistribution Objects from package
#'   **distr** with the inter-arrival and service-time laws.
#' @param H          Integer >= 1, size of the customer population.
#' @param staClients Warm-up customers discarded from statistics.
#' @param nClients   Customers collected for statistics.
#' @param historic   Logical, store the whole trajectory?
#' @export
GG1_Inf_HSIMHSim <- R6::R6Class(
  "GG1_Inf_HSIMHSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    H                   = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          H           = 5L,
                          staClients  = 100L,
                          nClients    = 1000L,
                          historic    = FALSE) {

      if (!belong(arrivalDistribution, "distr"))
        stop("'arrivalDistribution' must be a distr object", call. = FALSE)
      if (!belong(serviceDistribution, "distr"))
        stop("'serviceDistribution' must be a distr object", call. = FALSE)
      if (!is.numeric(H) || H < 1)
        stop("'H' must be >= 1", call. = FALSE)
      if (!is.numeric(staClients) || staClients < 0)
        stop("'staClients' must be >= 0", call. = FALSE)
      if (!is.numeric(nClients)  || nClients  <= 0)
        stop("'nClients' must be > 0", call. = FALSE)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$H                   <- as.integer(H)
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GG1_Inf_HSIMHSim  (G/G/1/Inf/", self$H, ")\n", sep = "")
      with(self$out, {
        cat(sprintf("L   = %.5f | Lq  = %.5f\n", l,  lq))
        cat(sprintf("W   = %.5f | Wq  = %.5f\n", w,  wq))
        cat(sprintf("rho = %.5f | Eff = %.5f\n", rho, eff))
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      H   <- self$H
      n_w <- self$staClients
      n_m <- self$nClients
      n_tot <- (n_w + n_m) * 2L                            # cushion

      tArr  <- r(self$arrivalDistribution)(n_tot)
      tServ <- r(self$serviceDistribution)(n_tot)
      if (any(tArr < 0) || any(tServ < 0))
        stop("negative times produced by distributions")

      ## initial arrival clocks for each source -----------------------
      nextArr <- tArr[1:H]
      iArr    <- H + 1L
      iServ   <- 1L

      sys  <- 0L; simC <- 0L; t <- 0
      nServed <- 0L
      cCum <- dCum <- 0
      pnVec <- numeric(n_m)

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = n_w + n_m, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      ### ---------- helper for next departure ------------------------
      next_dep_time <- function() if (sys) tServ[iServ] else Inf

      ### ---------- warm-up phase ------------------------------------
      while (simC < n_w) {
        dt <- min(min(nextArr, na.rm = TRUE), next_dep_time())
        t  <- t + dt
        if (sys) {
          cCum <- cCum + dt * sys
          dCum <- dCum + dt * (sys - 1)
          tServ[iServ] <- tServ[iServ] - dt
        }
        nextArr <- nextArr - dt

        if (abs(next_dep_time() - 0) < 1e-12) {  # departure
          sys <- sys - 1L
          if (sys) iServ <- iServ + 1L
        } else {                                 # arrival
          idx <- which.min(nextArr)
          simC <- simC + 1L
          if (sys == 0) {
            sys <- 1L
            nextArr[idx] <- NA
          } else {
            sys  <- sys + 1L
            nextArr[idx] <- NA
          }
        }

        # start new inter-arrival for the source just freed
        if (any(is.na(nextArr))) {
          nextArr[is.na(nextArr)][1] <- tArr[iArr]; iArr <- iArr + 1L
        }

        if (self$historic) {
          l  <- cCum / t
          lq <- dCum / t
          w  <- cCum / simC
          wq <- dCum / simC
          hist[simC, ] <- c(l,lq,w,wq,sys,l-lq,t)
        }
      }

      tWarm <- t
      cCum <- dCum <- 0
      simC <- 0L

      ### ---------- main statistics phase ----------------------------
      while (simC < n_m) {
        dt <- min(min(nextArr, na.rm = TRUE), next_dep_time())
        t  <- t + dt
        if (sys) {
          cCum <- cCum + dt * sys
          dCum <- dCum + dt * (sys - 1)
          tServ[iServ] <- tServ[iServ] - dt
        }
        nextArr <- nextArr - dt
        pnVec[min(sys, n_m) + 1L] <- pnVec[min(sys, n_m) + 1L] + dt

        if (abs(next_dep_time() - 0) < 1e-12) {  # departure
          sys <- sys - 1L
          if (sys) iServ <- iServ + 1L
        } else {                                 # arrival
          idx <- which.min(nextArr)
          simC <- simC + 1L
          if (sys == 0) {
            sys <- 1L
            nextArr[idx] <- NA
          } else {
            sys  <- sys + 1L
            nextArr[idx] <- NA
          }
        }

        if (any(is.na(nextArr))) {
          nextArr[is.na(nextArr)][1] <- tArr[iArr]; iArr <- iArr + 1L
        }

        if (self$historic && simC > 0L) {
          l  <- cCum / (t - tWarm)
          lq <- dCum / (t - tWarm)
          w  <- cCum / simC
          wq <- dCum / simC
          hist[simC + n_w, ] <- c(l,lq,w,wq,sys,l-lq,t)
        }
      }

      ### ---------- final statistics ---------------------------------
      L  <- cCum / (t - tWarm)
      Lq <- dCum / (t - tWarm)
      W  <- cCum / n_m
      Wq <- dCum / n_m
      rho <- L - Lq
      eff <- W / (W - Wq)

      out <- list(l = L, lq = Lq, w = W, wq = Wq,
                  rho = rho, eff = eff,
                  pn  = pnVec / (t - tWarm))
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Functional wrapper with replication / parallel                     ##
########################################################################

#' Simulate a G/G/1/Inf/H queue
#'
#' @inheritParams GG1_Inf_HSIMHSim
#' @param nsim  Integer, number of independent replications.
#' @param nproc Integer, CPU cores (1 = sequential).
#' @return A `GG1_Inf_HSIMHSim` object, or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
GG1_Inf_HSIMH <- function(arrivalDistribution = Exp(3),
                 serviceDistribution = Exp(6),
                 H             = 5L,
                 staClients    = 100L,
                 nClients      = 1000L,
                 historic      = FALSE,
                 nsim          = 10L,
                 nproc         = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  sim_fun <- function(arrivalDistribution,
                      serviceDistribution,
                      H,
                      staClients,
                      nClients,
                      historic) {
    GG1_Inf_HSIMHSim$new(arrivalDistribution,
                serviceDistribution,
                H,
                staClients,
                nClients,
                historic)
  }

  params <- list(arrivalDistribution,
                 serviceDistribution,
                 H,
                 staClients,
                 nClients,
                 historic)

  res <- ParallelizeSimulations(sim_fun, params, nsim, nproc)
  if (inherits(res, "GG1_Inf_HSIMHSim")) res else combineSimulations(res)
}

#' GGS_Inf_HSim – simulated G/G/s/Inf/H queue (R6)
#'
#' @description Discrete-event simulation of a multi-server queue
#'   (`s` identical servers) fed by a *finite* population of `H`
#'   sources.  Each source alternates between **inside** the system
#'   (being served or waiting) and **outside** where it generates its
#'   own inter-arrival time.  The total population is constant and no
#'   arrivals are lost.
#'
#' @param arrivalDistribution,serviceDistribution Objects from **distr**
#'   giving the inter-arrival and service-time laws.
#' @param servers Integer >= 1, number of servers (*s*).
#' @param H       Integer >= 1, customer population size.
#' @param staClients Warm-up customers discarded from statistics.
#' @param nClients   Customers collected for statistics.
#' @param historic   Logical, store the whole trajectory?
#' @export
GGS_Inf_HSim <- R6::R6Class(
  "GGS_Inf_HSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    servers             = NULL,
    H                   = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          servers      = 3L,
                          H            = 5L,
                          staClients   = 100L,
                          nClients     = 1000L,
                          historic     = FALSE) {

      stopifnot(belong(arrivalDistribution, "distr"),
                belong(serviceDistribution, "distr"),
                servers  >= 1, H >= 1,
                staClients >= 0, nClients > 0)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$servers             <- as.integer(servers)
      self$H                   <- as.integer(H)
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GGS_Inf_HSim  (G/G/", self$servers,
          "/Inf/", self$H, ")\n", sep = "")
      with(self$out, {
        cat(sprintf("L   = %.5f | Lq  = %.5f\n", l,  lq))
        cat(sprintf("W   = %.5f | Wq  = %.5f\n", w,  wq))
        cat(sprintf("rho = %.5f | Eff = %.5f\n", rho, eff))
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      s  <- self$servers
      H  <- self$H
      nw <- self$staClients
      nm <- self$nClients
      tot <- (nw + nm) * 2L                       # safety cushion

      tArr  <- r(self$arrivalDistribution)(tot)
      tServ <- r(self$serviceDistribution)(tot)
      if (any(tArr < 0) || any(tServ < 0))
        stop("negative times produced by distributions")

      ## arrival and service clocks -----------------------------------
      nextArr <- tArr[1:H]        # one clock per external source
      bServ   <- rep(NA_real_, s) # remaining service times per server
      iArr    <- H + 1L           # next inter-arrival index
      iServ   <- 1L               # next service-time index

      sys  <- 0L; simC <- 0L; t <- 0
      cCum <- dCum <- 0
      pnVec <- numeric(nm)

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = nw + nm, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      next_dep <- function() if (all(is.na(bServ))) Inf else min(bServ, na.rm = TRUE)

      advance <- function(dt) {
        t  <<- t + dt
        if (sys) {
          cCum <<- cCum + dt * sys
          dCum <<- dCum + dt * max(sys - s, 0)
          bServ[!is.na(bServ)] <<- bServ[!is.na(bServ)] - dt
        }
        nextArr[!is.na(nextArr)] <<- nextArr[!is.na(nextArr)] - dt
      }

      record <- function(row) {
        L  <- cCum / (t - tWarm)
        Lq <- dCum / (t - tWarm)
        W  <- cCum / simC
        Wq <- dCum / simC
        hist[row, ] <<- c(L,Lq,W,Wq,sys,(L-Lq)/s,t)
      }

      ## ----------------- warm-up phase ------------------------------
      while (simC < nw) {
        dt <- min(min(nextArr, na.rm = TRUE), next_dep())
        advance(dt)

        if (abs(next_dep()) < 1e-12) {            # departure
          sys <- sys - 1L
          bServ[bServ <= 1e-12] <- NA
          if (sys >= s) {
            free <- which(is.na(bServ))[1]
            bServ[free] <- tServ[iServ]; iServ <- iServ + 1L
          }
        } else {                                  # arrival
          idx <- which.min(nextArr)
          simC <- simC + 1L
          if (sys < s) {
            free <- which(is.na(bServ))[1]
            bServ[free] <- tServ[iServ]; iServ <- iServ + 1L
          }
          sys <- sys + 1L
          nextArr[idx] <- NA
        }

        if (any(is.na(nextArr))) {
          nextArr[is.na(nextArr)][1] <- tArr[iArr]; iArr <- iArr + 1L
        }

        if (self$historic) record(simC)
      }

      tWarm <- t
      cCum <- dCum <- 0
      simC <- 0L

      ## ----------------- main phase ---------------------------------
      while (simC < nm) {
        dt <- min(min(nextArr, na.rm = TRUE), next_dep())
        advance(dt)
        pnVec[min(sys, nm) + 1L] <- pnVec[min(sys, nm) + 1L] + dt

        if (abs(next_dep()) < 1e-12) {            # departure
          sys <- sys - 1L
          bServ[bServ <= 1e-12] <- NA
          if (sys >= s) {
            free <- which(is.na(bServ))[1]
            bServ[free] <- tServ[iServ]; iServ <- iServ + 1L
          }
        } else {                                  # arrival
          idx <- which.min(nextArr)
          simC <- simC + 1L
          if (sys < s) {
            free <- which(is.na(bServ))[1]
            bServ[free] <- tServ[iServ]; iServ <- iServ + 1L
          }
          sys <- sys + 1L
          nextArr[idx] <- NA
        }

        if (any(is.na(nextArr))) {
          nextArr[is.na(nextArr)][1] <- tArr[iArr]; iArr <- iArr + 1L
        }

        if (self$historic && simC > 0L) record(simC + nw)
      }

      ## ----------------- final stats --------------------------------
      L   <- cCum / (t - tWarm)
      Lq  <- dCum / (t - tWarm)
      W   <- cCum / nm
      Wq  <- dCum / nm
      rho <- (L - Lq) / s
      eff <- W / (W - Wq)

      out <- list(l = L, lq = Lq, w = W, wq = Wq,
                  rho = rho, eff = eff,
                  pn  = pnVec / (t - tWarm))
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Functional constructor with replication / parallel                 ##
########################################################################

#' Simulate a G/G/s/Inf/H queue
#'
#' Wrapper around `GGS_Inf_HSim$new()` plus optional parallel
#' replication.
#'
#' @inheritParams GGS_Inf_HSim
#' @param nsim  Integer, number of independent replications.
#' @param nproc Integer, CPU cores (1 = sequential).
#' @return A `GGS_Inf_HSim` object, or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
GGS_Inf_H <- function(arrivalDistribution = Exp(3),
                          serviceDistribution = Exp(6),
                          servers      = 3L,
                          H            = 5L,
                          staClients   = 100L,
                          nClients     = 1000L,
                          historic     = FALSE,
                          nsim         = 10L,
                          nproc        = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  make_one <- function(arrivalDistribution,
                       serviceDistribution,
                       servers,
                       H,
                       staClients,
                       nClients,
                       historic) {
    GGS_Inf_HSim$new(arrivalDistribution,
                         serviceDistribution,
                         servers,
                         H,
                         staClients,
                         nClients,
                         historic)
  }

  parms <- list(arrivalDistribution, serviceDistribution,
                servers, H, staClients, nClients, historic)

  res <- ParallelizeSimulations(make_one, parms, nsim, nproc)
  if (inherits(res, "GGS_Inf_HSim")) res else combineSimulations(res)
}

#' GGS_Inf_HYSim – simulated G/G/s/Inf/H queue with Y replacements (R6)
#'
#' @description Multi-server queue (`s` identical servers) fed by a
#'   *finite* population of `H` sources **plus** a replacement rule:
#'   when the number of customers *inside* the system is `<= Y`
#'   the source that just left is immediately replaced by a new one
#'   (fresh inter-arrival time is generated). For `Y >= H` the
#'   behaviour degenerates to the plain `G/G/s/Inf/H` case.
#'
#' @param arrivalDistribution,serviceDistribution Objects from **distr**
#'   defining the inter-arrival and service-time laws.
#' @param servers Integer >= 1, number of servers (*s*).
#' @param H       Integer >= 1, finite customer population.
#' @param Y       Integer >= 1, replacement threshold.
#' @param staClients Warm-up customers discarded from statistics.
#' @param nClients   Customers collected for statistics.
#' @param historic   Logical, store the whole trajectory?
#' @export
GGS_Inf_HYSim <- R6::R6Class(
  "GGS_Inf_HYSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    servers             = NULL,
    H                   = NULL,
    Y                   = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          servers      = 3L,
                          H            = 5L,
                          Y            = 3L,
                          staClients   = 100L,
                          nClients     = 1000L,
                          historic     = FALSE) {

      stopifnot(belong(arrivalDistribution, "distr"),
                belong(serviceDistribution, "distr"),
                servers  >= 1, H >= 1, Y >= 1,
                staClients >= 0, nClients > 0)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$servers             <- as.integer(servers)
      self$H                   <- as.integer(H)
      self$Y                   <- as.integer(Y)
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GGS_Inf_HYSim  (G/G/", self$servers,
          "/Inf/", self$H, ", Y = ", self$Y, ")\n", sep = "")
      with(self$out, {
        cat(sprintf("L   = %.5f | Lq  = %.5f\n", l,  lq))
        cat(sprintf("W   = %.5f | Wq  = %.5f\n", w,  wq))
        cat(sprintf("rho = %.5f | Eff = %.5f\n", rho, eff))
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      s  <- self$servers
      H  <- self$H
      Yt <- self$Y
      nw <- self$staClients
      nm <- self$nClients
      tot <- (nw + nm) * 2L                       # cushion

      tArr  <- r(self$arrivalDistribution)(tot)
      tServ <- r(self$serviceDistribution)(tot)
      if (any(tArr < 0) || any(tServ < 0))
        stop("negative times produced by distributions")

      ## arrival clocks per external source ---------------------------
      nextArr <- tArr[1:H]
      iArr    <- H + 1L
      ## service clocks per busy server -------------------------------
      bServ   <- rep(NA_real_, s)
      iServ   <- 1L

      sys  <- 0L; simC <- 0L; t <- 0
      cCum <- dCum <- 0
      pnVec <- numeric(nm)

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = nw + nm, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      next_dep <- function() if (all(is.na(bServ))) Inf else min(bServ, na.rm = TRUE)

      advance <- function(dt) {
        t  <<- t + dt
        if (sys) {
          cCum <<- cCum + dt * sys
          dCum <<- dCum + dt * max(sys - s, 0)
          bServ[!is.na(bServ)] <<- bServ[!is.na(bServ)] - dt
        }
        nextArr[!is.na(nextArr)] <<- nextArr[!is.na(nextArr)] - dt
      }

      maybe_new_source <- function() {
        if (any(is.na(nextArr))) {
          nextArr[is.na(nextArr)][1] <<- tArr[iArr]; iArr <<- iArr + 1L
        }
      }

      record <- function(row) {
        L  <- cCum / (t - tWarm)
        Lq <- dCum / (t - tWarm)
        W  <- cCum / simC
        Wq <- dCum / simC
        hist[row, ] <<- c(L,Lq,W,Wq,sys,(L-Lq)/s,t)
      }

      process_arrival <- function(idx) {
        simC <<- simC + 1L
        if (sys < s) {
          free <- which(is.na(bServ))[1]
          bServ[free] <<- tServ[iServ]; iServ <<- iServ + 1L
        }
        sys <<- sys + 1L
        nextArr[idx] <<- NA
        if (sys <= Yt) maybe_new_source()
      }

      process_departure <- function() {
        sys <<- sys - 1L
        bServ[bServ <= 1e-12] <<- NA
        if (sys >= s) {
          free <- which(is.na(bServ))[1]
          bServ[free] <<- tServ[iServ]; iServ <<- iServ + 1L
        }
        if (sys >= Yt) maybe_new_source()
      }

      ## ----------------- warm-up phase ------------------------------
      while (simC < nw) {
        dt <- min(min(nextArr, na.rm = TRUE), next_dep())
        advance(dt)

        if (abs(next_dep()) < 1e-12) {
          process_departure()
        } else {
          idx <- which.min(nextArr)
          process_arrival(idx)
        }

        if (self$historic) record(simC)
      }

      tWarm <- t
      cCum <- dCum <- 0
      simC <- 0L

      ## ----------------- statistics phase ---------------------------
      while (simC < nm) {
        dt <- min(min(nextArr, na.rm = TRUE), next_dep())
        advance(dt)
        pnVec[min(sys, nm) + 1L] <- pnVec[min(sys, nm) + 1L] + dt

        if (abs(next_dep()) < 1e-12) {
          process_departure()
        } else {
          idx <- which.min(nextArr)
          process_arrival(idx)
        }

        if (self$historic && simC > 0L) record(simC + nw)
      }

      ## ----------------- final stats -------------------------------
      L   <- cCum / (t - tWarm)
      Lq  <- dCum / (t - tWarm)
      W   <- cCum / nm
      Wq  <- dCum / nm
      rho <- (L - Lq) / s
      eff <- W / (W - Wq)

      out <- list(l = L, lq = Lq, w = W, wq = Wq,
                  rho = rho, eff = eff,
                  pn  = pnVec / (t - tWarm))
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Functional constructor with replication / parallel                 ##
########################################################################

#' Simulate a G/G/s/Inf/H/Y queue
#'
#' Wrapper around `GGS_Inf_HYSim$new()` plus optional parallel
#' replication.
#'
#' @inheritParams GGS_Inf_HYSim
#' @param nsim  Integer, number of independent replications.
#' @param nproc Integer, CPU cores (1 = sequential).
#' @return A `GGS_Inf_HYSim` object, or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
GGS_Inf_HY <- function(arrivalDistribution = Exp(3),
                       serviceDistribution = Exp(6),
                       servers      = 3L,
                       H            = 5L,
                       Y            = 3L,
                       staClients   = 100L,
                       nClients     = 1000L,
                       historic     = FALSE,
                       nsim         = 10L,
                       nproc        = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  make_one <- function(arrivalDistribution,
                       serviceDistribution,
                       servers,
                       H,
                       Y,
                       staClients,
                       nClients,
                       historic) {
    GGS_Inf_HYSim$new(arrivalDistribution,
                      serviceDistribution,
                      servers,
                      H,
                      Y,
                      staClients,
                      nClients,
                      historic)
  }

  parms <- list(arrivalDistribution, serviceDistribution,
                servers, H, Y, staClients, nClients, historic)

  res <- ParallelizeSimulations(make_one, parms, nsim, nproc)
  if (inherits(res, "GGS_Inf_HYSim")) res else combineSimulations(res)
}

#' GGInfSim – simulated G/G/Inf queue (R6)
#'
#' @description Discrete-event simulation of a queueing system with
#'   unlimited parallel servers: every arrival starts service
#'   immediately, therefore there is **no queue** and the waiting-time
#'   in queue is always 0.  The class keeps exactly the same interface
#'   used by the other simulation models in this package.
#'
#' @param arrivalDistribution,serviceDistribution  Distributions from
#'   package **distr** describing inter-arrival and service times.
#' @param staClients  Warm-up customers discarded from statistics.
#' @param nClients    Customers included in statistics.
#' @param historic    Logical, store the whole trajectory
#' @export
GGInfSim <- R6::R6Class(
  "GGInfSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    staClients          = NULL,
    nClients            = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          staClients = 100L,
                          nClients   = 1000L,
                          historic   = FALSE) {

      stopifnot(belong(arrivalDistribution, "distr"),
                belong(serviceDistribution, "distr"),
                staClients >= 0, nClients > 0)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$staClients          <- as.integer(staClients)
      self$nClients            <- as.integer(nClients)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: GGInfSim  (G/G/Inf)\n")
      with(self$out, {
        cat(sprintf("L   = %.5f\n",  l))
        cat(sprintf("W   = %.5f\n",  w))
        cat("Lq  = 0.00000 | Wq = 0.00000 | rho = 0\n")
        cat("Eff = 1.00000 (no queue)\n")
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      n_w <- self$staClients
      n_m <- self$nClients
      n_tot <- n_w + n_m                           # total arrivals

      tArr  <- r(self$arrivalDistribution)(n_tot)
      tServ <- r(self$serviceDistribution)(n_tot)
      if (any(tArr < 0) || any(tServ < 0))
        stop("negative times produced by distributions")

      next_arr <- tArr[1]; iArr <- 2L
      iServ <- 1L
      serv_clock <- numeric()                      # remaining service times

      t <- 0; simC <- 0L
      cCum <- 0
      pnVec <- numeric(n_m)

      if (self$historic) {
        hist <- matrix(NA_real_, nrow = n_tot, ncol = 7,
                       dimnames = list(NULL,
                                       c("L","Lq","W","Wq","Clients","Intensity","tClient")))
      }

      advance <- function(dt) {
        t <<- t + dt
        if (length(serv_clock))
          serv_clock <<- serv_clock - dt
        cCum <<- cCum + dt * length(serv_clock)
        next_arr <<- next_arr - dt
      }

      process_arrival <- function() {
        simC <<- simC + 1L
        serv_clock <<- c(serv_clock, tServ[iServ]); iServ <<- iServ + 1L
        if (iArr <= n_tot) {                       # schedule next arrival
          next_arr <<- tArr[iArr]; iArr <<- iArr + 1L
        } else {
          next_arr <<- Inf
        }
      }

      process_departure <- function() {
        idx <- which.min(serv_clock)
        serv_clock <<- serv_clock[-idx]
      }

      ## ---------------- warm-up ------------------------------------
      while (simC < n_w) {
        next_dep <- if (length(serv_clock)) min(serv_clock) else Inf
        dt <- min(next_arr, next_dep)
        advance(dt)

        if (next_arr <= next_dep) {
          process_arrival()
        } else {
          process_departure()
        }

        if (self$historic)
          hist[simC, ] <- c(length(serv_clock), 0, cCum/simC, 0,
                            length(serv_clock), 0, t)
      }

      tWarm <- t;  cCum <- 0; simC <- 0L

      ## ---------------- statistics phase ---------------------------
      while (simC < n_m) {
        next_dep <- if (length(serv_clock)) min(serv_clock) else Inf
        dt <- min(next_arr, next_dep)
        advance(dt)
        pnIdx <- min(length(serv_clock), n_m) + 1L
        pnVec[pnIdx] <- pnVec[pnIdx] + dt

        if (next_arr <= next_dep) {
          process_arrival()
        } else {
          process_departure()
        }

        if (self$historic && simC > 0L)
          hist[simC + n_w, ] <-
          c(length(serv_clock), 0, cCum/simC, 0,
            length(serv_clock), 0, t)
      }

      L  <- cCum / (t - tWarm)
      W  <- cCum / n_m

      out <- list(l = L, lq = 0, w = W, wq = 0,
                  rho = 0, eff = 1,
                  pn = pnVec / (t - tWarm))
      if (self$historic) out$historic <- hist
      self$out <- out
    }
  )
)

########################################################################
## Functional wrapper with replication / parallel                     ##
########################################################################

#' Simulate a G/G/Inf queue
#'
#' Wrapper around `GGInfSim$new()` plus optional parallel replications.
#'
#' @inheritParams GGInfSim
#' @param nsim  Integer, number of independent replications.
#' @param nproc Integer, CPU cores (1 = sequential).
#' @return A `GGInfSim` object, or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
GGInf <- function(arrivalDistribution = Exp(3),
                  serviceDistribution = Exp(6),
                  staClients = 100L,
                  nClients   = 1000L,
                  historic   = FALSE,
                  nsim       = 10L,
                  nproc      = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  make_one <- function(arrivalDistribution,
                       serviceDistribution,
                       staClients,
                       nClients,
                       historic) {
    GGInfSim$new(arrivalDistribution,
                 serviceDistribution,
                 staClients,
                 nClients,
                 historic)
  }

  parms <- list(arrivalDistribution, serviceDistribution,
                staClients, nClients, historic)

  res <- ParallelizeSimulations(make_one, parms, nsim, nproc)
  if (inherits(res, "GGInfSim")) res else combineSimulations(res)
}

#' ClosedNetSim – simulated closed Jackson-type network (R6)
#'
#' @description Pure simulation (event-driven) of a *closed* queueing
#'   network with a fixed population of `nClients` customers.  Each node
#'   is a `G/G/s` queue.  After completing service a customer is routed
#'   to the next node according to the user-supplied routing matrix `p`
#'   (rows must sum to 1).  The algorithm is a line-by-line port of the
#'   original S3 function, wrapped now in an R6 class.
#'
#' @param serviceDistribution List of `distr` objects (one per node).
#' @param s      Integer vector with the number of servers at each node.
#' @param p      Routing matrix (`length(s)` × `length(s)`), row-stochastic.
#' @param nClients  Total number of circulating customers `N`.
#' @param staClients  Warm-up completions discarded from statistics.
#' @param transitions Number of completed services counted for statistics.
#' @param historic  Logical, store the full trajectory?
#' @export
ClosedNetSim <- R6::R6Class(
  "ClosedNetSim",
  public = list(
    serviceDistribution = NULL,
    servers             = NULL,
    routing             = NULL,
    nClients            = NULL,
    staClients          = NULL,
    transitions         = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(serviceDistribution,
                          s,
                          p,
                          nClients    = 3L,
                          staClients  = 100L,
                          transitions = 1000L,
                          historic    = FALSE) {

      stopifnot(is.list(serviceDistribution),
                all(vapply(serviceDistribution, inherit, logical(1), what = "distr")),
                is.numeric(s), length(s)   == length(serviceDistribution),
                is.matrix(p),  nrow(p)     == length(s),
                all(abs(rowSums(p) - 1) < 1e-10),
                nClients   > 0,
                staClients >= 0,
                transitions > 0)

      self$serviceDistribution <- serviceDistribution
      self$servers             <- as.integer(s)
      self$routing             <- p
      self$nClients            <- as.integer(nClients)
      self$staClients          <- as.integer(staClients)
      self$transitions         <- as.integer(transitions)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: ClosedNetSim  (closed Jackson network)\n")
      with(self$out, {
        tab <- data.frame(L  = l,
                          Lq = lq,
                          W  = w,
                          Wq = wq,
                          row.names = paste0("Node", seq_along(l)))
        print(tab)
        cat("Total Lq =", lqt, "\n")
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      k  <- length(self$servers)
      s  <- self$servers
      N  <- self$nClients
      nw <- self$staClients
      nm <- self$transitions
      maxS <- max(s)

      ## --- pre-generate service times -------------------------------
      genTimes <- function(d) r(d)(nw + nm + 5L)  # cushion
      tServMat <- sapply(self$serviceDistribution, genTimes)

      ## --- initial allocation of customers --------------------------
      inNode <- floor(N * (s ^ -1) / sum(s ^ -1))
      inNode[ seq_len(N - sum(inNode)) ] <- inNode[ seq_len(N - sum(inNode)) ] + 1

      ## service-clock matrix (rows = servers, cols = nodes) ----------
      b <- matrix(NA_real_, nrow = maxS, ncol = k)
      iServ <- rep(1L, k)

      for (j in seq_len(k)) {
        m <- min(s[j], inNode[j])
        if (m) {
          idx <- seq_len(m)
          b[idx, j] <- tServMat[idx, j]
          iServ[j]  <- m + 1L
        }
      }

      ## --- state variables ------------------------------------------
      sys    <- inNode               # customers per node
      cCum   <- dCum <- cron <- numeric(k)
      comp   <- 0L                    # completed services counter
      randU  <- runif(nw + nm + 5L)   # routing RNG

      if (self$historic) {
        hist <- array(NA_real_, dim = c(k, 5, nw + nm),
                      dimnames = list(
                        paste0("Node", seq_len(k)),
                        c("L","Lq","W","Wq","tClient"),
                        NULL))
      }

      advance_one <- function(step) {
        ## locate node and service finishing next ---------------------
        flatIdx <- which.min(b)
        dt      <- b[flatIdx]
        node    <- ((flatIdx - 1L) %/% maxS) + 1L
        b       <<- b - dt

        ## update integrals ------------------------------------------
        cron[node]  <<- cron[node] + dt
        for (j in seq_len(k)) {
          cCum[j] <<- cCum[j] + dt * sys[j]
          if (sys[j] > s[j])
            dCum[j] <<- dCum[j] + dt * (sys[j] - s[j])
        }

        ## departure from node ---------------------------------------
        sys[node] <<- sys[node] - 1L
        if (sys[node] < s[node]) {
          b[flatIdx] <<- NA
        } else {
          b[flatIdx] <<- tServMat[iServ[node], node]
          iServ[node] <<- iServ[node] + 1L
        }

        ## routing ---------------------------------------------------
        target <- findInterval(randU[step], vec = cumsum(self$routing[node, ])) + 1L
        sys[target] <<- sys[target] + 1L

        if (sys[target] <= s[target]) {
          srvIdx <- which(is.na(b[, target]))[1]
          b[srvIdx, target] <<- tServMat[iServ[target], target]
          iServ[target] <<- iServ[target] + 1L
        }
      }

      ## ---------------- warm-up ------------------------------------
      while (comp < nw) {
        comp <- comp + 1L
        advance_one(comp)
        if (self$historic)
          hist[ , , comp] <- private$snapshot_state(comp, cCum, dCum, cron)
      }

      ## reset counters for main stats --------------------------------
      cCum[] <- dCum[] <- cron[] <- 0
      comp   <- 0L

      ## ---------------- main statistics loop ------------------------
      while (comp < nm) {
        comp <- comp + 1L
        advance_one(comp + nw)
        if (self$historic)
          hist[ , , comp + nw] <- private$snapshot_state(comp, cCum, dCum, cron)
      }

      ## ---------------- final metrics -------------------------------
      L   <- cCum / sum(cron)
      Lq  <- dCum / sum(cron)
      W   <- cCum / inNode
      Wq  <- dCum / inNode
      rho <- L - Lq
      eff <- W / (W - Wq)

      self$out <- list(l   = L,
                       lq  = Lq,
                       lqt = sum(Lq),
                       w   = W,
                       wq  = Wq,
                       pn  = NULL,   # omitted – rarely useful for large nets
                       rho = rho,
                       eff = eff)
      if (self$historic) self$out$historic <- hist
    },

    snapshot_state = function(step, cCum, dCum, cron) {
      L  <- cCum / sum(cron)
      Lq <- dCum / sum(cron)
      W  <- ifelse(step == 0, 0, cCum / step)
      Wq <- ifelse(step == 0, 0, dCum / step)
      array(c(L, Lq, W, Wq, cron),
            dim = c(length(L), 5L))
    }
  )
)

########################################################################
## Functional wrapper with replication / parallel                     ##
########################################################################

#' Simulate a closed Jackson-type network
#'
#' Wrapper around `ClosedNetSim$new()` plus optional parallel replicas.
#'
#' @inheritParams ClosedNetSim
#' @param nsim  Integer, number of independent replications.
#' @param nproc Integer, CPU cores (1 = sequential).
#' @return A `ClosedNetSim` object, or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
ClosedNet <- function(serviceDistribution,
                      s,
                      p,
                      nClients    = 3L,
                      staClients  = 100L,
                      transitions = 1000L,
                      historic    = FALSE,
                      nsim        = 10L,
                      nproc       = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  make_one <- function(serviceDistribution, s, p,
                       nClients, staClients, transitions, historic) {
    ClosedNetSim$new(serviceDistribution,
                     s,
                     p,
                     nClients,
                     staClients,
                     transitions,
                     historic)
  }

  args <- list(serviceDistribution, s, p,
               nClients, staClients, transitions, historic)

  res <- ParallelizeSimulations(make_one, args, nsim, nproc)
  if (inherits(res, "ClosedNetSim")) res else combineSimulations(res)
}

#' OpenNetSim – simulated open Jackson-type network (R6)
#'
#' @description Discrete-event simulation of an *open* queueing network
#'   with external arrivals at one or more nodes and probabilistic
#'   routing between nodes.  Each node behaves as a `G/G/s` queue.  The
#'   internal logic is delegated to the legacy helper
#'   `OpenNetwork_secuential()` so the numerical results remain
#'   identical to the original S3 code, but all outputs are exposed
#'   through the field `out` of the R6 object.
#'
#' @param arrivalDistribution List of `distr` objects (or `no_distr()`)
#'   giving the external arrival law for every node.
#' @param serviceDistribution List of `distr` objects with service-time
#'   laws.
#' @param s      Integer vector, servers per node.
#' @param p      Routing matrix; rows must sum to \eqn{\le} 1.  The
#'   extra probability \(1 - \eqn{\sum_j p_{ij}}\) is interpreted as leaving
#'   the network from node *i*.
#' @param staClients  Warm-up completions discarded from statistics.
#' @param transitions Number of completed services counted for
#'   statistics.
#' @param historic Logical, collect whole trajectory?
#' @export
OpenNetSim <- R6::R6Class(
  "OpenNetSim",
  public = list(
    arrivalDistribution = NULL,
    serviceDistribution = NULL,
    servers             = NULL,
    routing             = NULL,
    staClients          = NULL,
    transitions         = NULL,
    historic            = NULL,
    out                 = NULL,

    initialize = function(arrivalDistribution,
                          serviceDistribution,
                          s,
                          p,
                          staClients  = 100L,
                          transitions = 1000L,
                          historic    = FALSE) {

      k <- length(s)
      stopifnot(is.list(arrivalDistribution),
                is.list(serviceDistribution),
                length(arrivalDistribution) == k,
                length(serviceDistribution) == k,
                is.numeric(s),
                is.matrix(p), nrow(p) == k, ncol(p) == k,
                any(vapply(arrivalDistribution,
                           function(x) !inherits(x, "no_distr"), logical(1))),
                staClients >= 0, transitions > 0)

      self$arrivalDistribution <- arrivalDistribution
      self$serviceDistribution <- serviceDistribution
      self$servers             <- as.integer(s)
      self$routing             <- p
      self$staClients          <- as.integer(staClients)
      self$transitions         <- as.integer(transitions)
      self$historic            <- isTRUE(historic)

      private$simulate_core()
    },

    print = function(...) {
      cat("Model: OpenNetSim  (open Jackson network)\n")
      with(self$out, {
        tab <- data.frame(L   = l,
                          Lq  = lq,
                          W   = w,
                          Wq  = wq,
                          row.names = paste0("Node", seq_along(l)))
        print(tab)
        cat("Total Lq =", lqt, "\n")
      })
      invisible(self)
    }
  ),

  private = list(
    simulate_core = function() {
      res <- OpenNetwork_secuential(self$arrivalDistribution,
                                    self$serviceDistribution,
                                    self$servers,
                                    self$routing,
                                    self$staClients,
                                    self$transitions,
                                    self$historic)
      self$out <- res$out
    }
  )
)

########################################################################
## Functional wrapper with replication / parallel                     ##
########################################################################

#' Simulate an open Jackson-type network
#'
#' Wrapper around `OpenNetSim$new()` plus optional parallel replication.
#'
#' @inheritParams OpenNetSim
#' @param nsim  Integer, number of independent replications.
#' @param nproc Integer, CPU cores (1 = sequential).
#' @return A `OpenNetSim` object, or the aggregated result of
#'   `combineSimulations()` when `nsim > 1`.
#' @export
OpenNet <- function(arrivalDistribution,
                    serviceDistribution,
                    s,
                    p,
                    staClients  = 100L,
                    transitions = 1000L,
                    historic    = FALSE,
                    nsim        = 10L,
                    nproc       = 1L) {

  if (nsim  < 1L) stop("'nsim' must be >= 1")
  if (nproc < 1L) stop("'nproc' must be >= 1")

  make_one <- function(arrivalDistribution, serviceDistribution,
                       s, p, sta, trn, hist) {
    OpenNetSim$new(arrivalDistribution,
                   serviceDistribution,
                   s,
                   p,
                   staClients  = sta,
                   transitions = trn,
                   historic    = hist)
  }

  args <- list(arrivalDistribution,
               serviceDistribution,
               s,
               p,
               staClients,
               transitions,
               historic)

  res <- ParallelizeSimulations(make_one, args, nsim, nproc)
  if (inherits(res, "OpenNetSim")) res else combineSimulations(res)
}

# ====================================================================
# Global variables for R CMD check ------------------------------------
# ====================================================================
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("L", "Lq", "W", "Wq", "Clients",
                           "Intensity", "tClient"))
}

# ====================================================================
