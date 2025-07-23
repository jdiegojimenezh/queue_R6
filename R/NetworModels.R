#' @include SimpleModels.R
NULL

########################################################################
# NETWORK MODELS (R6) -- Open and Closed Jackson networks              #
# All code is ASCII‑only for portability.                              #
########################################################################

## ------------------------------------------------------------------- ##
## Routing‑matrix helpers                                              ##
## ------------------------------------------------------------------- ##
#' Generate a cyclic routing matrix.
#'
#' @param n Integer. Number of nodes.
#' @keywords internal
Ciclic <- function(n) {
  diag(n)[c(2:n, 1), ]
}

#' Generate a serial (tandem) routing matrix.
#'
#' @inheritParams Ciclic
#' @keywords internal
Tamdem <- function(n) {
  diag(n + 1)[2:(n + 1), 1:n]
}

## ------------------------------------------------------------------- ##
## OpenJackson class (open network)                                    ##
## ------------------------------------------------------------------- ##
#' OpenJackson‑class: Open Jackson network (R6).
#'
#' @description R6 class that represents an *open* Jackson network with
#'   independent external Poisson arrivals. Each node behaves as an
#'   `M/M/s` queue (implemented internally with `MMS`).
#'
#' @section Fields (read‑only):
#' * `lambda_vec` – external arrival rates (numeric vector).
#' * `mu_vec`     – service rates at nodes.
#' * `servers_vec`– servers per node.
#' * `routing`    – routing matrix `P`.
#'
#' @section Key methods:
#' * `$Pn(n)`  – joint steady‑state probability `Pr{N = n}`.
#' * `$node(i)`– returns the `MMS` model of node *i*.
#' * `$print()`– console summary (overrides default).
#'
#' @param lambda Numeric vector of external arrival rates.
#' @param mu     Numeric vector of service rates.
#' @param s      Integer vector with the number of servers per node.
#' @param p      Routing matrix (square, rows sum <= 1).
#'
#' @export
OpenJackson <- R6::R6Class(
  "OpenJackson",
  inherit = MarkovianModel,
  public = list(
    # public read‑only accessors (active bindings defined below)
    nodes = NULL,

    initialize = function(lambda, mu, s, p) {
      stopifnot(is.numeric(lambda), is.numeric(mu), is.numeric(s),
                is.matrix(p), length(lambda) == length(mu),
                length(lambda) == length(s),
                nrow(p) == length(lambda), ncol(p) == length(lambda))

      private$.lambda <- as.numeric(lambda)
      private$.mu     <- as.numeric(mu)
      private$.servers <- as.integer(s)
      private$.P       <- p

      # verify routing matrix validity
      if (any(p < 0 | p > 1) || any(rowSums(p) > 1)) {
        stop("Argument 'p' must have values in [0,1] and row sums <= 1",
             call. = FALSE)
      }

      # compute effective arrival rates barlambda solving (I - P^T) * x = lambda
      A <- diag(length(lambda)) - t(p)
      if (abs(det(A)) <= .Machine$double.eps^0.5) {
        stop("Not stationary system (I - P^T singular)", call. = FALSE)
      }
      private$.barlambda <- solve(A, lambda)
      if (any(private$.barlambda >= mu * s)) {
        stop("Not stationary system (utilisation >= 1)", call. = FALSE)
      }

      # build per‑node MMS models
      self$nodes <- mapply(M_M_S,
                           lambda = private$.barlambda,
                           mu     = private$.mu,
                           s      = private$.servers,
                           SIMPLIFY = FALSE)

      private$compute_system_metrics()
    },

    # joint probability
    Pn = function(n) {
      stopifnot(length(n) == length(private$.lambda), all(n >= 0),
                all(n == floor(n)))
      probs <- mapply(Pn, self$nodes, as.integer(n))
      prod(probs)
    },

    # convenience accessor
    node = function(i) {
      if (i < 1 || i > length(self$nodes))
        stop("node index out of range", call. = FALSE)
      self$nodes[[i]]
    },

    print = function(...) {
      cat("Model: OpenJackson\n")
      L  <- private$.L_vec
      Lq <- private$.Lq_vec
      W  <- private$.W_vec
      Wq <- private$.Wq_vec
      tab <- rbind(Node = seq_along(L),
                   L    = L,
                   Lq   = Lq,
                   W    = W,
                   Wq   = Wq)
      print(t(tab))
      cat("Total L  =", private$.L_total,
          " | Total Lq =", private$.Lq_total, "\n")
    }
  ),
  active = list(
    lambda_vec  = function() private$.lambda,
    barlambda_vec = function() private$.barlambda,
    mu_vec      = function() private$.mu,
    servers_vec = function() private$.servers,
    routing     = function() private$.P
  ),
  private = list(
    .lambda    = NULL,  # external rates
    .mu        = NULL,
    .servers   = NULL,
    .P         = NULL,
    .barlambda = NULL,
    .L_vec     = NULL,
    .Lq_vec    = NULL,
    .W_vec     = NULL,
    .Wq_vec    = NULL,
    .L_total   = NULL,
    .Lq_total  = NULL,

    compute_system_metrics = function() {
      lq <- sapply(self$nodes, function(x) x$out$lq)
      L  <- lq + (private$.barlambda / private$.mu)
      W  <- lq / private$.barlambda + 1 / private$.mu
      Wq <- lq / private$.barlambda

      private$.L_vec   <- L
      private$.Lq_vec  <- lq
      private$.W_vec   <- W
      private$.Wq_vec  <- Wq
      private$.L_total <- sum(L)
      private$.Lq_total<- sum(lq)
      self$out <- list(L  = L, Lq = lq, W = W, Wq = Wq,
                       L_total = private$.L_total,
                       Lq_total = private$.Lq_total)
    }
  )
)

#' Functional constructor: Open Jackson network.
#'
#' @inheritParams OpenJackson
#' @return An `OpenJackson` object.
#' @export
OPEN_JACKSON <- function(lambda, mu, s, p) {
  OpenJackson$new(lambda = lambda, mu = mu, s = s, p = p)
}

## ------------------------------------------------------------------- ##
## ClosedJackson class (closed network)                                ##
## ------------------------------------------------------------------- ##
#' ClosedJackson‑class: Closed Jackson network (R6).
#'
#' @description R6 class for a *closed* Jackson network with a fixed
#'   population of `n` circulating customers. Each node behaves as an
#'   `M/M/s` queue.
#'
#' @param mu Numeric vector of service rates.
#' @param s  Integer vector with the number of servers per node.
#' @param p  Routing matrix. Rows must sum to 1.
#' @param n  Integer. Total number of customers in the network.
#'
#' @export
ClosedJackson <- R6::R6Class(
  "ClosedJackson",
  inherit = MarkovianModel,
  public = list(
    g_matrix = NULL,

    initialize = function(mu, s, p, n) {
      stopifnot(is.numeric(mu), is.numeric(s), is.matrix(p),
                length(mu) == length(s),
                nrow(p) == length(mu), ncol(p) == length(mu),
                is.numeric(n), n > 0)

      private$.mu      <- as.numeric(mu)
      private$.servers <- as.integer(s)
      private$.P       <- p
      private$.n       <- as.integer(n)
      private$.k       <- length(mu)
      self$servers     <- private$.servers

      if (any(p < 0 | p > 1) ||
          any(abs(rowSums(p) - 1) > 1e-12))
        stop("Routing matrix rows must sum to 1")

      private$compute_lambda_relative()
      private$compute_metrics()
    },

    # joint probability of state vector n_vec
    Pn = function(n_vec) {
      if (length(n_vec) != private$.k)
        stop("length(n_vec) != nodes")
      if (sum(n_vec) != private$.n)
        stop("Sum(n_vec) must equal total customers n")
      prod(mapply(f_close, MoreArgs = list(ps = self),
                  i = seq_len(private$.k), n = n_vec)) /
        private$.gkn
    },

    # marginal distribution in selected node
    Pi = function(n, node) {
      if (node < 1 || node > private$.k)
        stop("node out of range")
      n <- as.integer(n)
      out <- numeric(length(n))
      inside <- n >= 0 & n <= private$.n
      if (any(inside)) {
        out[inside] <-
          (f_close(self, node, n[inside]) *
             self$g_matrix[node - 1, private$.n - n[inside] + 1]) /
          private$.gkn
      }
      out
    },

    print = function(...) {
      cat("Model: ClosedJackson\n")
      M <- cbind(L  = private$.L_vec,
                 Lq = private$.Lq_vec,
                 W  = private$.W_vec,
                 Wq = private$.Wq_vec)
      rownames(M) <- seq_len(private$.k)
      print(M)
      cat("n =", private$.n,
          "| Total L =", sum(private$.L_vec),
          "| Total Lq =", sum(private$.Lq_vec), "\n")
    }
  ),
  active = list(
    mu_vec      = function() private$.mu,
    servers_vec = function() private$.servers,
    routing     = function() private$.P,
    lambda_vec  = function() private$.lambda_rel,
    barlambda_vec = function() private$.barlambda
  ),
  private = list(
    .mu = NULL, .servers = NULL, .P = NULL,
    .n = NULL, .k = NULL, .lambda_rel = NULL,
    .L_vec = NULL, .Lq_vec = NULL, .W_vec = NULL, .Wq_vec = NULL,
    .gkn = NULL, .barlambda = NULL,

    compute_lambda_relative = function() {
      k <- private$.k
      P <- private$.P
      if (k == 1) {
        private$.lambda_rel <- 1
        return()
      }

      M <- diag(k) - t(P)
      D <- M[-1, -1, drop = FALSE]
      c <- M[-1,  1, drop = FALSE]
      x <- -solve(D, c)

      private$.lambda_rel <- c(1, as.numeric(x))
    },

    compute_metrics = function() {
      k   <- private$.k
      n   <- private$.n
      mu       <- private$.mu
      s        <- private$.servers
      P        <- private$.P
      lambda_r <- private$.lambda_rel


      L_vec  <- numeric(k)
      Lq_vec <- numeric(k)
      W_vec  <- numeric(k)
      Wq_vec <- numeric(k)
      barlam <- numeric(k)

      shift <- c(k, 1:(k-1))
      gkn   <- NA
      ccte  <- NA

      for (idx in k:1) {
        ps <- list(mu     = mu,
                   servers = s,
                   out     = list(rho = lambda_r/mu),
                   k       = length(mu),
                   n       = n)
        class(ps) <- "partial"

        g <- calculateG(ps)

        if (is.na(gkn))
          gkn <- g[length(mu), n + 1]

        i <- 1:n
        Pn_last <- (f_close(ps, length(mu), i) * g[length(mu)-1, n - i + 1]) / gkn
        L <- sum(i * Pn_last)

        if (idx == k) {
          mult <- pmin(i, s[length(mu)])
          barlambda_k <- mu[length(mu)] * sum(mult * Pn_last)
          ccte <- barlambda_k / lambda_r[length(mu)]
        }

        barlambda <- ccte * lambda_r[length(mu)]
        W  <- L / barlambda
        Wq <- W - 1/mu[length(mu)]
        Lq <- barlambda * Wq

        L_vec [idx] <- L
        Lq_vec[idx] <- Lq
        W_vec [idx] <- W
        Wq_vec[idx] <- Wq
        barlam [idx] <- barlambda

        if (length(mu) > 1) {
          mu       <- mu[shift]
          s        <- s[shift]
          lambda_r <- lambda_r[shift]
          P        <- P[shift, shift]
        }
      }

      private$.gkn        <- gkn
      private$.L_vec      <- L_vec
      private$.Lq_vec     <- Lq_vec
      private$.W_vec      <- W_vec
      private$.Wq_vec     <- Wq_vec
      private$.barlambda  <- barlam

      ps0 <- list(mu = private$.mu, servers = private$.servers,
                  out = list(rho = private$.barlambda /(private$.mu * private$.servers)),
                  k = k, n = n)
      class(ps0) <- "partial"
      self$g_matrix <- calculateG(ps0)

      self$out <- list(
        rho = private$.barlambda / (private$.mu * private$.servers),
        L   = L_vec,
        Lq  = Lq_vec,
        W   = W_vec,
        Wq  = Wq_vec
      )
    }
  )
)

#' Functional constructor: Closed Jackson network.
#'
#' @inheritParams ClosedJackson
#' @return A `ClosedJackson` object.
#' @export
CLOSED_JACKSON <- function(mu, s, p, n) {
  ClosedJackson$new(mu = mu, s = s, p = p, n = n)
}

########################################################################
## Helpers reused from original S3 code but ASCII‑only ---------------- ##
########################################################################

f_close <- function(ps, i, n) {
  rho_i <- ps$out$rho[i]
  s_i   <- ps$servers[i]

  f <- function(m) {
    if (m == 0) return(1)
    if (m <= s_i) {
      return(rho_i ^ m / factorial(m))
    }
    rho_i ^ m / (factorial(s_i) * s_i ^ (m - s_i))
  }
  vapply(n, f, numeric(1))
}

calculateG <- function(ps) {
  k <- ps$k; n <- ps$n
  G <- matrix(0, nrow = k, ncol = n + 1)
  G[1, ] <- c(1, f_close(ps, 1, 1:n))           # primera fila

  if (k > 1) {
    for (j in 2:k) {                            # nodos
      G[j, 1] <- 1
      for (m in 2:(n + 1)) {                    # clientes 0..n
        G[j, m] <- sum(G[j - 1, m:1] *
                         f_close(ps, j, 0:(m - 1)))
      }
    }
  }
  G
}


pnlast <- function(ps, i) {
  if (ps$k == 1) return(1)
  if (min(i) < 0 || max(i) > ps$n) stop("i out of range in pnlast")
  (f_close(ps, ps$k, i) * ps$g[ps$k - 1, ps$n - i + 1]) / ps$g[ps$k, ps$n + 1]
}
