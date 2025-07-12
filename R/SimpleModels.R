#' @include MarkovianModel.R
NULL

# --------------------------------------------------------------------- #
# Utils ----------------------------------------------------------------#
# --------------------------------------------------------------------- #
.rate <- function(d) distr::rate(d)   # acceso rapido

# ===================================================================== #
#  MM1  --  M/M/1                                                       #
# ===================================================================== #
#' @title MM1-class: Class MM1
#' @description R6 class for the single-server exponential queue (M/M/1).
#' @inheritSection MarkovianModel Public methods
#' @param lambda Arrival rate (> 0).
#' @param mu     Service rate (> 0).
#' @examples
#' mm1 <- MM1$new(5, 8)
#' mm1$out$lq
#' Pn(mm1, 0:2)
#' FW(mm1, 1)
#' @usage NULL
#' @format An \\code{R6Class} generator.
#' @export
MM1 <- R6::R6Class(
  "MM1",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 3, mu = 6) {
      stopifnot(lambda > 0, mu > 0)
      private$lambda <- lambda
      private$mu     <- mu
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- 1L
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n >= 0L), all(n == floor(n)))
      rho <- self$out$rho
      (1 - rho) * rho^n
    },

    FW = function(x) {
      stopifnot(all(x >= 0))
      1 - exp((private$lambda - private$mu) * x)
    },

    FWq = function(x) {
      stopifnot(all(x >= 0))
      1 - (private$lambda / private$mu) *
          exp((private$lambda - private$mu) * x)
    }
  ),
  private = list(
    lambda = NA_real_,
    mu     = NA_real_,
    compute = function() {
      rho <- private$lambda / private$mu
      if (rho >= 1) stop("rho >= 1 (unstable queue).", call. = FALSE)
      l   <- private$lambda / (private$mu - private$lambda)
      wq  <- private$lambda / (private$mu * (private$mu - private$lambda))
      lq  <- private$lambda * wq
      w   <- l / private$lambda
      eff <- private$mu * w
      self$out <- list(rho = rho, barrho = rho,
                       l = l, lq = lq, wq = wq, w = w, eff = eff)
    }
  )
)

#' Functional constructor for MM1
#' @inheritParams MM1
#' @export
M_M_1 <- function(lambda = 3, mu = 6) MM1$new(lambda, mu)

# ===================================================================== #
#  MMS  --  M/M/s                                                       #
# ===================================================================== #
#' @title MMS-class: Class MMS
#' @description Multi-server exponential queue (M/M/s).
#' @inheritSection MarkovianModel Public methods
#' @param lambda Arrival rate (> 0).
#' @param mu     Service rate (> 0).
#' @param s      Number of servers (integer >= 1).
#' @export
MMS <- R6::R6Class(
  "MMS",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 3, mu = 6, s = 2L) {
      stopifnot(lambda > 0, mu > 0, s >= 1L)
      private$lambda <- lambda
      private$mu     <- mu
      private$s      <- as.integer(s)
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- private$s
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n >= 0L), all(n == floor(n)))
      n   <- as.integer(n)
      s   <- private$s
      rho <- self$out$rho
      cn  <- self$out$cn
      p0  <- self$out$p0
      pn  <- numeric(length(n))
      below <- n < s
      pn[below] <- p0 * cn[n[below] + 1L]
      if (any(!below)) {
        idx <- n[!below] - s
        pn[!below] <- p0 * cn[s + 1L] * rho^idx
      }
      pn
    },

    FWq = function(x) {
      stopifnot(all(x >= 0))
      lambda <- private$lambda; mu <- private$mu; s <- private$s
      cn_s1  <- self$out$cn[s + 1L]; p0 <- self$out$p0
      1 - cn_s1 * p0 * exp(-(s * mu - lambda) * x)
    },

    FW = function(x) {
      stopifnot(all(x >= 0))
      lambda <- private$lambda; mu <- private$mu; s <- private$s
      rho    <- lambda / (s * mu)
      cn_s1  <- self$out$cn[s + 1L]; p0 <- self$out$p0
      if (abs(rho - (s - 1)) < .Machine$double.eps^0.5) {
        1 - (1 + cn_s1 * p0 * x * mu) * exp(-mu * x)
      } else {
        a <- (lambda - s * mu + mu * FWq(self, 0)) /
             (s * mu - lambda - mu)
        b <- (cn_s1 * mu * p0) / (s * mu - lambda - mu)
        1 + a * exp(-mu * x) + b * exp(-(s * mu - lambda) * x)
      }
    }
  ),
  private = list(
    lambda = NA_real_, mu = NA_real_, s = NA_integer_,
    compute = function() {
      lambda <- private$lambda; mu <- private$mu; s <- private$s
      rho <- lambda / (s * mu)
      if (rho >= 1) stop("rho >= 1 (unstable queue).", call. = FALSE)
      cn <- if (s == 1)
               c(1, lambda / (s * mu - lambda))
            else
               c(1, lambda / ((1:(s - 1)) * mu),
                 lambda / (s * mu - lambda))
      cn <- cumprod(cn)
      p0 <- 1 / sum(cn)
      lq <- (cn[s + 1L] * lambda * p0) / (s * mu - lambda)
      wq <- lq / lambda
      w  <- wq + 1 / mu
      l  <- lambda * w
      eff <- mu * w
      self$out <- list(rho = rho, barrho = rho,
                       cn = cn, p0 = p0,
                       l = l, lq = lq, wq = wq, w = w, eff = eff)
    }
  )
)

#' Functional constructor for MMS
#' @inheritParams MMS
#' @export
M_M_S <- function(lambda = 3, mu = 6, s = 2L) MMS$new(lambda, mu, s)

# ===================================================================== #
#  MM1K --  M/M/1/K                                                     #
# ===================================================================== #
#' @title MM1K-class: Class MM1K
#' @description Single-server queue with finite capacity K (M/M/1/K).
#' @inheritSection MarkovianModel Public methods
#' @param lambda Arrival rate (> 0).
#' @param mu     Service rate (> 0).
#' @param k      Buffer size K (integer >= 0).
#' @export
MM1K <- R6::R6Class(
  "MM1K",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 3, mu = 6, k = 2L) {
      stopifnot(lambda > 0, mu > 0, k >= 0L)
      private$lambda <- lambda; private$mu <- mu; private$k <- as.integer(k)
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- 1L
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n >= 0L), all(n == floor(n)))
      n   <- as.integer(n)
      rho <- self$out$rho; k <- private$k
      ifelse(
        n > (k + 1L), 0,
        if (abs(rho - 1) < .Machine$double.eps^0.5)
          1 / (k + 2)
        else
          ((rho - 1) / (rho^(k + 2) - 1)) * rho^n
      )
    },

    Qn = function(n) {
      stopifnot(all(n >= 0L), all(n == floor(n)))
      k <- private$k
      ifelse(n > k, 0,
             Pn(self, n) / (1 - Pn(self, k + 1L)))
    },

    maxCustomers = function() private$k + 1L,

    FWq = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu; k <- private$k
      A <- S <- rep(1, length(x))
      B <- rep(Qn(self, 0), length(x))
      if (k >= 1) {
        for (n in 1:k) {
          A <- A * ((mu * x) / n)
          S <- S + A
          B <- B + Qn(self, n) * S
        }
      }
      1 - B * exp(-mu * x)
    },

    FW = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu
      vapply(x, function(t) {
        fwaux <- function(u) FWq(self, t - u) *
                              mu * exp(-mu * u)
        stats::integrate(fwaux, lower = 0, upper = t)$value
      }, numeric(1))
    }
  ),
  private = list(
    lambda = NA_real_, mu = NA_real_, k = NA_integer_,
    compute = function() {
      lambda <- private$lambda; mu <- private$mu; k <- private$k
      rho <- lambda / mu
      if (abs(rho - 1) < .Machine$double.eps^0.5) {
        l <- (k + 1) / 2
        barlambda <- lambda * (k + 1) / (k + 2)
      } else {
        l <- (rho / (1 - rho)) -
             ((k + 2) * rho^(k + 2) /
              (1 - rho^(k + 2)))
        barlambda <- lambda *
                     (rho^(k + 1) - 1) /
                     (rho^(k + 2) - 1)
      }
      w  <- l / barlambda
      wq <- w - 1 / mu
      lq <- barlambda * wq
      eff <- mu * w
      barrho <- barlambda / mu
      self$out <- list(rho = rho, barrho = barrho, barlambda = barlambda,
                       l = l, lq = lq, wq = wq, w = w, eff = eff)
    }
  )
)

#' Functional constructor for MM1K
#' @inheritParams MM1K
#' @export
M_M_1_K <- function(lambda = 3, mu = 6, k = 2L) MM1K$new(lambda, mu, k)

# ===================================================================== #
#  MMSK --  M/M/s/K                                                     #
# ===================================================================== #
#' @title MMSK-class: Class MMSK
#' @description Multi-server queue with capacity K (M/M/s/K).
#' @inheritSection MarkovianModel Public methods
#' @param lambda Arrival rate (> 0).
#' @param mu     Service rate (> 0).
#' @param s      Servers (integer >= 1).
#' @param k      Capacity K (integer >= 0).
#' @export
MMSK <- R6::R6Class(
  "MMSK",
  inherit = MarkovianModel,

  public = list(
    initialize = function(lambda = 3, mu = 6,
                          s = 2L, k = 3L) {
      stopifnot(lambda > 0, mu > 0, s >= 1L, k >= 0L)
      private$lambda <- lambda
      private$mu     <- mu
      private$s      <- as.integer(s)
      private$k      <- as.integer(k)
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- private$s
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n >= 0L), all(n == floor(n)))
      lambda <- private$lambda; mu <- private$mu
      s      <- private$s;      k  <- private$k
      p0     <- self$out$p0
      n      <- as.integer(n)

      pn <- numeric(length(n))
      below <- n < s
      if (any(below)) {
        pn[below] <- p0 * (lambda / mu) ^ n[below] /
          factorial(n[below])
      }
      if (any(!below & n <= s + k)) {
        idx           <- (!below) & (n <= s + k)
        nn            <- n[idx]
        pn[idx] <- p0 * (lambda / mu) ^ nn /
          (factorial(s) * s ^ (nn - s))
      }
      pn
    },

  Qn = function(n) {                       # ←  MÉTODO QUE FALTABA
    stopifnot(all(n >= 0L), all(n == floor(n)))
    limit <- private$s + private$k - 1L
    ifelse(n > limit, 0,
           Pn(self, n) / (1 - self$out$p_block))
  },

  maxCustomers = function() private$s + private$k,

  FWq = function(x) {
    stopifnot(all(x >= 0))
    mu <- private$mu; s <- private$s; k <- private$k
    # distribución Erlang truncada (mismo esquema que antes)
    A <- S <- rep(1, length(x))
    B <- rep(self$Qn(s), length(x))          # P(cola >= 1)
    if ((s + 1) <= (s + k - 1)) {
      for (n in (s + 1):(s + k - 1)) {
        A <- A * ((s * mu * x) / (n - s))
        S <- S + A
        B <- B + self$Qn(n) * S
      }
    }
    1 - B * exp(-s * mu * x)
  },

  FW = function(x) {
    stopifnot(all(x >= 0))
    mu <- private$mu
    vapply(x, function(t) {
      fwaux <- function(u) self$FWq(t - u) * mu * exp(-mu * u)
      stats::integrate(fwaux, lower = 0, upper = t)$value
    }, numeric(1))
  }
  ),


  private = list(
    lambda = NA_real_, mu = NA_real_,
    s = NA_integer_,   k = NA_integer_,

    ## --- compute() totalmente reescrito ----------------------------- ##
    compute = function() {
      lambda <- private$lambda; mu <- private$mu
      s      <- private$s;      k  <- private$k
      rho    <- lambda / (s * mu)              # utilización por servidor

      ## 1. Constante de normalización p0 ----------------------------- ##
      terms <- numeric(s + k + 1)              # n = 0 .. s+k
      for (n in 0:(s + k)) {
        if (n < s) {
          terms[n + 1] <- (lambda / mu) ^ n / factorial(n)
        } else {
          terms[n + 1] <- (lambda / mu) ^ n /
            (factorial(s) * s ^ (n - s))
        }
      }
      p0  <- 1 / sum(terms)
      Pn0 <- p0 * terms                         # vector de probabilidades

      ## 2. Probabilidad de bloqueo (s + k) --------------------------- ##
      p_block <- Pn0[s + k + 1]                # estado lleno

      ## 3. Llegada efectiva y métricas globales ---------------------- ##
      bar_lambda <- lambda * (1 - p_block)     # λ efectivo
      L  <- sum(0:(s + k) * Pn0)               # clientes en sistema
      Ls <- bar_lambda / mu                    # clientes siendo servidos
      Lq <- L - Ls                             # clientes en cola
      W  <- L  / bar_lambda                    # Little
      Wq <- Lq / bar_lambda
      eff <- mu * W

      self$out <- list(rho      = rho,
                       p0       = p0,
                       p_block  = p_block,
                       l        = L,
                       lq       = Lq,
                       w        = W,
                       wq       = Wq,
                       eff      = eff)
    }
  )
)

## wrapper sin cambios
M_M_S_K <- function(lambda = 3, mu = 6, s = 2L, k = 3L) {
  MMSK$new(lambda, mu, s, k)
}

