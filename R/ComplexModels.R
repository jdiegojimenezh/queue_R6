#' @include MarkovianModel.R
NULL

# --------------------------------------------------------------------- #
# Helper ----------------------------------------------------------------
.rate <- function(d) distr::rate(d)   # shorthand

# ===================================================================== #
#  MM1InfH  --  M/M/1/Inf/H                                            #
# ===================================================================== #
#' @title MM1InfH-class: M/M/1/Inf/H finite-population queue
#' @description R6 class for a single-server Markovian queue with a finite
#'   population of size \\eqn{H}.  Kendall notation: \\emph{M/M/1/Inf/H}.
#' @inheritSection MarkovianModel Public methods
#' @param lambda Mean arrival rate (> 0).
#' @param mu     Mean service rate  (> 0).
#' @param h      Population size (integer >= 1).
#' @examples
#' mm <- MM1InfH$new(0.5, 12, 5)
#' mm$out$l; Pn(mm, 0:3); FWq(mm, 1)
#' @format An \\code{R6Class} generator.
#' @export
MM1InfH <- R6::R6Class(
  "MM1InfH",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 0.5, mu = 12, h = 5L) {
      stopifnot(lambda > 0, mu > 0, h >= 1L)
      private$lambda <- lambda; private$mu <- mu; private$h <- as.integer(h)
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- 1L
      private$compute()
    },

    # ------------------ Steady-state probabilities ------------------ #
    Pn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      n  <- as.integer(n)
      h  <- private$h
      cn <- self$out$cn; p0 <- self$out$p0
      ifelse(n > h, 0, p0 * cn[n + 1L])
    },

    Qn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      h <- private$h
      ifelse(n > (h - 1L), 0,
             ((h - n) * Pn(self, n)) / (h - self$out$l))
    },

    # --------------- Waiting-time distributions --------------------- #
    FWq = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu; h <- private$h
      A <- S <- rep(1, length(x))
      B <- rep(Qn(self, 0), length(x))
      if (h > 1) {
        for (n in 1:(h - 1L)) {
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
        fwaux <- function(u) FWq(self, t - u) * mu * exp(-mu * u)
        stats::integrate(fwaux, lower = 0, upper = t)$value
      }, numeric(1))
    },

    maxCustomers = function() private$h
  ),
  private = list(
    lambda = NA_real_, mu = NA_real_, h = NA_integer_,
    compute = function() {
      lambda <- private$lambda; mu <- private$mu; h <- private$h
      rho <- lambda / mu
      cn <- c(1, rho * (h - (1:h) + 1))
      cn <- cumprod(cn)
      p0 <- 1 / sum(cn)
      l  <- sum((1:h) * cn[-1] * p0)
      barlambda <- lambda * (h - l)
      w  <- l / barlambda
      wq <- w - 1 / mu
      lq <- barlambda * wq
      eff <- mu * w
      barrho <- rho * (h - l)
      self$out <- list(rho = rho, barrho = barrho, barlambda = barlambda,
                       cn = cn, p0 = p0,
                       l = l, lq = lq, wq = wq, w = w, eff = eff)
    }
  )
)

#' Functional constructor for MM1InfH
#' @inheritParams MM1InfH
#' @export
M_M_1_INF_H <- function(lambda = 0.5, mu = 12, h = 5L)
  MM1InfH$new(lambda, mu, h)

# ===================================================================== #
#  MMSInfH  --  M/M/s/Inf/H                                            #
# ===================================================================== #
#' @title MMSInfH-class: M/M/s/Inf/H finite-population queue
#' @inheritSection MarkovianModel Public methods
#' @param lambda Mean arrival rate (> 0).
#' @param mu     Mean service rate  (> 0).
#' @param s      Servers (integer >= 1).
#' @param h      Population size (>= s).
#' @export
MMSInfH <- R6::R6Class(
  "MMSInfH",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 0.5, mu = 12, s = 2L, h = 5L) {
      stopifnot(lambda > 0, mu > 0, s >= 1L, h >= s)
      private$lambda <- lambda; private$mu <- mu
      private$s <- as.integer(s); private$h <- as.integer(h)
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- private$s
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      n <- as.integer(n)
      h <- private$h; s <- private$s
      cn <- self$out$cn; p0 <- self$out$p0; rho <- self$out$rho
      ifelse(n > h, 0, {
        below <- n < s
        pn <- numeric(length(n))
        pn[below] <- p0 * cn[n[below] + 1L]
        if (any(!below))
          pn[!below] <- p0 * cn[s + 1L] * rho^(n[!below] - s)
        pn
      })
    },

    Qn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      h <- private$h
      ifelse(n > (h - 1L), 0,
             ((h - n) * Pn(self, n)) / (h - self$out$l))
    },

    FWq = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu; s <- private$s; h <- private$h
      A <- S <- rep(1, length(x))
      B <- rep(Qn(self, s), length(x))
      if (s <= (h - 1L)) {
        for (n in (s + 1):(h - 1L)) {
          A <- A * ((s * mu * x) / (n - s))
          S <- S + A
          B <- B + Qn(self, n) * S
        }
      }
      1 - B * exp(-s * mu * x)
    },

    FW = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu
      vapply(x, function(t) {
        fwaux <- function(u) FWq(self, t - u) * mu * exp(-mu * u)
        stats::integrate(fwaux, lower = 0, upper = t)$value
      }, numeric(1))
    },

    maxCustomers = function() private$h
  ),
  private = list(
    lambda = NA_real_, mu = NA_real_, s = NA_integer_, h = NA_integer_,
    compute = function() {
      lambda <- private$lambda; mu <- private$mu
      s <- private$s; h <- private$h
      rho <- lambda / (s * mu)

      cn <- c(
        1,
        (lambda / ((1:s) * mu)) * (h - (1:s) + 1),
        rho * (h - ((s + 1):h) + 1)
      )
      cn <- cumprod(cn)
      p0 <- 1 / sum(cn)
      l  <- sum((1:h) * cn[-1] * p0)
      barlambda <- lambda * (h - l)
      w  <- l / barlambda
      wq <- w - 1 / mu
      lq <- barlambda * wq
      eff <- mu * w
      barrho <- rho * (h - l)
      self$out <- list(rho = rho, barrho = barrho, barlambda = barlambda,
                       cn = cn, p0 = p0,
                       l = l, lq = lq, wq = wq, w = w, eff = eff)
    }
  )
)

#' Functional constructor for MMSInfH
#' @inheritParams MMSInfH
#' @export
M_M_S_INF_H <- function(lambda = 0.5, mu = 12, s = 2L, h = 5L)
  MMSInfH$new(lambda, mu, s, h)

# ===================================================================== #
#  MMSInfHY --  M/M/s/Inf/H/Y                                           #
# ===================================================================== #
#' @title MMSInfHY-class: M/M/s/Inf/H with Y replacements
#' @inheritSection MarkovianModel Public methods
#' @param lambda Mean arrival rate (> 0).
#' @param mu     Mean service rate  (> 0).
#' @param s      Servers (>= 1).
#' @param h      Population size (> 0).
#' @param y      Number of replacements (>= 1).
#' @export
MMSInfHY <- R6::R6Class(
  "MMSInfHY",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 3, mu = 6,
                          s = 3L, h = 5L, y = 3L) {
      stopifnot(lambda > 0, mu > 0, s >= 1L,
                h > 0, y > 0, s <= y + h)
      private$lambda <- lambda; private$mu <- mu
      private$s <- as.integer(s); private$h <- as.integer(h)
      private$y <- as.integer(y)
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- private$s
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      n <- as.integer(n)
      h <- private$h; y <- private$y
      cn <- self$out$cn; p0 <- self$out$p0
      ifelse(n > (y + h), 0, p0 * cn[n + 1L])
    },

    Qn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      h <- private$h; y <- private$y
      emes <- y:(y + h)
      sumemes <- sum((emes - y) * Pn(self, emes))
      ifelse(n > (y + h - 1L), 0, {
        if (n <= (y - 1L))
          (h * Pn(self, n)) / (h - sumemes)
        else
          ((h + y - n) * Pn(self, n)) / (h - sumemes)
      })
    },

    FWq = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu; s <- private$s
      Nmax <- private$y + private$h - 1L
      A <- S <- rep(1, length(x))
      B <- rep(Qn(self, s), length(x))
      if (s <= Nmax) {
        for (n in (s + 1):Nmax) {
          A <- A * ((s * mu * x) / (n - s))
          S <- S + A
          B <- B + Qn(self, n) * S
        }
      }
      1 - B * exp(-s * mu * x)
    },

    FW = function(x) {
      stopifnot(all(x >= 0))
      mu <- private$mu
      vapply(x, function(t) {
        fwaux <- function(u) FWq(self, t - u) * mu * exp(-mu * u)
        stats::integrate(fwaux, lower = 0, upper = t)$value
      }, numeric(1))
    },

    maxCustomers = function() private$y + private$h
  ),
  private = list(
    lambda = NA_real_, mu = NA_real_, s = NA_integer_,
    h = NA_integer_, y = NA_integer_,
    compute = function() {
      lambda <- private$lambda; mu <- private$mu
      s <- private$s; h <- private$h; y <- private$y
      rho <- lambda / (s * mu)

      if (s <= y) {
        cn <- c(
          1,
          (h * lambda) / ((1:s) * mu),
          rep(rho * h, y - s),
          rho * (h + y - ((y + 1):(y + h)) + 1)
        )
      } else {
        cn <- c(
          1,
          (h * lambda) / ((1:y) * mu),
          ((h + y - ((y + 1):s) + 1) * lambda) /
            (((y + 1):s) * mu),
          (h + y - ((s + 1):(y + h)) + 1) * rho
        )
      }
      cn <- cumprod(cn)
      p0 <- 1 / sum(cn)

      n_vec <- y:(y + h)
      barlambda <- lambda *
        (h - sum((n_vec - y) * p0 * cn[n_vec + 1L]))
      barrho <- barlambda / (s * mu)
      n_full <- 1:(y + h)
      l <- sum(n_full * p0 * cn[n_full + 1L])
      w  <- l / barlambda
      wq <- w - 1 / mu
      lq <- barlambda * wq
      eff <- mu * w

      self$out <- list(rho = rho, barrho = barrho, barlambda = barlambda,
                       cn = cn, p0 = p0,
                       l = l, lq = lq, wq = wq, w = w, eff = eff)
    }
  )
)

#' Functional constructor for MMSInfHY
#' @inheritParams MMSInfHY
#' @export
M_M_S_INF_H_Y <- function(lambda = 3, mu = 6,
                          s = 3L, h = 5L, y = 3L)
  MMSInfHY$new(lambda, mu, s, h, y)

# ===================================================================== #
#  MMInf  --  M/M/Inf                                                  #
# ===================================================================== #
#' @title MMInf-class: M/M/Inf infinite-server queue
#' @inheritSection MarkovianModel Public methods
#' @param lambda Mean arrival rate (> 0).
#' @param mu     Mean service rate  (> 0).
#' @export
MMInf <- R6::R6Class(
  "MMInf",
  inherit = MarkovianModel,
  public = list(
    initialize = function(lambda = 3, mu = 6) {
      stopifnot(lambda > 0, mu > 0)
      private$lambda <- lambda; private$mu <- mu
      super$initialize(distr::Exp(lambda), distr::Exp(mu))
      self$servers <- Inf
      private$compute()
    },

    Pn = function(n) {
      stopifnot(all(n == floor(n)), all(n >= 0))
      n <- as.integer(n)
      lambda <- private$lambda; mu <- private$mu
      p0 <- self$out$p0
      cn <- if (max(n) > 0)
        cumprod(c(1, lambda / ((1:max(n)) * mu)))
      else 1
      p0 * cn[n + 1L]
    },

    Qn  = function(n) rep(0, length(n)),   # no queue
    FWq = function(x) rep(0, length(x)),
    FW  = function(x) stats::pexp(x, rate = private$mu),
    maxCustomers = function() Inf
  ),
  private = list(
    lambda = NA_real_, mu = NA_real_,
    compute = function() {
      lambda <- private$lambda; mu <- private$mu
      rho <- lambda / mu
      p0 <- exp(-rho)
      l  <- rho
      w  <- 1 / mu
      self$out <- list(rho = rho, barrho = 0, p0 = p0,
                       l = l, lq = 0, wq = 0, w = w, eff = 1)
    }
  )
)

#' Functional constructor for MMInf
#' @inheritParams MMInf
#' @export
M_M_INF <- function(lambda = 3, mu = 6) MMInf$new(lambda, mu)
