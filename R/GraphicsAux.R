# ====================================================================
#  plot_helpers.R  –  generic plotting & small S3 utilities
# --------------------------------------------------------------------
#  These helpers work with any object that stores a `$out$historic`
#  matrix/array following the convention used by the new R6 models.
# ====================================================================

#' @importFrom stats aggregate
#' @importFrom utils globalVariables
NULL

#' Internal helper: single-queue history
#' @noRd
.plot_single_history <- function(hist, var = "L", minrange = 1L,
                                 maxrange = nrow(hist), depth = maxrange - minrange + 1L) {
  rng  <- seq.int(minrange, maxrange, length.out = depth)
  dat  <- data.frame(n = rng, value = hist[rng, var])
  ggplot2::ggplot(dat, ggplot2::aes(x = n, y = value)) +
    ggplot2::geom_line() +
    ggplot2::xlab("n") + ggplot2::ylab(var) +
    ggplot2::ggtitle(paste("Evolution of", var))
}

#' Internal helper: network history (array [node, var, time])
#' @noRd
.plot_network_history <- function(hist, var = "L", minrange = 1L,
                                  maxrange = dim(hist)[3L],
                                  depth = maxrange - minrange + 1L) {
  rng   <- seq.int(minrange, maxrange, length.out = depth)
  nodes <- dim(hist)[1L]
  dat   <- do.call(rbind, lapply(seq_len(nodes), function(i) {
    data.frame(n     = rng,
               value = hist[i, var, rng],
               node  = factor(i))
  }))
  ggplot2::ggplot(dat, ggplot2::aes(x = n, y = value, colour = node, group = node)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::xlab("n") + ggplot2::ylab(var) +
    ggplot2::ggtitle(paste("Evolution of", var)) +
    ggplot2::scale_colour_discrete(name = "Node")
}

#' Internal helper: list of simulations
#' @noRd
.plot_list_history <- function(obj_list, var = "L", minrange = 1L,
                               maxrange = NULL,
                               depth     = NULL,
                               showMean  = TRUE,
                               showValues = TRUE) {

  first_hist <- obj_list[[1L]]$out$historic
  if (is.null(maxrange)) {
    maxrange <- if (is.matrix(first_hist))
      nrow(first_hist) else dim(first_hist)[3L]
  }
  if (is.null(depth)) depth <- maxrange - minrange + 1L
  rng <- seq.int(minrange, maxrange, length.out = depth)

  # build long data frame --------------------------------------------
  dat <- do.call(rbind, Map(function(obj, sim_id) {
    h <- obj$out$historic
    if (is.matrix(h)) {
      data.frame(n     = rng,
                 value = h[rng, var],
                 sim   = factor(sim_id))
    } else {
      # network
      do.call(rbind, lapply(seq_len(dim(h)[1L]), function(node) {
        data.frame(n     = rng,
                   value = h[node, var, rng],
                   node  = factor(node),
                   sim   = factor(sim_id))
      }))
    }
  }, obj_list, seq_along(obj_list)))

  g <- ggplot2::ggplot()

  if (showValues)
    g <- g + ggplot2::geom_line(
      data = dat,
      ggplot2::aes(x = n, y = value, group = interaction(sim, node),
                   colour = node, alpha = 0.8),
      na.rm = TRUE
    ) +
    ggplot2::scale_alpha_continuous(name = NULL, breaks = NULL, labels = NULL)

  if (showMean) {
    if ("node" %in% names(dat)) {
      mean_dat <- aggregate(value ~ n + node, dat, mean, na.rm = TRUE)
      g <- g + ggplot2::geom_line(
        data = mean_dat,
        ggplot2::aes(x = n, y = value, colour = node),
        size = 1
      )
    } else {
      mean_dat <- aggregate(value ~ n, dat, mean, na.rm = TRUE)
      g <- g + ggplot2::geom_line(
        data = mean_dat,
        ggplot2::aes(x = n, y = value),
        colour = "red", size = 1
      )
    }
  }

  g + ggplot2::xlab("n") + ggplot2::ylab(var) +
    ggplot2::ggtitle(paste("Evolution of", var)) +
    ggplot2::scale_colour_discrete(name = "Node")
}

#' Quick historic plot for any simulated queue or network
#'
#' @param x   A single simulation object (`GG1Sim`, `OpenNetSim`, …) **or**
#'            a list of such objects.
#' @param var One of \code{"L"}, \code{"Lq"}, \code{"W"}, \code{"Wq"},
#'            \code{"Clients"}, \code{"Intensity"}.
#' @param ... Extra parameters forwarded to the internal helpers
#'            (e.g., \code{minrange}, \code{maxrange}, \code{depth},
#'            \code{showMean}, \code{showValues}).
#' @return A \pkg{ggplot2} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_histogram geom_vline
#'             xlab ylab ggtitle scale_colour_discrete theme element_blank
#'             scale_alpha_continuous
plot_history <- function(x, var = "L", ...) {
  stopifnot(var %in% c("L", "Lq", "W", "Wq", "Clients", "Intensity"))

  if (is.list(x) && !inherits(x, "SimulatedModel")) {
    .plot_list_history(x, var, ...)
  } else if (inherits(x, "SimulatedNetwork")) {
    .plot_network_history(x$out$historic, var, ...)
  } else if (inherits(x, "SimulatedModel")) {
    .plot_single_history(x$out$historic, var, ...)
  } else {
    stop("Unsupported object type for 'plot_history()'")
  }
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("n", "value", "sim", "node"))
}

