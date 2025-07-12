#' @title Ajuste y contraste de distribuciones continuas
#' @description Conjunto de funciones auxiliares para:
#' * Ajustar varias distribuciones mediante \pkg{fitdistrplus}
#' * Evaluar bondad de ajuste (chi-cuadrado y Kolmogorov–Smirnov)
#' * Resumir gráficamente densidad, CDF y Q-Q plot con \pkg{ggplot2}
#' @name DistributionAnalysis-internal
#' @keywords internal
#' @section Dependencias:
#' Estas funciones requieren los paquetes:
#' \itemize{
#'   \item \pkg{fitdistrplus} – ajuste MLE y gráficos base (`denscomp`, …)
#'   \item \pkg{ggplot2} y \pkg{grid} – visualizaciones opcionales
#' }
#'
#' @importFrom fitdistrplus fitdist denscomp cdfcomp qqcomp
#' @importFrom stats pchisq ks.test
#' @importFrom graphics hist curve layout
#' @import ggplot2
#' @import grid
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("x", "dist", "theor"))
}

#' Ajusta varias distribuciones a un vector numérico
#'
#' @param data Vector numérico con las observaciones
#' @param ldistr `character` con nombres abreviados de distribuciones
#'               (ej. `"exp"`, `"norm"`, `"weibull"`, …)
#' @return Lista de objetos \code{fitdist}.  Clase extra: `"FitList"`
#' @family DistributionAnalysis
#' @export
fitData <- function(data,
                    ldistr = c("exp", "norm", "weibull",
                               "unif", "lnorm", "gamma", "beta")) {

  if (is.null(data))
    stop("Argument 'data' must be a numeric vector.")

  options(warn = -1)               # silenciamos avisos de fitdistrplus
  res <- lapply(ldistr, function(dn) {
    tryCatch(
      fitdist(data, dn, method = "mle"),
      error = function(e) NA      # devolvemos NA si el ajuste falla
    )
  })
  names(res) <- ldistr
  res <- res[!is.na(res)]          # quitamos fallos
  class(res) <- c("FitList", class(res))
  res
}

#' Pruebas chi-cuadrado y KS para cada ajuste
#'
#' @param lfitdata Lista producida por \code{\link{fitData}}
#' @return \code{data.frame} con: distribución, estadísticos y *p-values*
#' @family DistributionAnalysis
#' @export
goodnessFit <- function(lfitdata) {

  out <- data.frame(distrnames = character(),
                    chisq      = numeric(),
                    chisq.pvalue = character(),
                    ks         = numeric(),
                    ks.pvalue  = character(),
                    stringsAsFactors = FALSE)

  for (fit in lfitdata) {

    ## --- Kolmogorov–Smirnov -----------------------------------------
    ks_call <- as.list(fit$estimate)
    ks_fun  <- paste0("p", fit$distname)
    ks_res  <- ks.test(fit$data, ks_fun, !!!ks_call)

    ## --- Chi-cuadrado personalizado ---------------------------------
    chisq_res <- do.call(
      chisq.test.cont,
      c(list(x          = fit$data,
             distribution = fit$distname),
        as.list(fit$estimate))
    )

    out <- rbind(out,
                 data.frame(distrnames   = fit$distname,
                            chisq        = chisq_res$statistic,
                            chisq.pvalue = format(round(chisq_res$p.value, 3),
                                                  nsmall = 3),
                            ks           = ks_res$statistic,
                            ks.pvalue    = format(round(ks_res$p.value, 3),
                                                  nsmall = 3)))
  }
  class(out) <- c("GoodnessFit", class(out))
  out
}

#' Resumen gráfico (densidad, CDF y Q-Q plot)
#'
#' @param lfitdata Salida de \code{\link{fitData}}
#' @param graphics `"graphics"` o `"ggplot2"`
#' @param show     `"all"`, `"dens"`, `"cdf"` o `"qq"`
#' @family DistributionAnalysis
#' @export
summaryFit <- function(lfitdata,
                       graphics = c("ggplot2", "graphics"),
                       show     = c("all", "dens", "cdf", "qq")) {

  graphics <- match.arg(graphics)
  show     <- match.arg(show)

  labs <- if (names(lfitdata)[1] == "estimate")
    lfitdata$distname else names(lfitdata)

  if (graphics == "graphics") {
    ## Base R ---------------------------------------------------------
    if (show == "all") layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    if (show %in% c("all","dens")) denscomp(lfitdata, legendtext = labs)
    if (show %in% c("all","cdf"))  cdfcomp (lfitdata, legendtext = labs)
    if (show %in% c("all","qq"))   qqcomp  (lfitdata, legendtext = labs)

  } else {
    ## ggplot2 --------------------------------------------------------
    switch(show,
           "all"  = {
             grid::grid.newpage()
             lay <- grid::grid.layout(2, 2)
             grid::pushViewport(grid::viewport(layout = lay))
             print(denscompggplot2(lfitdata),
                   vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2))
             print(cdfcompggplot2(lfitdata),
                   vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
             print(qqcompggplot2(lfitdata),
                   vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
           },
           "dens" = denscompggplot2(lfitdata),
           "cdf"  = cdfcompggplot2(lfitdata),
           "qq"   = qqcompggplot2(lfitdata)
    )
  }
}

denscompggplot2 <- function(lfitdata) {
  y <- if (names(lfitdata)[1] == "estimate") lfitdata$data else lfitdata[[1]]$data
  aux <- hist(y, plot = FALSE)
  dfh <- data.frame(x = aux$mids, y = aux$density)

  p <- ggplot(dfh, aes(x, y)) +
    geom_histogram(stat = "identity", binwidth = diff(aux$breaks)[1],
                   fill = "grey80", colour = NA)

  xseq <- seq(min(y), max(y), length.out = 512)
  for (fit in lfitdata) {
    dens <- do.call(paste0("d", fit$distname),
                    c(list(xseq), as.list(fit$estimate)))
    p <- p + geom_line(data = data.frame(x = xseq, y = dens,
                                         dist = fit$distname),
                       aes(colour = dist), linewidth = .3)
  }
  p + labs(x = NULL, y = "Density", colour = "Dist.")
}

cdfcompggplot2 <- function(lfitdata) {
  y <- if (names(lfitdata)[1] == "estimate") lfitdata$data else lfitdata[[1]]$data
  y <- sort(y); n <- length(y)
  df_ecdf <- data.frame(x = y, y = (1:n)/n)

  p <- ggplot(df_ecdf, aes(x, y)) + geom_step()

  xseq <- seq(min(y), max(y), length.out = 512)
  for (fit in lfitdata) {
    cdf <- do.call(paste0("p", fit$distname),
                   c(list(xseq), as.list(fit$estimate)))
    p <- p + geom_line(data = data.frame(x = xseq, y = cdf,
                                         dist = fit$distname),
                       aes(colour = dist), linewidth = .3)
  }
  p + labs(x = NULL, y = "F(x)", colour = "Dist.")
}

qqcompggplot2 <- function(lfitdata) {
  y <- if (names(lfitdata)[1] == "estimate") lfitdata$data else lfitdata[[1]]$data
  y <- sort(y); n <- length(y)
  p <- ggplot(data.frame(sample = y), aes(sample = sample)) +
    geom_point(aes(sample, sample), colour = "grey60", size = .1)

  probs <- (1:n)/(n + 1)
  for (fit in lfitdata) {
    qth <- do.call(paste0("q", fit$distname),
                   c(list(probs), as.list(fit$estimate)))
    p <- p + geom_line(data = data.frame(sample = y, theor = qth,
                                         dist = fit$distname),
                       aes(sample, theor, colour = dist), linewidth = .3)
  }
  p + labs(x = "Empirical quantiles", y = "Theoretical quantiles",
           colour = "Dist.")
}

#' @keywords internal
chisq.test.cont <- function(x, distribution = "norm",
                            nclasses = min(100L, max(3L, floor(length(x)/5))),
                            output = FALSE, nestpar = 0, ...) {

  qfun <- get(paste0("q", distribution), mode = "function")
  dfun <- get(paste0("d", distribution), mode = "function")

  breaks <- c(-Inf,
              qfun((1:(nclasses - 1))/nclasses, ...),
              Inf)

  O <- hist(x, breaks = breaks, plot = FALSE)$counts
  E <- length(x) / nclasses

  STAT <- sum((O - E)^2/E)
  df   <- nclasses - nestpar - 1
  pval <- stats::pchisq(STAT, df, lower.tail = FALSE)

  structure(list(statistic = c("X-squared" = STAT),
                 parameter = c(df = df),
                 p.value   = pval,
                 method    = "Pearson's Chi-squared test",
                 data.name = deparse(substitute(x))),
            class = "htest")
}
