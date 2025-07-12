#' @import doParallel
#' @import foreach
#' @importFrom foreach %dopar%
#' @importFrom stats sd
NULL

#' Comprueba si un objeto pertenece a un paquete concreto
#'
#' @param obj  Objeto a inspeccionar
#' @param packagename  Nombre del paquete
#' @return `TRUE` si **obj** proviene de *packagename*, `FALSE` en otro caso
#' @keywords internal
belong <- function(obj, packagename) {
  if (length(obj) > 1) {
    packageobj <- sapply(lapply(obj, class), attr, "package")
  } else {
    packageobj <- attr(class(obj), "package")
  }
  all(packageobj == packagename)
}

#' Agrega un listado de simulaciones independientes
#'
#' Combina varias réplicas de un modelo de colas calculando promedio,
#' desviación típica y un resumen rápido para cada métrica de interés.
#'
#' @param listsims `list` con objetos de clase (S3 o R6) simulada
#' @return Un único objeto de la misma clase que `listsims[[1]]`,
#'   pero con los campos `$out` reemplazados por listas `mean/sd/summary`
#' @export
combineSimulations <- function(listsims) {

  if (!inherits(listsims, "list")) {
    # Si el usuario nos pasó directamente un solo objeto, lo devolvemos tal cual
    if (inherits(listsims, "SimulatedModel")) return(listsims)
    stop("'listsims' debe ser una lista de modelos simulados")
  }
  if (length(listsims) == 1L) return(listsims[[1]])

  # ------------------------------------------------------------------- #
  #  Pequeña función helper para extraer un atributo dentro de $out     #
  # ------------------------------------------------------------------- #
  get_out_attr <- function(qm, attr) qm$out[[attr]]

  metrics <- c("l", "lq", "w", "wq", "rho", "eff")
  mat <- lapply(metrics, function(m) sapply(listsims, get_out_attr, m))

  # --- Construimos objeto resultado clonando el primero --------------
  res <- listsims[[1]]

  for (k in seq_along(metrics)) {
    m   <- metrics[k]
    val <- mat[[k]]
    if (is.null(nrow(val))) {
      res$out[[m]] <- list(mean = mean(val),
                           sd   = sd(val),
                           summary = summary(val))
    } else {
      res$out[[m]] <- list(mean    = apply(val, 1, mean),
                           sd      = apply(val, 1, sd),
                           summary = apply(val, 1, summary))
    }
  }

  # --- Probabilidades en estado estacionario (pn) --------------------
  first <- listsims[[1]]
  if (is.null(nrow(first$out$wq))) {
    maxlen <- max(vapply(listsims, function(x) length(x$out$pn), 1L))
    pn_mat <- vapply(listsims, function(x) {
      len <- length(x$out$pn)
      c(x$out$pn, rep(NA_real_, maxlen - len))
    }, numeric(maxlen))
    res$out$pn <- rowMeans(pn_mat, na.rm = TRUE)
  } else {
    dims   <- vapply(listsims, function(x) dim(x$out$pn)[1], 1L)
    maxlen <- max(dims)
    s_dim  <- length(first$s)          # nº de servidores del nodo (para Redes)
    arr    <- array(NA_real_,
                    dim = c(maxlen, s_dim, length(listsims)))
    for (i in seq_along(listsims)) {
      pn_i <- listsims[[i]]$out$pn
      len  <- dim(pn_i)[1]
      arr[seq_len(len), , i] <- pn_i
    }
    res$out$pn <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  }

  res
}

#' Ejecuta en paralelo `nsim` simulaciones de un modelo
#'
#' @param modelfunction  Función que genera **una** réplica (debe retornar
#'   un objeto con campo `$out`)
#' @param parameters     `list` con los argumentos que recibe `modelfunction`
#' @param nsim           Número de réplicas a lanzar
#' @param nproc          Núcleos a utilizar (>= 1).  Si vale 1 → secuencial
#' @return Una lista con las réplicas, o la única réplica si `nsim == 1`
#' @export
ParallelizeSimulations <- function(modelfunction,
                                   parameters,
                                   nsim  = 1L,
                                   nproc = 1L) {

  stopifnot(nsim  >= 1, nproc >= 1)

  # --- Caso paralelo -------------------------------------------------- #
  if (nproc > 1L) {
    cl <- parallel::makeCluster(nproc)
    doParallel::registerDoParallel(cl)

    res <- foreach::foreach(i = seq_len(nsim),
                            .packages = c("distr")) %dopar% {
                              do.call(modelfunction, parameters)
                            }

    parallel::stopCluster(cl)
  } else {
    # --- Caso secuencial ---------------------------------------------- #
    res <- lapply(seq_len(nsim), function(i) do.call(modelfunction, parameters))
  }

  if (length(res) == 1L) res[[1]] else res
}
