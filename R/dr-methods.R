#' Default UMAP configuration for interval UMAP
#'
#' Returns a list of UMAP hyperparameters with default settings.
#'
#' @return A named list of UMAP configuration parameters.
#' @export
#' @examples
#' cfg <- umap_config_default()
umap_config_default <- function() {
  list(
    n_components = 2,
    n_neighbors = 10,
    min_dist = 0.01,
    spread = 1.0
  )
}


#' Strong UMAP configuration for interval UMAP
#'
#' Returns a list of UMAP hyperparameters with stronger optimization settings,
#' emphasizing local structure preservation.
#'
#' @return A named list of UMAP configuration parameters.
#' @export
#' @examples
#' cfg <- umap_config_strong()
umap_config_strong <- function() {
  list(
    n_components = 2,
    n_neighbors = 10,
    min_dist = 0.01,
    spread = 1.0,
    n_epochs = 2000,
    negative_sample_rate = 10,
    local_connectivity = 2,
    set_op_mix_ratio = 0.1,
    init = "spectral"
  )
}


#' Run interval dimensionality reduction methods
#'
#' Applies multiple interval DR methods (C-PCA, V-PCA, MR-PCA, SPCA, IMDS,
#' Int-UMAP) to an \code{interval_data} object and returns the projections.
#'
#' @param x An \code{interval_data} object (should be standardized).
#' @param methods Character vector of methods to run. Default includes all six:
#'   \code{"C-PCA"}, \code{"V-PCA"}, \code{"MR-PCA"}, \code{"SPCA"},
#'   \code{"IMDS"}, \code{"Int-UMAP"}.
#' @param labels Optional factor of class labels. If \code{NULL}, uses
#'   \code{x$labels}.
#' @param no_dim Integer number of target dimensions (default 2).
#' @param umap_config List of UMAP hyperparameters. Default uses
#'   \code{umap_config_default()}.
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#' @return An object of class \code{idr_projections}: a named list where each
#'   element contains \code{C} (centers matrix), \code{R} (radii matrix), and
#'   \code{type} ("Point" or "Interval").
#' @export
#' @examples
#' \dontrun{
#' data(cars_mm)
#' x <- standardize(cars_mm)
#' proj <- run_idr(x)
#' }
run_idr <- function(x,
                    methods = c("C-PCA", "V-PCA", "MR-PCA", "SPCA",
                                "IMDS", "Int-UMAP"),
                    labels = NULL,
                    no_dim = 2,
                    umap_config = umap_config_default(),
                    verbose = TRUE) {

  stopifnot(inherits(x, "interval_data"))

  if (is.null(labels)) labels <- x$labels

  # Build symbolicDA array
  data_MM_array <- .to_sym_array(x)
  N_OBS <- nrow(x$centers)

  # Build RSDA object if needed
  data_rsda <- NULL
  if ("Int-UMAP" %in% methods) {
    if (!requireNamespace("RSDA", quietly = TRUE)) {
      stop("Package 'RSDA' is required for Int-UMAP. ",
           "Install with install.packages('RSDA').")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      stop("Package 'dplyr' is required for Int-UMAP. ",
           "Install with install.packages('dplyr').")
    }
    tmp <- .to_mm_df(data_MM_array)
    if (!is.null(labels)) {
      tmp$class <- labels
      data_rsda <- .to_rsda(tmp, "class")
      data_rsda["class"] <- NULL
    } else {
      data_rsda <- .to_rsda(tmp)
    }
  }

  projections <- list()

  for (m in methods) {
    if (verbose) {
      cat("- Running:", m, "\n")
    }

    proj <- .run_single_dr(m, data_MM_array, data_rsda,
                           no_dim = no_dim,
                           umap_config = umap_config,
                           N_OBS = N_OBS)
    if (!is.null(proj)) {
      projections[[m]] <- proj
    }
  }

  structure(projections, class = "idr_projections")
}


#' Run a single DR method (internal)
#' @noRd
.run_single_dr <- function(method, s_arr, s_tbl,
                           no_dim = 2,
                           umap_config = umap_config_default(),
                           N_OBS = NULL) {

  if (method %in% c("C-PCA", "V-PCA", "MR-PCA", "SPCA", "IMDS")) {
    if (!requireNamespace("symbolicDA", quietly = TRUE)) {
      warning("Package 'symbolicDA' not available; skipping ", method)
      return(NULL)
    }
  }

  res <- tryCatch({
    capture.output({
      out <- switch(method,
        "C-PCA"  = symbolicDA::PCA.centers.SDA(s_arr, pc.number = no_dim),
        "V-PCA"  = symbolicDA::PCA.vertices.SDA(s_arr, pc.number = no_dim),
        "MR-PCA" = symbolicDA::PCA.mrpca.SDA(s_arr, pc.number = no_dim),
        "SPCA"   = symbolicDA::PCA.spca.SDA(s_arr, pc.number = no_dim),
        "IMDS"   = symbolicDA::interscal.SDA(s_arr, d = no_dim,
                                              calculateDist = TRUE),
        "Int-UMAP" = {
          if (!requireNamespace("RSDA", quietly = TRUE)) {
            stop("Package 'RSDA' is required for Int-UMAP.")
          }
          RSDA::sym.umap(s_tbl,
                         n_components = umap_config$n_components %||% 2,
                         n_neighbors  = umap_config$n_neighbors %||% 10,
                         min_dist     = umap_config$min_dist %||% 0.01,
                         spread       = umap_config$spread %||% 1.0)
        }
      )
    })
    out
  }, error = function(e) {
    warning("DR method '", method, "' failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(res)) return(NULL)

  # Post-process Int-UMAP
  if (method == "Int-UMAP") {
    l <- length(attr(res, "names_umap"))
    res$group <- sort(rep(seq_len(l), (nrow(res) / l)))
    x2 <- stats::aggregate(
      cbind(V1, V2) ~ group, data = res,
      FUN = function(v) c(min = min(v), max = max(v))
    )
    n_obj <- if (!is.null(N_OBS)) N_OBS else l
    res_arr <- array(0, dim = c(n_obj, 2, 2))
    res_arr[, , 1] <- cbind(x2$V1[, "min"], x2$V2[, "min"])
    res_arr[, , 2] <- cbind(x2$V1[, "max"], x2$V2[, "max"])
    dimnames(res_arr)[[3]] <- c("Min", "Max")
    res <- list(layout = res_arr)
  }

  # Extract projection
  p <- if (is.list(res)) {
    if (!is.null(res$ex)) res$ex
    else if (!is.null(res$xprim)) res$xprim
    else if (!is.null(res$pc)) res$pc
    else if (!is.null(res$layout)) res$layout
    else res
  } else {
    res
  }

  # Classify output
  if (length(dim(p)) == 3) {
    list(C = (p[, , 1] + p[, , 2]) / 2,
         R = (p[, , 2] - p[, , 1]) / 2,
         type = "Interval")
  } else if (method == "V-PCA" && is.matrix(p) &&
             !is.null(N_OBS) && nrow(p) > N_OBS) {
    nv <- nrow(p) / N_OBS
    Cv <- matrix(0, N_OBS, no_dim)
    Rv <- matrix(0, N_OBS, no_dim)
    for (i in seq_len(N_OBS)) {
      idx <- ((i - 1) * nv + 1):(i * nv)
      pts <- p[idx, , drop = FALSE]
      mn <- apply(pts, 2, min)
      mx <- apply(pts, 2, max)
      Cv[i, ] <- (mn + mx) / 2
      Rv[i, ] <- (mx - mn) / 2
    }
    list(C = Cv, R = Rv, type = "Interval")
  } else {
    list(C = as.matrix(p),
         R = matrix(0, nrow(p), no_dim),
         type = "Point")
  }
}


#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' @export
print.idr_projections <- function(x, ...) {
  methods <- names(x)
  cat(sprintf("idr_projections: %d methods\n", length(methods)))
  for (m in methods) {
    n <- nrow(x[[m]]$C)
    cat(sprintf("  %s: %d obs, type = %s\n", m, n, x[[m]]$type))
  }
  invisible(x)
}
