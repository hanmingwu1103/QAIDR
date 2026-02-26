#' Assess quality of interval DR projections
#'
#' Computes quality and behavior indices for each combination of DR method
#' and interval distance metric. Optionally performs permutation tests for
#' statistical significance.
#'
#' @param x An \code{interval_data} object (standardized).
#' @param projections An \code{idr_projections} object from \code{run_idr()},
#'   or a named list with the same structure.
#' @param K Integer neighborhood size (default 5).
#' @param metrics Character vector of distance metrics to evaluate.
#' @param perm_test Logical; perform permutation tests (default \code{FALSE}).
#' @param n_perm Integer number of permutations (default 1000).
#' @return An object of class \code{qaidr_assessment} containing:
#'   \describe{
#'     \item{results}{Data frame with columns IDR, Metric, Q_TC, B_TC,
#'       Q_RE, B_RE, Q_LC, B_LC.}
#'     \item{pvalues}{Data frame of p-values (if \code{perm_test = TRUE}).}
#'     \item{K}{The neighborhood size used.}
#'   }
#' @export
#' @examples
#' \dontrun{
#' data(cars_mm)
#' x <- standardize(cars_mm)
#' proj <- run_idr(x)
#' result <- assess_quality(x, proj, K = 5, perm_test = TRUE)
#' print(result)
#' }
assess_quality <- function(x,
                           projections,
                           K = 5,
                           metrics = c("Int-Euclidean", "Hausdorff",
                                       "Ichino-Yaguchi", "Wasserstein"),
                           perm_test = FALSE,
                           n_perm = 1000) {

  if (inherits(x, "interval_data")) {
    Centers <- x$centers
    Radii <- x$radii
  } else {
    stop("x must be an interval_data object")
  }

  N_OBS <- nrow(Centers)
  results <- data.frame()
  pvalues <- data.frame()

  for (m in names(projections)) {
    out <- projections[[m]]
    for (met in metrics) {
      Dh <- idist(Centers, Radii, met)

      Dl <- if (out$type == "Point") {
        as.matrix(dist(out$C))
      } else {
        idist(out$C, out$R, met)
      }

      obs <- coranking_indices(Dh, Dl, K)

      if (perm_test) {
        cnt <- rep(0, 6)
        for (i in seq_len(n_perm)) {
          idx <- sample(N_OBS)
          null_val <- coranking_indices(Dh, Dl[idx, idx], K)
          cnt[1] <- cnt[1] + (null_val[1] >= obs[1])
          cnt[2] <- cnt[2] + (abs(null_val[2]) >= abs(obs[2]))
          cnt[3] <- cnt[3] + (null_val[3] >= obs[3])
          cnt[4] <- cnt[4] + (abs(null_val[4]) >= abs(obs[4]))
          cnt[5] <- cnt[5] + (null_val[5] >= obs[5])
          cnt[6] <- cnt[6] + (abs(null_val[6]) >= abs(obs[6]))
        }
        p <- (cnt + 1) / (n_perm + 1)
      } else {
        p <- rep(NA_real_, 6)
      }

      results <- rbind(results, data.frame(
        IDR = m, Metric = met,
        Q_TC = obs[1], B_TC = obs[2],
        Q_RE = obs[3], B_RE = obs[4],
        Q_LC = obs[5], B_LC = obs[6],
        stringsAsFactors = FALSE
      ))

      pvalues <- rbind(pvalues, data.frame(
        IDR = m, Metric = met,
        p_Q_TC = p[1], p_B_TC = p[2],
        p_Q_RE = p[3], p_B_RE = p[4],
        p_Q_LC = p[5], p_B_LC = p[6],
        stringsAsFactors = FALSE
      ))
    }
  }

  rownames(results) <- NULL
  rownames(pvalues) <- NULL

  structure(
    list(results = results,
         pvalues = if (perm_test) pvalues else NULL,
         K = K),
    class = "qaidr_assessment"
  )
}


#' Permutation test for a single DR projection
#'
#' Performs a permutation test to assess statistical significance of
#' quality and behavior indices.
#'
#' @param D_high High-dimensional distance matrix (n x n).
#' @param D_low Low-dimensional distance matrix (n x n).
#' @param K Integer neighborhood size.
#' @param n_perm Number of permutations (default 1000).
#' @return A list with elements:
#'   \describe{
#'     \item{vals}{Named vector of observed index values.}
#'     \item{pQ}{P-values for quality indices (one-tailed).}
#'     \item{pB}{P-values for behavior indices (two-tailed).}
#'   }
#' @export
#' @examples
#' set.seed(42)
#' Dh <- as.matrix(dist(matrix(rnorm(50), 10, 5)))
#' Dl <- as.matrix(dist(matrix(rnorm(20), 10, 2)))
#' pt <- perm_test(Dh, Dl, K = 3, n_perm = 99)
perm_test <- function(D_high, D_low, K, n_perm = 1000) {
  obs <- coranking_indices(D_high, D_low, K)
  null_stats <- matrix(0, n_perm, 6)
  N <- nrow(D_low)

  for (i in seq_len(n_perm)) {
    perm_idx <- sample(N)
    D_null <- D_low[perm_idx, perm_idx]
    null_stats[i, ] <- coranking_indices(D_high, D_null, K)
  }

  pQ <- c(
    TC = mean(null_stats[, 1] >= obs["Q_TC"]),
    RE = mean(null_stats[, 3] >= obs["Q_RE"]),
    LC = mean(null_stats[, 5] >= obs["Q_LC"])
  )
  pB <- c(
    TC = mean(abs(null_stats[, 2]) >= abs(obs["B_TC"])),
    RE = mean(abs(null_stats[, 4]) >= abs(obs["B_RE"])),
    LC = mean(abs(null_stats[, 6]) >= abs(obs["B_LC"]))
  )

  list(vals = obs, pQ = pQ, pB = pB)
}


#' Compute quality/behavior index profiles over K
#'
#' Computes co-ranking indices for all neighborhood sizes from 1 to
#' \code{K_max}, for each combination of DR method and distance metric.
#'
#' @param x An \code{interval_data} object (standardized).
#' @param projections An \code{idr_projections} object.
#' @param K_max Maximum neighborhood size (default \code{nrow(x$centers) - 2}).
#' @param metrics Character vector of distance metrics.
#' @return A data frame with columns Method, Metric, K, Q_TC, B_TC,
#'   Q_RE, B_RE, Q_LC, B_LC.
#' @export
#' @examples
#' \dontrun{
#' data(cars_mm)
#' x <- standardize(cars_mm)
#' proj <- run_idr(x)
#' profiles <- k_profiles(x, proj, K_max = 10)
#' }
k_profiles <- function(x,
                       projections,
                       K_max = NULL,
                       metrics = c("Int-Euclidean", "Hausdorff",
                                   "Ichino-Yaguchi", "Wasserstein")) {

  if (inherits(x, "interval_data")) {
    Centers <- x$centers
    Radii <- x$radii
  } else {
    stop("x must be an interval_data object")
  }

  N_OBS <- nrow(Centers)
  if (is.null(K_max)) K_max <- N_OBS - 2

  plot_data <- data.frame()

  for (met in metrics) {
    Dh <- idist(Centers, Radii, met)

    for (m in names(projections)) {
      out <- projections[[m]]
      Dl <- if (out$type == "Point") {
        as.matrix(dist(out$C))
      } else {
        idist(out$C, out$R, met)
      }

      for (k in seq_len(K_max)) {
        val <- coranking_indices(Dh, Dl, k)
        plot_data <- rbind(
          plot_data,
          data.frame(Method = m, Metric = met, K = k,
                     Q_TC = val["Q_TC"], B_TC = val["B_TC"],
                     Q_RE = val["Q_RE"], B_RE = val["B_RE"],
                     Q_LC = val["Q_LC"], B_LC = val["B_LC"],
                     stringsAsFactors = FALSE)
        )
      }
    }
  }

  rownames(plot_data) <- NULL
  plot_data
}


#' @export
print.qaidr_assessment <- function(x, ...) {
  cat(sprintf("QAIDR Assessment (K = %d)\n\n", x$K))

  res <- x$results
  if (!is.null(x$pvalues)) {
    # Merge stars into display
    pv <- x$pvalues
    for (col in c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")) {
      pcol <- paste0("p_", col)
      res[[col]] <- mapply(
        .fmt_pval, x$results[[col]], pv[[pcol]],
        MoreArgs = list(alpha = 0.05, digits = 3)
      )
    }
  } else {
    for (col in c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")) {
      res[[col]] <- sprintf("%.3f", x$results[[col]])
    }
  }

  print(res, row.names = FALSE)
  invisible(x)
}


#' @export
summary.qaidr_assessment <- function(object, ...) {
  print.qaidr_assessment(object, ...)
}
