#' Interval Euclidean distance
#'
#' Computes the interval Euclidean distance between interval-valued
#' observations, defined as a weighted combination of the maximum and
#' minimum distances between hyperrectangles.
#'
#' @param centers Numeric matrix of interval midpoints (n x p).
#' @param radii Numeric matrix of interval half-widths (n x p).
#' @param lambda Weight parameter in \eqn{[0, 1]} (default 0.5).
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist_euclidean(C, R)
idist_euclidean <- function(centers, radii, lambda = 0.5) {
  centers <- as.matrix(centers)
  radii <- as.matrix(radii)
  N <- nrow(centers)
  D <- matrix(0, N, N)
  for (i in seq_len(N - 1)) {
    for (j in (i + 1):N) {
      d_c <- abs(centers[i, ] - centers[j, ])
      sum_r <- radii[i, ] + radii[j, ]

      d_max <- sqrt(sum((d_c + sum_r)^2))
      diff_r <- d_c - sum_r
      diff_r[diff_r < 0] <- 0
      d_min <- sqrt(sum(diff_r^2))

      D[i, j] <- D[j, i] <- lambda * d_max + (1 - lambda) * d_min
    }
  }
  D
}


#' Hausdorff distance for intervals
#'
#' Computes the L2-Hausdorff distance between interval-valued observations
#' represented as hyperrectangles.
#'
#' @param centers Numeric matrix of interval midpoints (n x p).
#' @param radii Numeric matrix of interval half-widths (n x p).
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist_hausdorff(C, R)
idist_hausdorff <- function(centers, radii) {
  centers <- as.matrix(centers)
  radii <- as.matrix(radii)
  L <- centers - radii
  U <- centers + radii
  N <- nrow(centers)
  D <- matrix(0, N, N)
  for (i in seq_len(N - 1)) {
    for (j in (i + 1):N) {
      diff_L <- abs(L[i, ] - L[j, ])
      diff_U <- abs(U[i, ] - U[j, ])
      D[i, j] <- D[j, i] <- sqrt(sum(pmax(diff_L, diff_U)^2))
    }
  }
  D
}


#' Ichino-Yaguchi dissimilarity for intervals
#'
#' Computes the Ichino-Yaguchi dissimilarity between interval-valued
#' observations based on intersection and union of intervals.
#'
#' @param centers Numeric matrix of interval midpoints (n x p).
#' @param radii Numeric matrix of interval half-widths (n x p).
#' @param gamma Weighting parameter (default 0.5).
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist_ichino_yaguchi(C, R)
idist_ichino_yaguchi <- function(centers, radii, gamma = 0.5) {
  centers <- as.matrix(centers)
  radii <- as.matrix(radii)
  L <- centers - radii
  U <- centers + radii
  N <- nrow(centers)
  P <- ncol(centers)
  D <- matrix(0, N, N)

  for (i in seq_len(N - 1)) {
    for (j in (i + 1):N) {
      phi_sum <- 0
      for (k in seq_len(P)) {
        I_min <- max(L[i, k], L[j, k])
        I_max <- min(U[i, k], U[j, k])
        len_I <- max(0, I_max - I_min)

        U_min <- min(L[i, k], L[j, k])
        U_max <- max(U[i, k], U[j, k])
        len_U <- U_max - U_min

        len_A <- U[i, k] - L[i, k]
        len_B <- U[j, k] - L[j, k]

        phi <- len_U - len_I + gamma * (2 * len_I - len_A - len_B)
        phi_sum <- phi_sum + phi^2
      }
      D[i, j] <- D[j, i] <- sqrt(phi_sum)
    }
  }
  D
}


#' L2-Wasserstein (Mallows) distance for intervals
#'
#' Computes the L2-Wasserstein distance between interval-valued observations,
#' assuming uniform distributions within each interval.
#'
#' @param centers Numeric matrix of interval midpoints (n x p).
#' @param radii Numeric matrix of interval half-widths (n x p).
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist_wasserstein(C, R)
idist_wasserstein <- function(centers, radii) {
  centers <- as.matrix(centers)
  radii <- as.matrix(radii)
  D_c <- as.matrix(dist(centers))^2
  D_r <- as.matrix(dist(radii))^2
  sqrt(D_c + (1 / 3) * D_r)
}


#' Compute interval distance matrix
#'
#' Dispatcher that computes a distance matrix using a named interval metric.
#'
#' @param centers Numeric matrix of interval midpoints (n x p), or an
#'   \code{interval_data} object.
#' @param radii Numeric matrix of interval half-widths (n x p). Ignored if
#'   \code{centers} is an \code{interval_data} object.
#' @param metric Character string: one of \code{"Int-Euclidean"},
#'   \code{"Hausdorff"}, \code{"Ichino-Yaguchi"}, or \code{"Wasserstein"}.
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist(C, R, "Wasserstein")
idist <- function(centers, radii = NULL, metric = "Wasserstein") {
  if (inherits(centers, "interval_data")) {
    radii <- centers$radii
    centers <- centers$centers
  }
  stopifnot(!is.null(radii))

  switch(metric,
    "Int-Euclidean" = idist_euclidean(centers, radii),
    "Hausdorff" = idist_hausdorff(centers, radii),
    "Ichino-Yaguchi" = idist_ichino_yaguchi(centers, radii),
    "Wasserstein" = idist_wasserstein(centers, radii),
    stop("Unknown metric: ", metric,
         ". Use one of: Int-Euclidean, Hausdorff, Ichino-Yaguchi, Wasserstein")
  )
}
