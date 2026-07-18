#' Interval Euclidean distance
#'
#' Computes the interval Euclidean distance between interval-valued
#' observations, defined as a weighted combination of the maximum and
#' minimum distances between hyperrectangles.
#'
#' @param centers Numeric matrix of interval midpoints (n x p).
#' @param radii Numeric matrix of interval half-widths (n x p).
#' @param lambda Weight parameter: one finite numeric scalar in \eqn{[0, 1]}
#'   (default 0.5); validated here as well as in the \code{idist()}
#'   dispatcher.
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist_euclidean(C, R)
idist_euclidean <- function(centers, radii, lambda = 0.5) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0 || lambda > 1) {
    stop("idist_euclidean(): 'lambda' must be one finite numeric scalar in [0, 1]")
  }
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


#' Product-Hausdorff distance for intervals
#'
#' Computes the coordinatewise product-Hausdorff dissimilarity between
#' interval-valued observations: the interval Hausdorff distance is taken
#' per coordinate and aggregated in \eqn{\ell_2}. This is the manuscript's
#' product-Hausdorff form; it is generally NOT the set-Hausdorff metric
#' between hyperrectangles under the Euclidean point metric. The public
#' metric key \code{"Hausdorff"} is unchanged.
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
#' observations based on the interval join (hull) and meet (intersection).
#'
#' @param centers Numeric matrix of interval midpoints (n x p).
#' @param radii Numeric matrix of interval half-widths (n x p).
#' @param nu Span weight in the canonical Ichino-Yaguchi range \eqn{[0, 0.5]}
#'   (default 0.5); named \eqn{\nu} in the manuscript.
#' @param gamma Deprecated alias for \code{nu}.
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist_ichino_yaguchi(C, R)
idist_ichino_yaguchi <- function(centers, radii, nu = 0.5, gamma = NULL) {
  if (!is.null(gamma)) {
    warning("argument 'gamma' is deprecated; use 'nu'")
    nu <- gamma
  }
  stopifnot(is.numeric(nu), length(nu) == 1, nu >= 0, nu <= 0.5)
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

        phi <- len_U - len_I + nu * (2 * len_I - len_A - len_B)
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
#' @param lambda Optimism index in \eqn{[0, 1]} for the Interval Euclidean
#'   scalarization (default 0.5); ignored by the other dissimilarities.
#' @param nu Span weight for the Ichino-Yaguchi dissimilarity, canonical
#'   range \eqn{[0, 0.5]} (default 0.5); ignored by the other dissimilarities.
#' @return A symmetric n x n distance matrix.
#' @export
#' @examples
#' C <- matrix(rnorm(12), 4, 3)
#' R <- matrix(runif(12, 0.1, 0.5), 4, 3)
#' D <- idist(C, R, "Wasserstein")
idist <- function(centers, radii = NULL, metric = "Wasserstein",
                  lambda = 0.5, nu = 0.5) {
  if (inherits(centers, "interval_data")) {
    radii <- centers$radii
    centers <- centers$centers
  }
  stopifnot(!is.null(radii))
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda >= 0, lambda <= 1)

  switch(metric,
    "Int-Euclidean" = idist_euclidean(centers, radii, lambda = lambda),
    "Hausdorff" = idist_hausdorff(centers, radii),
    "Ichino-Yaguchi" = idist_ichino_yaguchi(centers, radii, nu = nu),
    "Wasserstein" = idist_wasserstein(centers, radii),
    stop("Unknown metric: ", metric,
         ". Use one of: Int-Euclidean, Hausdorff, Ichino-Yaguchi, Wasserstein")
  )
}
