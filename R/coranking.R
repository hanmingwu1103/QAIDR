#' Compute co-ranking quality and behavior indices
#'
#' Given high-dimensional and low-dimensional distance matrices, computes
#' quality (Q) and behavior (B) indices based on the co-ranking matrix
#' framework. Three index families are computed: Trustworthiness & Continuity
#' (T&C), Mean Relative Rank Errors (MRRE), and Local Continuity Meta-Criterion
#' (LCMC).
#'
#' @param Dh Numeric matrix of high-dimensional pairwise distances (n x n).
#' @param Dl Numeric matrix of low-dimensional pairwise distances (n x n).
#' @param K Integer neighborhood size.
#' @return A named numeric vector with elements \code{Q_TC}, \code{B_TC},
#'   \code{Q_RE}, \code{B_RE}, \code{Q_LC}, \code{B_LC}.
#' @export
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(50), 10, 5)
#' Y <- matrix(rnorm(20), 10, 2)
#' Dh <- as.matrix(dist(X))
#' Dl <- as.matrix(dist(Y))
#' coranking_indices(Dh, Dl, K = 3)
coranking_indices <- function(Dh, Dl, K) {
  N <- nrow(Dh)
  stopifnot(K >= 1, K <= N - 2)

  # Compute ranks
  Rh <- t(apply(Dh, 1, rank, ties.method = "first"))
  Rl <- t(apply(Dl, 1, rank, ties.method = "first"))

  # Construct co-ranking matrix Q
  Q <- matrix(0L, N - 1, N - 1)
  for (i in seq_len(N)) {
    rh <- Rh[i, -i] - 1L
    rl <- Rl[i, -i] - 1L
    valid <- rh > 0 & rh < N & rl > 0 & rl < N
    if (any(valid)) {
      idx <- cbind(rh[valid], rl[valid])
      Q[idx] <- Q[idx] + 1L
    }
  }

  near <- seq_len(K)
  far  <- (K + 1):(N - 1)

  # --- T&C ---
  intr <- 0
  for (k in far) {
    if (any(Q[k, near] > 0))
      intr <- intr + sum(Q[k, near]) * (k - K)
  }
  extr <- 0
  for (l in far) {
    if (any(Q[near, l] > 0))
      extr <- extr + sum(Q[near, l]) * (l - K)
  }
  if (K < N / 2) {
    G <- N * K * (2 * N - 3 * K - 1)
  } else {
    G <- N * (N - K) * (N - K - 1)
  }
  if (G == 0) G <- 1
  Q_TC <- ((1 - 2 * intr / G) + (1 - 2 * extr / G)) / 2
  B_TC <- (1 - 2 * extr / G) - (1 - 2 * intr / G)

  # --- MRRE ---
  Wn <- 0
  for (k in far) {
    vals <- Q[k, near]
    if (any(vals > 0)) {
      Wn <- Wn + sum(vals * abs(k - near) / near)
    }
  }
  Wv <- 0
  for (k in near) {
    vals <- Q[k, far]
    if (any(vals > 0)) {
      Wv <- Wv + sum(vals * abs(k - far) / k)
    }
  }
  kv <- seq_len(K)
  H <- N * sum(abs(N - 2 * kv + 1) / kv)
  if (H <= 0) H <- 1
  Q_RE <- 1 - (Wn + Wv) / (2 * H)
  B_RE <- (Wn - Wv) / H

  # --- LCMC ---
  overlap <- sum(Q[near, near])
  Q_LC <- overlap / (N * K)
  B_LC <- .b_lc_vectorized(Q, K)

  # Clamp
  Q_TC <- max(0, min(1, Q_TC))
  B_TC <- max(-1, min(1, B_TC))
  Q_RE <- max(0, min(1, Q_RE))
  B_RE <- max(-1, min(1, B_RE))
  Q_LC <- max(0, min(1, Q_LC))
  B_LC <- max(-1, min(1, B_LC))

  c(Q_TC = Q_TC, B_TC = B_TC,
    Q_RE = Q_RE, B_RE = B_RE,
    Q_LC = Q_LC, B_LC = B_LC)
}


#' Vectorized B_LC computation from co-ranking matrix
#' @noRd
.b_lc_vectorized <- function(Q, K) {
  n_minus_1 <- nrow(Q)
  N <- n_minus_1 + 1

  lt_mask <- matrix(FALSE, n_minus_1, n_minus_1)
  for (i in 2:K) {
    lt_mask[i, 1:(i - 1)] <- TRUE
  }

  ut_mask <- matrix(FALSE, n_minus_1, n_minus_1)
  for (i in 1:(K - 1)) {
    if ((i + 1) <= n_minus_1) {
      ut_mask[i, (i + 1):n_minus_1] <- TRUE
    }
  }

  sum_LT <- sum(Q[lt_mask])
  sum_UT <- sum(Q[ut_mask])

  U_N <- sum_LT / (N * K)
  U_X <- sum_UT / (N * K)
  U_X - U_N
}
