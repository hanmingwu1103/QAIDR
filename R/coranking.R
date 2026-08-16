#' Rank matrix with self-exclusion and tie detection
#'
#' Computes, for each observation, the ranks of all other observations by
#' distance, excluding the self-distance before ranking. Ranks take values in
#' \code{1, ..., N-1}. Off-diagonal ties are detected within a scale-aware
#' tolerance; under the primary (default) policy \code{ties = "error"} any tie
#' aborts with an error, so that ordinary tie-free cells are computed on
#' strict rank systems as assumed by the co-ranking theory. The opt-in
#' \code{ties = "random"} resolves each tolerance-connected tie block by a
#' uniform random ordering of the whole block (seed the session RNG for
#' reproducibility); it is used only in explicitly labelled tie audits and
#' daggered structural-tie cells, where results are summarized over seeded
#' repeated resolutions, and is never used silently as if the ranks were
#' tie-free. Randomization covers every detected block, including near-ties
#' whose members are not exactly equal (base R's random tie method alone
#' would randomize only exact equals).
#'
#' @param D Numeric symmetric distance/dissimilarity matrix (n x n),
#'   non-negative, zero diagonal.
#' @param ties Character, one of \code{"error"} (primary) or \code{"random"}
#'   (documented fallback).
#' @param tol Numeric relative tolerance used for tie detection (default
#'   \code{1e-12}): two consecutive sorted distances are a tie when their gap
#'   is at most \code{tol} times the larger of the two. The default catches
#'   structural ties (duplicate objects, coincident interval geometry) while
#'   ignoring floating-point near-misses between genuinely distinct
#'   distances.
#' @return Integer matrix (n x n) of ranks with 0 on the diagonal and each
#'   row's off-diagonal entries a permutation of \code{1:(n-1)}. Attribute
#'   \code{"n_ties"} is the total number of adjacent tolerance links in the
#'   sorted off-diagonal rows; a tolerance-connected block of size b
#'   contributes b - 1.
#' @export
#' @examples
#' D <- as.matrix(dist(matrix(rnorm(30), 10, 3)))
#' R <- rank_matrix(D)
rank_matrix <- function(D, ties = c("error", "random"), tol = 1e-12) {
  ties <- match.arg(ties)
  .check_dist(D)
  N <- nrow(D)
  R <- matrix(0L, N, N)
  n_ties <- 0L

  for (i in seq_len(N)) {
    d <- D[i, -i]                          # self-exclusion first
    o <- order(d)
    ds <- d[o]
    gap_scale <- pmax(ds[-1], .Machine$double.xmin)
    tie_link <- diff(ds) <= tol * gap_scale
    tie_pairs <- sum(tie_link)             # adjacent links, not all pairs

    if (tie_pairs > 0L) {
      n_ties <- n_ties + tie_pairs
      if (ties == "error") {
        stop("rank_matrix(): ", tie_pairs,
             " tied off-diagonal distance(s) detected for observation ", i,
             " (relative tolerance ", tol, "). Primary tie-free analyses ",
             "require strict ranks; use ties = \"random\" (seeded) only for ",
             "documented, explicitly labelled structural-tie analyses.")
      }
    }

    rk <- integer(length(d))
    rk[o] <- seq_along(d)                  # exact pre-patch tie-free ranks

    if (ties == "random" && tie_pairs > 0L) {
      bid <- cumsum(c(TRUE, !tie_link))    # tolerance-connected components
      for (b in unique(bid[duplicated(bid)])) {
        mem <- which(bid == b)             # occupied sorted rank positions
        rk[o[mem]] <- mem[sample.int(length(mem))]
      }
    }
    R[i, -i] <- rk
  }

  attr(R, "n_ties") <- n_ties
  R
}


#' Co-ranking matrix from two rank matrices
#'
#' Builds the (N-1) x (N-1) co-ranking matrix \eqn{q_{kl} = |\{(i,j): \rho_{ij}
#' = k, \gamma_{ij} = l\}|} from self-excluded rank matrices. Row index k is
#' the high-dimensional rank, column index l the low-dimensional rank.
#'
#' @param Rh Integer rank matrix for the high-dimensional space (from
#'   \code{rank_matrix()}).
#' @param Rl Integer rank matrix for the low-dimensional space.
#' @return Integer (N-1) x (N-1) matrix of class \code{"coranking"} with
#'   \code{sum(Q) == N * (N - 1)}.
#' @export
coranking_matrix <- function(Rh, Rl) {
  N <- nrow(Rh)
  stopifnot(identical(dim(Rh), dim(Rl)))
  off <- row(Rh) != col(Rh)
  k <- Rh[off]; l <- Rl[off]
  code <- (k - 1L) * (N - 1L) + l
  Q <- matrix(tabulate(code, nbins = (N - 1L)^2), N - 1L, N - 1L, byrow = TRUE)
  if (sum(Q) != N * (N - 1L)) {
    stop("coranking_matrix(): mass invariant violated (sum(Q) = ", sum(Q),
         ", expected ", N * (N - 1L), "). Check rank matrices.")
  }
  class(Q) <- c("coranking", class(Q))
  Q
}


#' Quality and behavior indices from a co-ranking matrix
#'
#' Computes the six adjudicated indices from a co-ranking matrix at
#' neighborhood size K, following Lee & Verleysen (2009, Neurocomputing
#' 72:1431-1443): Trustworthiness & Continuity (Eqs. 7-9, 26-27), Mean
#' Relative Rank Errors (Eqs. 10-12, 28-29), the Local Continuity
#' Meta-Criterion (Eq. 13), and the behavior index \eqn{B_{LC} \equiv B_{NX} =
#' U_X - U_N} (Eqs. 22-23) computed on the bounded triangles of the K x K
#' upper-left block. Positive behavior values indicate intrusion-dominated
#' embeddings for all three families.
#'
#' @param Q Co-ranking matrix from \code{coranking_matrix()}.
#' @param K Integer neighborhood size, \code{1 <= K <= N - 2}.
#' @return Named numeric vector: \code{Q_TC, B_TC, Q_RE, B_RE, Q_LC, B_LC}.
#' @export
indices_from_coranking <- function(Q, K) {
  Nm1 <- nrow(Q)
  N <- Nm1 + 1L
  stopifnot(K >= 1, K <= N - 2)
  kk <- .row_idx(Nm1); ll <- .col_idx(Nm1)

  ## --- T&C (regions LL_K / UR_K, weights (k-K) / (l-K)) ---
  G <- if (K < N / 2) N * K * (2 * N - 3 * K - 1) else N * (N - K) * (N - K - 1)
  intr <- sum((kk - K) * Q * (kk > K & ll <= K))
  extr <- sum((ll - K) * Q * (kk <= K & ll > K))
  M_T <- 1 - 2 * intr / G
  M_C <- 1 - 2 * extr / G
  Q_TC <- (M_T + M_C) / 2
  B_TC <- M_C - M_T

  ## --- MRRE (W_n over {l<=K} weight |k-l|/l; W_v over {k<=K} weight |k-l|/k) ---
  H <- N * sum(abs(N - 2 * seq_len(K) + 1) / seq_len(K))
  W_n <- sum(abs(kk - ll) / ll * Q * (ll <= K)) / H
  W_v <- sum(abs(kk - ll) / kk * Q * (kk <= K)) / H
  Q_RE <- 1 - (W_n + W_v) / 2
  B_RE <- W_n - W_v

  ## --- LCMC / B_NX (UL block and its strict triangles, both bounded by K) ---
  Q_LC <- sum(Q * (kk <= K & ll <= K)) / (N * K)
  U_N <- sum(Q * (kk <= K & ll < kk)) / (N * K)
  U_X <- sum(Q * (kk < ll & ll <= K)) / (N * K)
  B_LC <- U_X - U_N

  out <- c(Q_TC = Q_TC, B_TC = B_TC, Q_RE = Q_RE, B_RE = B_RE,
           Q_LC = Q_LC, B_LC = B_LC)
  .check_ranges(out, K)
  out
}


#' Compute co-ranking quality and behavior indices from distance matrices
#'
#' Convenience wrapper: ranks both distance matrices (self-excluded,
#' tie-checked), forms the co-ranking matrix, and evaluates the six indices.
#'
#' @inheritParams rank_matrix
#' @param Dh Numeric matrix of high-dimensional pairwise dissimilarities.
#' @param Dl Numeric matrix of low-dimensional pairwise dissimilarities.
#' @param K Integer neighborhood size.
#' @return Named numeric vector \code{Q_TC, B_TC, Q_RE, B_RE, Q_LC, B_LC},
#'   with the co-ranking matrix attached as attribute \code{"Q"} and the tie
#'   counts as attribute \code{"n_ties"}.
#' @export
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(50), 10, 5)
#' Y <- matrix(rnorm(20), 10, 2)
#' coranking_indices(as.matrix(dist(X)), as.matrix(dist(Y)), K = 3)
coranking_indices <- function(Dh, Dl, K, ties = c("error", "random"),
                              tol = 1e-12) {
  ties <- match.arg(ties)
  Rh <- rank_matrix(Dh, ties = ties, tol = tol)
  Rl <- rank_matrix(Dl, ties = ties, tol = tol)
  Q <- coranking_matrix(Rh, Rl)
  out <- indices_from_coranking(Q, K)
  attr(out, "Q") <- Q
  attr(out, "n_ties") <- attr(Rh, "n_ties") + attr(Rl, "n_ties")
  out
}


#' Permutation null distribution of the indices
#'
#' Generates the joint permutation null for the six indices under the uniform
#' random-correspondence null: a single bijection \eqn{\Pi} is drawn uniformly
#' and applied to both endpoints of the embedded ranks
#' (\eqn{\gamma^\Pi_{ij} = \gamma_{\Pi(i)\Pi(j)}}), while the high-dimensional
#' ranks stay fixed. Because ranks are label-equivariant, no re-ranking is
#' needed; each draw permutes the rows and columns of \code{Rl}.
#'
#' @param Rh,Rl Rank matrices from \code{rank_matrix()}.
#' @param K Integer neighborhood size.
#' @param m Integer number of permutations (default 999).
#' @param perms Optional m x n integer matrix whose rows are permutations of
#'   \code{1:n}. When supplied, these joint draws are used instead of fresh
#'   \code{sample.int()} draws. Required for joint-draw min-P families:
#'   every family member (e.g., every metric for one DR method) must be
#'   evaluated under the SAME permutation draws so that the joint null
#'   dependence is preserved.
#' @return An m x 6 matrix of null index values (columns as in
#'   \code{indices_from_coranking()}).
#' @export
perm_null_indices <- function(Rh, Rl, K, m = 999, perms = NULL) {
  N <- nrow(Rh)
  if (!is.null(perms)) {
    stopifnot(is.matrix(perms), ncol(perms) == N)
    m <- nrow(perms)
  }
  null_stats <- matrix(NA_real_, m, 6L)
  colnames(null_stats) <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

  ## Fast path: at fixed K every index depends only on pairs whose high-space
  ## rank is <= K or whose low-space rank is <= K (about 2NK pairs instead of
  ## N(N-1)). The low-neighbor pair set is fixed in the unpermuted
  ## coordinates: under a bijection Pi, gamma^Pi_ij = gamma_{Pi(i)Pi(j)}, so
  ## the pairs with permuted low rank = l are the preimages of the fixed
  ## pairs {(a,b): Rl[a,b] = l}. Each draw therefore only remaps two fixed
  ## pair lists. Verified against the direct co-ranking computation in the
  ## test suite.
  hi <- which(Rh <= K & Rh > 0L, arr.ind = TRUE)      # pairs with k <= K
  hi_k <- Rh[hi]
  lo <- which(Rl <= K & Rl > 0L, arr.ind = TRUE)      # pairs with l <= K (low coords)
  lo_l <- Rl[lo]

  G <- if (K < N / 2) N * K * (2 * N - 3 * K - 1) else N * (N - K) * (N - K - 1)
  H <- N * sum(abs(N - 2 * seq_len(K) + 1) / seq_len(K))

  for (s in seq_len(m)) {
    p <- if (is.null(perms)) sample.int(N) else perms[s, ]
    pinv <- integer(N); pinv[p] <- seq_len(N)

    ## low-rank pairs after permutation: (i,j) = (pinv[a], pinv[b]), rank l fixed
    i_lo <- pinv[lo[, 1L]]; j_lo <- pinv[lo[, 2L]]
    k_lo <- Rh[cbind(i_lo, j_lo)]
    ## high-rank pairs: k fixed, low rank = Rl[p(i), p(j)]
    l_hi <- Rl[cbind(p[hi[, 1L]], p[hi[, 2L]])]

    intr <- sum((k_lo - K)[k_lo > K])
    extr <- sum((l_hi - K)[l_hi > K])
    M_T <- 1 - 2 * intr / G
    M_C <- 1 - 2 * extr / G
    W_n <- sum(abs(k_lo - lo_l) / lo_l) / H
    W_v <- sum(abs(hi_k - l_hi) / hi_k) / H
    ul <- k_lo <= K
    n_ul <- sum(ul)
    U_N <- sum(ul & lo_l < k_lo) / (N * K)
    U_X <- sum(ul & lo_l > k_lo) / (N * K)

    null_stats[s, ] <- c((M_T + M_C) / 2, M_C - M_T,
                         1 - (W_n + W_v) / 2, W_n - W_v,
                         n_ul / (N * K), U_X - U_N)
  }
  null_stats
}


#' Reference (slow) permutation null via full co-ranking reconstruction
#'
#' Used to validate the fast pair-list engine in \code{perm_null_indices()}.
#' @inheritParams perm_null_indices
#' @export
perm_null_indices_reference <- function(Rh, Rl, K, m = 999) {
  N <- nrow(Rh)
  null_stats <- matrix(NA_real_, m, 6L)
  colnames(null_stats) <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  for (s in seq_len(m)) {
    p <- sample.int(N)
    Qs <- coranking_matrix(Rh, Rl[p, p, drop = FALSE])
    null_stats[s, ] <- indices_from_coranking(Qs, K)
  }
  null_stats
}


## ---- internal helpers -----------------------------------------------------

#' @noRd
.row_idx <- function(n) matrix(seq_len(n), n, n)

#' @noRd
.col_idx <- function(n) matrix(seq_len(n), n, n, byrow = TRUE)

#' @noRd
.check_dist <- function(D) {
  if (!is.matrix(D) || !is.numeric(D)) stop("distance input must be a numeric matrix")
  if (nrow(D) != ncol(D)) stop("distance matrix must be square")
  if (nrow(D) < 3) stop("need at least 3 observations")
  if (any(!is.finite(D))) stop("distance matrix contains non-finite values")
  if (any(D < 0)) stop("distance matrix contains negative values")
  if (max(abs(diag(D))) > 0) stop("distance matrix must have a zero diagonal")
  if (max(abs(D - t(D))) > 1e-8 * (1 + max(abs(D))))
    stop("distance matrix must be symmetric")
  invisible(TRUE)
}

#' @noRd
.check_ranges <- function(x, K) {
  eps <- 1e-10
  ## Q_RE is bounded above by 1 but is NOT bounded below by 0. The MRRE
  ## normalizer H_K of Lee & Verleysen (2009, Eqs. (10)-(12)) does not dominate
  ## the attainable weighted rank error when K/n is large, so Q_RE can take
  ## small negative values on entirely valid input. The normalizer is retained
  ## as published, so the lower bound is left unenforced here rather than
  ## reported as an error. Q_TC and Q_LC remain in [0, 1] unconditionally.
  rng <- rbind(Q_TC = c(0, 1), B_TC = c(-1, 1), Q_RE = c(-Inf, 1), B_RE = c(-1, 1),
               Q_LC = c(0, 1), B_LC = c(-(K - 1) / K, (K - 1) / K))
  for (nm in rownames(rng)) {
    if (x[[nm]] < rng[nm, 1] - eps || x[[nm]] > rng[nm, 2] + eps) {
      warning("index ", nm, " = ", format(x[[nm]]),
              " outside its theoretical range [", rng[nm, 1], ", ", rng[nm, 2],
              "]; this indicates an input or implementation problem ",
              "(values are NOT clamped).")
    }
  }
  if (x[["Q_RE"]] < -eps) {
    message("index Q_RE = ", format(x[["Q_RE"]]),
            " is negative; this is expected for large K/n under the retained ",
            "Lee & Verleysen (2009) normalizer and is not an input or ",
            "implementation error (values are NOT clamped).")
  }
  invisible(TRUE)
}
