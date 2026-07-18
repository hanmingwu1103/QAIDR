# E10b published tie-impact audit library (Prompt E, 2026-07-17).
#
# For every tie-affected manuscript cell, classifies each detected tolerance
# tie block (in either rank space) as an exact structural duplicate versus a
# non-exact tolerance near-tie, and compares the production (pre-patch) tie
# resolution against the corrected (post-patch) resolution over the original
# 50 seeded resolutions. The pre-patch ranking routine is VENDORED here for
# the audit only; production code uses the corrected rank_matrix().
suppressMessages(library(QAIDR))

E10B_TOL <- 1e-12

# --- vendored pre-patch rank_matrix (QAIDR 0.2.0 before Prompt E) ----------
rank_matrix_prepatch <- function(D, ties = c("error", "random"), tol = 1e-12) {
  ties <- match.arg(ties)
  N <- nrow(D)
  R <- matrix(0L, N, N)
  n_ties <- 0L
  for (i in seq_len(N)) {
    d <- D[i, -i]
    ds <- sort(d)
    gap_scale <- pmax(ds[-1], .Machine$double.xmin)
    tie_pairs <- sum(diff(ds) <= tol * gap_scale)
    if (tie_pairs > 0L) {
      n_ties <- n_ties + tie_pairs
      if (ties == "error") stop("tie")
    }
    R[i, -i] <- as.integer(rank(d, ties.method = if (ties == "random") "random" else "first"))
  }
  attr(R, "n_ties") <- n_ties
  R
}

coranking_indices_prepatch <- function(Dh, Dl, K) {
  Rh <- rank_matrix_prepatch(Dh, ties = "random")
  Rl <- rank_matrix_prepatch(Dl, ties = "random")
  indices_from_coranking(coranking_matrix(Rh, Rl), K)
}

# 50-resolution protocol exactly as assess_cell_tie_audit() runs it.
res50 <- function(Dh, Dl, K, seed, fun) {
  vals <- matrix(NA_real_, 50L, 6L)
  for (s in seq_len(50L)) {
    set.seed(seed + s)
    vals[s, ] <- fun(Dh, Dl, K)
  }
  colnames(vals) <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
  vals
}

post_fun <- function(Dh, Dl, K) coranking_indices(Dh, Dl, K, ties = "random")

# --- block extraction (identical detection rule to rank_matrix) ------------
extract_blocks <- function(D, geom = NULL, tol = E10B_TOL) {
  n <- nrow(D); out <- list(); oi <- 0L
  for (i in seq_len(n)) {
    d <- D[i, -i]
    cand <- setdiff(seq_len(n), i)
    o <- order(d); ds <- d[o]
    gs <- pmax(ds[-1], .Machine$double.xmin)
    link <- diff(ds) <= tol * gs
    if (!any(link)) next
    bid <- cumsum(c(TRUE, !link))
    for (b in unique(bid[duplicated(bid)])) {
      mem <- which(bid == b)
      ids <- cand[o[mem]]
      vals <- ds[mem]
      gaps <- abs(outer(vals, vals, "-"))[upper.tri(matrix(0, length(mem), length(mem)))]
      dup <- FALSE
      if (!is.null(geom) && max(gaps) == 0) {
        # exact structural duplicate: some pair of tied CANDIDATES has
        # bit-identical geometry (centers and radii rows), or a candidate
        # duplicates the ANCHOR (zero distance).
        rows_eq <- function(a, b) identical(geom$C[a, ], geom$C[b, ]) &&
                                  identical(geom$R[a, ], geom$R[b, ])
        prs <- utils::combn(ids, 2)
        dup <- any(vapply(seq_len(ncol(prs)), function(k)
          rows_eq(prs[1, k], prs[2, k]), logical(1))) ||
          any(vapply(ids, function(a) rows_eq(a, i), logical(1)))
      }
      oi <- oi + 1L
      out[[oi]] <- data.frame(anchor = i, block_id = b,
        block_size = length(mem),
        candidate_ids = paste(ids, collapse = ";"),
        min_within_block_gap = min(gaps), max_within_block_gap = max(gaps),
        classification = if (max(gaps) > 0) "tolerance_near_tie"
                         else if (dup) "exact_duplicate"
                         else "exact_equal_nonduplicate")
    }
  }
  if (oi) do.call(rbind, out) else NULL
}

# Audit one flagged cell: block records for both spaces + pre/post value
# comparison over the production 50-resolution protocol.
audit_cell <- function(scen, rep, method, metric, Dh, Dl, geom_h, geom_l,
                       K, tie_seed) {
  bh <- extract_blocks(Dh, geom_h)
  bl <- extract_blocks(Dl, geom_l)
  blocks <- rbind(
    if (!is.null(bh)) cbind(space = "high", bh),
    if (!is.null(bl)) cbind(space = "low", bl))
  if (!is.null(blocks)) {
    blocks <- cbind(scenario = scen, replication = rep, Method = method,
                    Metric = metric, blocks)
  }
  pre <- res50(Dh, Dl, K, tie_seed, coranking_indices_prepatch)
  post <- res50(Dh, Dl, K, tie_seed, post_fun)
  vals <- data.frame(scenario = scen, replication = rep, Method = method,
    Metric = metric,
    t(setNames(colMeans(pre), paste0("pre_", colnames(pre)))),
    t(setNames(colMeans(post), paste0("post_", colnames(post)))),
    max_abs_mean_delta = max(abs(colMeans(pre) - colMeans(post))),
    pre_max_spread = max(apply(pre, 2, function(v) diff(range(v)))),
    post_max_spread = max(apply(post, 2, function(v) diff(range(v)))))
  list(blocks = blocks, vals = vals)
}

dl_for <- function(pr, metric) {
  if (metric == "Centers-Euclidean" || pr$type == "Point")
    as.matrix(stats::dist(pr$C)) else idist(pr$C, pr$R, metric)
}

dh_for <- function(xs, metric) {
  if (metric == "Centers-Euclidean") as.matrix(stats::dist(xs$centers))
  else idist(xs$centers, xs$radii, metric)
}
