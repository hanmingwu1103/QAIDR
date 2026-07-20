#' Assess quality of interval DR projections
#'
#' Computes the adjudicated quality and behavior indices for each combination
#' of DR method and interval dissimilarity, optionally with permutation
#' inference. P-values use the finite-sample Monte Carlo estimator
#' \code{(c + 1) / (m + 1)} throughout. When \code{perm_test = TRUE},
#' family-wise adjusted p-values are additionally computed by a symmetric
#' single-step randomization min-P construction using the shared joint
#' permutation draws: the observed row and all \code{m} draws enter one
#' \code{B = m + 1} row pool with a common denominator and inclusive
#' comparisons, giving finite-sample weak FWER control at attainable levels
#' under the complete random-correspondence null. The multiplicity family
#' is, per DR method and index family (T&C, MRRE, LCMC) at the fixed K, the
#' set of \{quality, behavior\} tests across all requested metrics. Quality
#' tests are one-sided (upper); behavior tests are two-sided.
#'
#' A center-only Euclidean baseline (\code{baseline = TRUE}) evaluates every
#' projection with plain Euclidean distances between interval centers in both
#' spaces, quantifying what interval-aware evaluation adds.
#'
#' @param x An \code{interval_data} object (standardized or raw; the caller
#'   controls preprocessing).
#' @param projections An \code{idr_projections} object from \code{run_idr()}.
#' @param K Integer neighborhood size (default 5).
#' @param metrics Character vector of interval dissimilarities.
#' @param lambda Optimism index for the Interval Euclidean score (default 0.5).
#' @param nu Ichino-Yaguchi span weight (default 0.5).
#' @param baseline Logical; add the center-only Euclidean evaluation row
#'   (default \code{TRUE}).
#' @param perm_test Logical; perform permutation inference (default
#'   \code{FALSE}).
#' @param n_perm Integer number of permutations (default 999).
#' @param ties Tie policy passed to \code{rank_matrix()} (\code{"error"} for
#'   primary analyses).
#' @param seed Optional integer seed for the permutation stream.
#' @return A \code{qaidr_assessment} object: \code{results} (indices),
#'   \code{pvalues} (unadjusted), \code{pvalues_adj} (min-P adjusted),
#'   \code{null_draws} (joint permutation draws, for reproducible
#'   multiplicity adjustment), \code{K}, \code{params}.
#' @export
assess_quality <- function(x,
                           projections,
                           K = 5,
                           metrics = c("Int-Euclidean", "Hausdorff",
                                       "Ichino-Yaguchi", "Wasserstein"),
                           lambda = 0.5,
                           nu = 0.5,
                           baseline = TRUE,
                           perm_test = FALSE,
                           n_perm = 999,
                           ties = c("error", "random"),
                           seed = NULL) {
  ties <- match.arg(ties)
  stopifnot(inherits(x, "interval_data"))
  if (!is.null(seed)) set.seed(seed)

  Centers <- x$centers
  Radii <- x$radii
  N_OBS <- nrow(Centers)
  idx_names <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

  eval_metrics <- metrics
  if (baseline) eval_metrics <- c(eval_metrics, "Centers-Euclidean")

  ## Pre-rank the high-dimensional space once per metric
  Rh_list <- list()
  for (met in eval_metrics) {
    Dh <- if (met == "Centers-Euclidean") as.matrix(stats::dist(Centers))
          else idist(Centers, Radii, met, lambda = lambda, nu = nu)
    Rh_list[[met]] <- rank_matrix(Dh, ties = ties)
  }

  results <- data.frame()
  pvalues <- data.frame()
  pvalues_adj <- data.frame()
  null_draws <- list()

  for (m_name in names(projections)) {
    out <- projections[[m_name]]
    obs_mat <- matrix(NA_real_, length(eval_metrics), 6,
                      dimnames = list(eval_metrics, idx_names))
    null_arr <- if (perm_test)
      array(NA_real_, c(n_perm, length(eval_metrics), 6),
            dimnames = list(NULL, eval_metrics, idx_names)) else NULL
    ## Joint-draw requirement: ONE set of joint permutation draws per
    ## method, applied identically to every metric in the family, so the
    ## min-P adjustment sees the true joint null dependence across metrics.
    perms <- if (perm_test)
      t(vapply(seq_len(n_perm), function(s) sample.int(N_OBS),
               integer(N_OBS))) else NULL

    for (g in seq_along(eval_metrics)) {
      met <- eval_metrics[g]
      Dl <- if (met == "Centers-Euclidean" || out$type == "Point") {
        as.matrix(stats::dist(out$C))
      } else {
        idist(out$C, out$R, met, lambda = lambda, nu = nu)
      }
      Rl <- rank_matrix(Dl, ties = ties)
      Qm <- coranking_matrix(Rh_list[[met]], Rl)
      obs_mat[g, ] <- indices_from_coranking(Qm, K)
      if (perm_test) {
        null_arr[, g, ] <- perm_null_indices(Rh_list[[met]], Rl, K,
                                             perms = perms)
      }
    }

    for (g in seq_along(eval_metrics)) {
      results <- rbind(results, data.frame(
        IDR = m_name, Metric = eval_metrics[g], t(obs_mat[g, ]),
        stringsAsFactors = FALSE))
    }

    if (perm_test) {
      null_draws[[m_name]] <- null_arr
      praw <- .p_unadjusted(obs_mat, null_arr)
      padj <- .p_minP(obs_mat, null_arr)
      for (g in seq_along(eval_metrics)) {
        pvalues <- rbind(pvalues, data.frame(
          IDR = m_name, Metric = eval_metrics[g], t(praw[g, ]),
          stringsAsFactors = FALSE))
        pvalues_adj <- rbind(pvalues_adj, data.frame(
          IDR = m_name, Metric = eval_metrics[g], t(padj[g, ]),
          stringsAsFactors = FALSE))
      }
    }
  }

  rownames(results) <- NULL
  structure(
    list(results = results,
         pvalues = if (perm_test) `rownames<-`(pvalues, NULL) else NULL,
         pvalues_adj = if (perm_test) `rownames<-`(pvalues_adj, NULL) else NULL,
         null_draws = if (perm_test) null_draws else NULL,
         K = K,
         params = list(lambda = lambda, nu = nu, n_perm = if (perm_test) n_perm else NA,
                       ties = ties, baseline = baseline, seed = seed)),
    class = "qaidr_assessment"
  )
}


#' Permutation test for a single pair of distance matrices
#'
#' Shares the permutation engine of \code{assess_quality()}: ranks are
#' computed once and permuted (label-equivariance), the null uses a single
#' uniform bijection applied to both endpoints, and p-values use
#' \code{(c + 1) / (m + 1)} so the smallest attainable p-value is
#' \code{1 / (m + 1)}.
#'
#' @param D_high High-dimensional dissimilarity matrix.
#' @param D_low Low-dimensional dissimilarity matrix.
#' @param K Integer neighborhood size.
#' @param n_perm Number of permutations (default 999).
#' @param ties Tie policy (see \code{rank_matrix()}).
#' @param seed Optional integer seed.
#' @return List with \code{vals} (observed indices), \code{pQ} (one-sided
#'   upper p-values for quality indices), \code{pB} (two-sided p-values for
#'   behavior indices), and \code{null_stats} (the m x 6 null draws).
#' @export
perm_test <- function(D_high, D_low, K, n_perm = 999,
                      ties = c("error", "random"), seed = NULL) {
  ties <- match.arg(ties)
  if (!is.null(seed)) set.seed(seed)
  Rh <- rank_matrix(D_high, ties = ties)
  Rl <- rank_matrix(D_low, ties = ties)
  obs <- indices_from_coranking(coranking_matrix(Rh, Rl), K)
  null_stats <- perm_null_indices(Rh, Rl, K, m = n_perm)

  p1 <- function(j) (1 + sum(null_stats[, j] >= obs[[j]])) / (n_perm + 1)
  p2 <- function(j) (1 + sum(abs(null_stats[, j]) >= abs(obs[[j]]))) / (n_perm + 1)
  pQ <- c(TC = p1("Q_TC"), RE = p1("Q_RE"), LC = p1("Q_LC"))
  pB <- c(TC = p2("B_TC"), RE = p2("B_RE"), LC = p2("B_LC"))
  list(vals = obs, pQ = pQ, pB = pB, null_stats = null_stats)
}


#' Compute quality/behavior index profiles over K
#'
#' Ranks each space once per metric and evaluates the indices for every K on
#' the same co-ranking matrix (the co-ranking matrix does not depend on K),
#' reducing the cost from the previous per-K re-ranking.
#'
#' @inheritParams assess_quality
#' @param K_max Maximum neighborhood size (default \code{n - 2}).
#' @return Data frame with columns Method, Metric, K and the six indices.
#' @export
k_profiles <- function(x,
                       projections,
                       K_max = NULL,
                       metrics = c("Int-Euclidean", "Hausdorff",
                                   "Ichino-Yaguchi", "Wasserstein"),
                       lambda = 0.5,
                       nu = 0.5,
                       baseline = TRUE,
                       ties = c("error", "random")) {
  ties <- match.arg(ties)
  stopifnot(inherits(x, "interval_data"))
  Centers <- x$centers; Radii <- x$radii
  N_OBS <- nrow(Centers)
  if (is.null(K_max)) K_max <- N_OBS - 2
  stopifnot(K_max >= 1, K_max <= N_OBS - 2)

  eval_metrics <- metrics
  if (baseline) eval_metrics <- c(eval_metrics, "Centers-Euclidean")

  plot_data <- vector("list", length(eval_metrics) * length(projections))
  ii <- 0L
  for (met in eval_metrics) {
    Dh <- if (met == "Centers-Euclidean") as.matrix(stats::dist(Centers))
          else idist(Centers, Radii, met, lambda = lambda, nu = nu)
    Rh <- rank_matrix(Dh, ties = ties)
    for (m_name in names(projections)) {
      out <- projections[[m_name]]
      Dl <- if (met == "Centers-Euclidean" || out$type == "Point") {
        as.matrix(stats::dist(out$C))
      } else {
        idist(out$C, out$R, met, lambda = lambda, nu = nu)
      }
      Qm <- coranking_matrix(Rh, rank_matrix(Dl, ties = ties))
      vals <- t(vapply(seq_len(K_max), function(k) indices_from_coranking(Qm, k),
                       numeric(6)))
      ii <- ii + 1L
      plot_data[[ii]] <- data.frame(Method = m_name, Metric = met,
                                    K = seq_len(K_max), vals,
                                    stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, plot_data)
  rownames(out) <- NULL
  out
}


## ---- internal: p-value machinery ------------------------------------------

#' Unadjusted Monte Carlo p-values, (c+1)/(m+1)
#' @noRd
.p_unadjusted <- function(obs_mat, null_arr) {
  m <- dim(null_arr)[1]
  p <- obs_mat; p[] <- NA_real_
  qcols <- c("Q_TC", "Q_RE", "Q_LC"); bcols <- c("B_TC", "B_RE", "B_LC")
  for (g in seq_len(nrow(obs_mat))) {
    for (j in qcols) p[g, j] <- (1 + sum(null_arr[, g, j] >= obs_mat[g, j])) / (m + 1)
    for (j in bcols) p[g, j] <- (1 + sum(abs(null_arr[, g, j]) >= abs(obs_mat[g, j]))) / (m + 1)
  }
  p
}

#' Symmetric single-step randomization min-P adjusted p-values.
#' Family: per index family (TC/RE/LC), the {Q, B} x metrics tests, using the
#' shared joint draws. The observed row and all m randomization rows use one
#' common B = m + 1 denominator and inclusive comparisons: with rows
#' b = 0..m (b = 0 observed), p[b,h] = (1/B) sum_r 1{T[r,h] >= T[b,h]},
#' M[b] = min_h p[b,h], padj[h] = (1/B) sum_b 1{M[b] <= p[0,h]}. Under row
#' exchangeability (complete random-correspondence null) this gives
#' finite-sample weak FWER control at attainable levels.
#' @noRd
.p_minP <- function(obs_mat, null_arr) {
  m <- dim(null_arr)[1]
  B <- m + 1L
  G <- nrow(obs_mat)
  fams <- list(TC = c("Q_TC", "B_TC"),
               RE = c("Q_RE", "B_RE"),
               LC = c("Q_LC", "B_LC"))
  padj <- obs_mat
  padj[] <- NA_real_

  for (fam in fams) {
    ## Symmetric randomization p-values for observed + all m joint rows.
    pall <- array(NA_real_, c(B, G, length(fam)))
    for (g in seq_len(G)) {
      for (j in seq_along(fam)) {
        v <- c(obs_mat[g, fam[j]], null_arr[, g, fam[j]])
        if (startsWith(fam[j], "B")) v <- abs(v)
        ## ties.method = "max" gives #{r: v[r] >= v[b]} inclusively.
        pall[, g, j] <- rank(-v, ties.method = "max") / B
      }
    }
    minp <- apply(pall, 1L, min)
    for (g in seq_len(G)) {
      for (j in seq_along(fam)) {
        padj[g, fam[j]] <- sum(minp <= pall[1L, g, j]) / B
      }
    }
  }
  padj
}


#' @export
print.qaidr_assessment <- function(x, ...) {
  cat(sprintf("QAIDR Assessment (K = %d%s)\n\n", x$K,
              if (!is.null(x$pvalues_adj)) ", min-P adjusted stars" else ""))
  res <- x$results
  if (!is.null(x$pvalues_adj)) {
    pv <- x$pvalues_adj
    for (col in c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")) {
      res[[col]] <- mapply(.fmt_pval, x$results[[col]], pv[[col]],
                           MoreArgs = list(alpha = 0.05, digits = 3))
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
