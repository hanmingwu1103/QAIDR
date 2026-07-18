# Data-generating mechanisms for the three manuscript simulation scenarios.
# Regenerated designs are identical to the manuscript's stated designs
# (tex Section 5); only the execution protocol (seeding, replication,
# corrected indices) changes. Sourced after _common.R.

#' Scenario I: three groups separated only by interval width (n=300, p=5).
gen_scenario1 <- function(seed) {
  set.seed(seed)
  N <- 300L; P <- 5L
  centers <- matrix(rnorm(N * P), N, P)
  radii <- rbind(
    matrix(runif(100 * P, 0.1, 0.2), 100, P),
    matrix(runif(100 * P, 1.0, 1.2), 100, P),
    matrix(runif(100 * P, 5.0, 6.0), 100, P)
  )
  colnames(centers) <- colnames(radii) <- paste0("V", seq_len(P))
  interval_data(centers, radii,
                labels = factor(rep(c("small", "medium", "large"), each = 100)))
}

#' Scenario II: interval Swiss roll with variance-dominated y (n=800, p=3).
gen_scenario2 <- function(seed, n = 800L) {
  set.seed(seed)
  t <- sort(runif(n, 1.5 * pi, 4.5 * pi))
  x <- t * cos(t); z <- t * sin(t)
  y <- runif(n, 0, 70)
  centers <- cbind(x = x, y = y, z = z)
  radii <- matrix(runif(n * 3, 0.2, 0.5), n, 3,
                  dimnames = list(NULL, c("x", "y", "z")))
  interval_data(centers, radii)
}

#' Scenario III: dense overlapping regime (n=200, p=5).
gen_scenario3 <- function(seed) {
  set.seed(seed)
  N <- 200L; P <- 5L
  centers <- matrix(rnorm(N * P, 0, 0.05), N, P)
  radii <- matrix(runif(N * P, 1.2, 2.5), N, P)
  colnames(centers) <- colnames(radii) <- paste0("V", seq_len(P))
  interval_data(centers, radii)
}

#' Scenario IV (B2, preregistered A-D2): heterogeneous/adversarial stress.
#' Order of operations and all constants exactly as frozen in
#' B2_PREREGISTRATION.md (SHA-256 2c25e3b2...). Returns interval_data plus
#' the per-replication record sets as attributes.
gen_scenario4 <- function(seed) {
  set.seed(seed)
  N <- 300L; P <- 5L
  sizes <- c(200L, 80L, 20L)
  donors <- c(1L, 3L, 5L, 7L, 9L, 11L)
  recipients <- c(2L, 4L, 6L, 8L, 10L, 12L)
  g <- rep(1:3, times = sizes)
  mu <- matrix(rnorm(3 * P, 0, 2), 3, P)
  centers <- mu[g, ] + matrix(rnorm(N * P), N, P)
  dist_core <- sqrt(rowSums((centers - mu[g, ])^2))
  u <- matrix(runif(N * P, 0.8, 1.2), N, P)
  radii <- 0.15 * dist_core * u                       # multiplicative only
  eligible <- setdiff(13:N, integer(0))
  contam <- sort(sample(eligible, 15L))
  centers[contam, ] <- matrix(runif(15L * P, -3, 3), 15L, P)
  radii[contam, ] <- radii[contam, ] * 2
  degen_pool <- setdiff(eligible, contam)
  degen <- sort(sample(degen_pool, 30L))
  radii[degen, ] <- 0
  M <- matrix(rnorm(P * P), P, P)
  qr_ <- qr(M)
  Qr <- qr.Q(qr_) %*% diag(sign(diag(qr.R(qr_))))     # seeded Haar rotation
  centers <- centers %*% Qr                            # centers ONLY (declared
                                                       # misspecified representation)
  centers[recipients, ] <- centers[donors, ]           # duplication LAST
  radii[recipients, ] <- radii[donors, ]
  colnames(centers) <- colnames(radii) <- paste0("V", seq_len(P))
  x <- interval_data(centers, radii, labels = factor(g))
  attr(x, "contaminated") <- contam
  attr(x, "degenerate") <- degen
  attr(x, "donors") <- donors
  attr(x, "recipients") <- recipients
  x
}


#' Run all DR methods on (optionally standardized) data, timing each method.
#' Returns list(projections, timings, failures).
run_methods_timed <- function(x, methods = METHODS, umap_config = umap_config_default()) {
  proj <- list(); times <- numeric(0); fails <- character(0)
  for (m in methods) {
    t0 <- proc.time()[3]
    p <- tryCatch(
      run_idr(x, methods = m, umap_config = umap_config, verbose = FALSE)[[m]],
      error = function(e) { fails <<- c(fails, paste0(m, ": ", conditionMessage(e))); NULL })
    times[m] <- proc.time()[3] - t0
    if (!is.null(p)) proj[[m]] <- p
  }
  list(projections = structure(proj, class = "idr_projections"),
       timings = times, failures = fails)
}
