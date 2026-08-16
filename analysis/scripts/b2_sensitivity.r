# B2 lambda/nu sensitivity (preregistered A-D3, frozen 2026-07-16).
# Face arm: EXACT breakpoint enumeration (union of high-space and per-method
# low-space exceptional values) for lambda on [0,1] (Int-Euclidean) and nu on
# [0,0.5] (Ichino-Yaguchi), all six methods, K=5. Scenario I arm: C-PCA only,
# frozen grids (101 pts lambda / 51 pts nu), DESCRIPTIVE lower bound.
#
# Prompt E (2026-07-17) exceptional-value contract (three-provider agreed):
# an exceptional value is ANY zero in the closed domain of a
# non-identically-zero comparison difference (affine in lambda; degree <= 2
# in nu, a zero-discriminant double root/tangency counting once), in any
# parameter-dependent rank system. Classification:
#   internal   - zeros in the open domain (published n_breakpoints = count of
#                deduplicated internal union keys; unchanged convention);
#   boundary   - zeros at a domain endpoint (metadata; audited separately);
#   tangency   - flagged subset of internal (double roots);
#   persistent - identically-zero comparisons (recorded, never a root value).
# Published index RANGES use one representative per tie-free open cell plus
# NONexceptional domain endpoints only; exceptional endpoints are evaluated
# under the documented 50-resolution tie audit and reported separately.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("b2_sensitivity")

idxn <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
TOL_ROOT <- 1e-10; TOL_IDX <- 1e-9
TOL_DISC <- 128 * .Machine$double.eps
TOL_ENDPOINT <- 8 * TOL_DISC

source(file.path(SCRIPT_DIR, "b2_breaks_lib.r"))

eval_profile <- function(xs, proj, K, par, values, metric) {
  out <- matrix(NA_real_, length(values), 6, dimnames = list(NULL, idxn))
  rk_hash <- character(length(values))
  tie_aff <- logical(length(values)); tie_spr <- numeric(length(values))
  for (v in seq_along(values)) {
    lam <- if (par == "lambda") values[v] else 0.5
    nu <- if (par == "nu") values[v] else 0.5
    Dh <- idist(xs$centers, xs$radii, metric, lambda = lam, nu = nu)
    Dl <- if (proj$type == "Point") as.matrix(stats::dist(proj$C))
          else idist(proj$C, proj$R, metric, lambda = lam, nu = nu)
    cell <- assess_cell_tie_audit(Dh, Dl, K, seed = 1L)
    out[v, ] <- cell$vals
    tie_aff[v] <- cell$tie_affected; tie_spr[v] <- cell$max_spread
    rk <- tryCatch(digest::digest(list(rank_matrix(Dh), rank_matrix(Dl))),
                   error = function(e) NA_character_)
    rk_hash[v] <- rk
  }
  list(vals = out, hash = rk_hash, tie_affected = tie_aff, tie_spread = tie_spr)
}

## ---------------- Face arm (exact) -----------------------------------------
if (!requireNamespace("dataSDA", quietly = TRUE)) {
  stop("Package 'dataSDA' is required for the Face sensitivity analysis.")
}
utils::data("face.iGAP", package = "dataSDA", envir = environment())
face_labels <- sub("[[:digit:]]+$", "", rownames(face.iGAP))
face_data <- interval_data_from_dataSDA(face.iGAP, labels = face_labels)
xf <- standardize(face_data)
set.seed(SEED_BASE + 4000L)
rmf <- run_methods_timed(xf)
face_rows <- list(); fr <- 0L
exc_summary <- list(); es <- 0L
exc_records <- list()
endpoint_audit <- list(); ea <- 0L
for (par in c("lambda", "nu")) {
  metric <- if (par == "lambda") "Int-Euclidean" else "Ichino-Yaguchi"
  dom_hi <- if (par == "lambda") 1 else 0.5
  brk_fun <- if (par == "lambda") lambda_breaks else nu_breaks
  br_hi <- brk_fun(xf$centers, xf$radii, "high")
  for (m in names(rmf$projections)) {
    pr <- rmf$projections[[m]]
    br_lo <- if (pr$type == "Point") NULL else brk_fun(pr$C, pr$R, "low")
    br_all <- if (is.null(br_lo)) br_hi else union_breaks(br_hi, br_lo)
    br <- br_all$internal               # published count convention
    cuts <- c(0, br, dom_hi)
    mids <- (head(cuts, -1) + tail(cuts, -1)) / 2
    pts <- sort(unique(pmin(pmax(c(0, mids, dom_hi), 0), dom_hi)))
    ev <- eval_profile(xf, pr, K = 5, par, pts, metric)

    ## endpoint classification: exceptional iff a boundary record hits it
    ep_exc <- vapply(c(0, dom_hi), function(e)
      any(abs(br_all$boundary - e) <= TOL_ENDPOINT * max(1, abs(e))),
      logical(1))
    in_range <- rep(TRUE, length(pts))
    for (k in 1:2) {
      e <- c(0, dom_hi)[k]
      if (ep_exc[k]) in_range[which(abs(pts - e) < .Machine$double.eps)] <- FALSE
    }
    rng <- apply(ev$vals[in_range, , drop = FALSE], 2,
                 function(v) diff(range(v)))
    rng_all_pts <- apply(ev$vals, 2, function(v) diff(range(v)))

    ## nonexceptional endpoints must agree with their adjacent cell value
    for (k in 1:2) {
      e <- c(0, dom_hi)[k]
      vi <- which(abs(pts - e) < .Machine$double.eps)
      adj <- if (k == 1) 2L else length(pts) - 1L
      ea <- ea + 1L
      endpoint_audit[[ea]] <- data.frame(arm = "Face", par = par, Method = m,
        endpoint = e, exceptional = ep_exc[k],
        tie_affected = ev$tie_affected[vi], tie_spread = ev$tie_spread[vi],
        max_abs_diff_vs_adjacent_cell =
          if (!ep_exc[k]) max(abs(ev$vals[vi, ] - ev$vals[adj, ])) else NA_real_,
        t(setNames(ev$vals[vi, ], paste0("val_", idxn))))
    }

    fr <- fr + 1L
    face_rows[[fr]] <- data.frame(arm = "Face", par = par, Method = m,
                                  n_breakpoints = length(br),
                                  n_distinct_rankings = length(unique(na.omit(ev$hash))),
                                  n_distinct_idx = nrow(unique(round(ev$vals / TOL_IDX))),
                                  t(setNames(rng, paste0("range_", idxn))))
    es <- es + 1L
    exc_summary[[es]] <- data.frame(arm = "Face", par = par, Method = m,
      n_internal_keys = length(br),
      n_boundary_values = length(br_all$boundary),
      n_tangency_keys = length(br_all$tangency),
      n_persistent_comparisons = nrow(br_all$persistent),
      range_convention = "open-cell midpoints + nonexceptional endpoints",
      max_range_all_pts_minus_published = max(rng_all_pts - rng),
      tol_root = TOL_ROOT, tol_disc = TOL_DISC)
    exc_records[[paste(par, m, sep = "_")]] <-
      list(internal_meta = br_all$internal_meta,
           boundary = br_all$boundary, tangency = br_all$tangency,
           persistent = br_all$persistent,
           records = br_all$records)
    cat(sprintf("Face %s %s: %d breakpoints (%d boundary, %d tangency, %d persistent), %d rankings, max range %.4f\n",
                par, m, length(br), length(br_all$boundary),
                length(br_all$tangency), nrow(br_all$persistent),
                length(unique(na.omit(ev$hash))), max(rng)))
  }
}

## ---------------- Scenario I arm (frozen grid, descriptive) ----------------
x1 <- standardize(gen_scenario1(SEED_BASE + 1001L))
set.seed(SEED_BASE + 1001L + 500000L)
rm1 <- run_methods_timed(x1, methods = "C-PCA")
sim_rows <- list(); sr <- 0L
for (par in c("lambda", "nu")) {
  metric <- if (par == "lambda") "Int-Euclidean" else "Ichino-Yaguchi"
  grid <- if (par == "lambda") seq(0, 1, length.out = 101) else seq(0, 0.5, length.out = 51)
  ev <- eval_profile(x1, rm1$projections[["C-PCA"]], K = 10, par, grid, metric)
  rng <- apply(ev$vals, 2, function(v) diff(range(v)))
  sr <- sr + 1L
  sim_rows[[sr]] <- data.frame(arm = "ScenarioI-C-PCA", par = par, Method = "C-PCA",
                               n_breakpoints = NA, n_distinct_rankings =
                                 length(unique(na.omit(ev$hash))),
                               n_distinct_idx = nrow(unique(round(ev$vals / TOL_IDX))),
                               t(setNames(rng, paste0("range_", idxn))))
  saveRDS(list(grid = grid, vals = ev$vals),
          file.path(OUTPUT_DIR, paste0("b2_sensitivity_grid_", par, ".rds")))
}
out <- rbind(do.call(rbind, face_rows), do.call(rbind, sim_rows))
save_result(out, "b2_sensitivity",
            extra = list(prereg_sha256 = "2c25e3b2d39ee824043fdf5261e299f13773c5a7734799f7f35ba6257f2df67a",
                         guidance_rule = "adequate iff max range over ALL methods < 0.02 on exact Face profiles; grid arm descriptive only",
                         exceptional_contract = "Prompt E 2026-07-17: internal/boundary/tangency/persistent classification; ranges over tie-free open cells + nonexceptional endpoints"))
save_result(do.call(rbind, exc_summary), "b2_sensitivity_exceptional_summary")
save_result(do.call(rbind, endpoint_audit), "b2_sensitivity_endpoint_audit")
saveRDS(exc_records, file.path(OUTPUT_DIR, "b2_sensitivity_exceptional_records.rds"))
cat("B2 sensitivity done:", format(Sys.time()), "\n")
