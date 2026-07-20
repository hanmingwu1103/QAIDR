# Prompt H tie-policy verification (H_PREREGISTRATION.md sec. 4 + A1.7 + A2;
# audit, not a study). On the tie-affected production cells:
#   (i)  the existing 50 uniform strict refinements (production estimand),
#   (ii) 500 refinements where runtime permits (USHCN daggered cells + the
#        flagged cells of Scenario I replications 1..10),
#   (iii) a FIXED coupled-priority strict refinement: priority = one seeded
#        uniform permutation of object indices per cell, seed
#        20260719 + 9000 + cell_index, cell_index = lexicographic order of
#        the frozen cell key (scenario, replication, method, metric); the
#        priority vector is attached to object identities and permuted WITH
#        the objects under any relabelling (A1.7); label-equivariance is
#        verified per cell by a relabelling test.
# Midranks/fractional ranks are prohibited (H.6). Material threshold: any
# displayed-value change at manuscript precision (3 decimals) or conclusion
# change -> prominent disclosure + three-provider review.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.R"))
save_sessioninfo("h_tie_policy_audit")

stopifnot(utils::packageVersion("QAIDR") >= "0.3.0")
SEED_H <- 20260719L
K_S <- 10L
pcix <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
t_start <- Sys.time()

## Priority-refined strict ranks: rank_matrix's exact tolerance-block
## detection, blocks ordered by the carried per-object priority vector.
rank_matrix_priority <- function(D, priority, tol = 1e-12) {
  N <- nrow(D)
  R <- matrix(0L, N, N)
  for (i in seq_len(N)) {
    idx <- seq_len(N)[-i]
    d <- D[i, -i]
    o <- order(d)
    ds <- d[o]
    gap_scale <- pmax(ds[-1], .Machine$double.xmin)
    tie_link <- diff(ds) <= tol * gap_scale
    rk <- integer(length(d))
    rk[o] <- seq_along(d)
    if (any(tie_link)) {
      bid <- cumsum(c(TRUE, !tie_link))
      for (b in unique(bid[duplicated(bid)])) {
        mem <- which(bid == b)
        ## strict refinement: the block member with the k-th smallest
        ## carried priority receives the block's k-th rank (direct
        ## assignment; the inverse form is not label-equivariant for
        ## blocks of size >= 3)
        pr <- priority[idx[o[mem]]]
        rk[o[mem][order(pr)]] <- mem
      }
    }
    R[i, -i] <- rk
  }
  R
}

indices_pair <- function(Dh, Dl, K, rank_fun_h, rank_fun_l) {
  Qm <- coranking_matrix(rank_fun_h(Dh), rank_fun_l(Dl))
  indices_from_coranking(Qm, K)
}

audit_cell <- function(Dh, Dl, K, cell_key, cell_index, n_500) {
  ## (i) 50 uniform refinements at the production tie seed protocol is the
  ## stored estimand; recomputed here from the same machinery for delta
  ## consistency (seeded identically to assess_cell_tie_audit usage).
  res <- list()
  for (nr in c(50L, if (n_500) 500L)) {
    vals <- matrix(NA_real_, nr, 6)
    for (s in seq_len(nr)) {
      set.seed(cell_key$tie_seed + s)
      vals[s, ] <- coranking_indices(Dh, Dl, K, ties = "random")
    }
    res[[paste0("unif", nr)]] <- colMeans(vals)
  }
  ## (iii) coupled-priority refinement
  set.seed(SEED_H + 9000L + cell_index)
  priority <- sample.int(nrow(Dh))
  rf <- function(D) rank_matrix_priority(D, priority)
  res$priority <- indices_pair(Dh, Dl, K, rf, rf)
  ## label-equivariance check: permute objects AND carried priority jointly
  set.seed(SEED_H + 990000L + cell_index)
  pmt <- sample.int(nrow(Dh))
  rf2 <- function(D) rank_matrix_priority(D, priority[pmt])
  v2 <- indices_pair(Dh[pmt, pmt], Dl[pmt, pmt], K, rf2, rf2)
  res$equivariant <- isTRUE(all.equal(unname(unlist(res$priority)),
                                      unname(unlist(v2)), tolerance = 1e-12))
  res
}

rows <- list(); rr <- 0L
add_row <- function(scenario, rep, method, metric, cell_index, res) {
  rr <<- rr + 1L
  u50 <- res$unif50
  u500 <- if (!is.null(res$unif500)) res$unif500 else rep(NA_real_, 6)
  pri <- unlist(res$priority)
  rows[[rr]] <<- data.frame(
    scenario = scenario, rep = rep, Method = method, Metric = metric,
    cell_index = cell_index,
    stats::setNames(as.list(u50), paste0(pcix, "_u50")),
    stats::setNames(as.list(u500), paste0(pcix, "_u500")),
    stats::setNames(as.list(pri), paste0(pcix, "_pri")),
    d50_500 = max(abs(u50 - u500)),
    d50_pri = max(abs(u50 - pri)),
    equivariant = res$equivariant,
    stringsAsFactors = FALSE)
}

## ---- enumerate frozen tie-affected cells (outcome-independent keys; -------
## replication subsets fixed by Amendment A3 BEFORE any audit result) --------
cells <- list()
per1 <- readRDS(file.path(OUTPUT_DIR, "scenario1_per_rep.rds"))
for (r in 1:10) {
  fl <- per1[per1$rep == r & per1$tie_affected, c("Method", "Metric")]
  if (nrow(fl)) for (i in seq_len(nrow(fl)))
    cells[[length(cells) + 1L]] <- list(scenario = "S1", rep = r,
                                        method = fl$Method[i],
                                        metric = fl$Metric[i], n500 = TRUE)
}
n_drop_s1 <- sum(per1$tie_affected & per1$rep > 10)
per2 <- readRDS(file.path(OUTPUT_DIR, "scenario2_per_rep.rds"))
for (r in intersect(1:10, unique(per2$rep[per2$tie_affected]))) {
  fl <- per2[per2$rep == r & per2$tie_affected, c("Method", "Metric")]
  for (i in seq_len(nrow(fl)))
    cells[[length(cells) + 1L]] <- list(scenario = "S2", rep = r,
                                        method = fl$Method[i],
                                        metric = fl$Metric[i], n500 = FALSE)
}
n_drop_s2 <- sum(per2$tie_affected & per2$rep > 10)
per4 <- readRDS(file.path(OUTPUT_DIR, "scenario4_per_rep.rds"))
sp4 <- grep("_spread$", names(per4), value = TRUE)
p4k <- per4[per4$K == 10 & per4$rep <= 25, ]
aff4 <- rowSums(as.matrix(p4k[sp4]) > 0) > 0
for (r in 1:5) {
  fl <- p4k[aff4 & p4k$rep == r, c("Method", "Metric")]
  for (i in seq_len(nrow(fl)))
    cells[[length(cells) + 1L]] <- list(scenario = "S4", rep = r,
                                        method = fl$Method[i],
                                        metric = fl$Metric[i], n500 = FALSE)
}
n_drop_s4 <- sum(aff4 & p4k$rep > 5)
mid <- readRDS(file.path(OUTPUT_DIR, "midscale_indices.rds"))
midf <- mid[mid$tie_affected, c("Method", "Metric")]
for (i in seq_len(nrow(midf)))
  cells[[length(cells) + 1L]] <- list(scenario = "USHCN", rep = 0L,
                                      method = midf$Method[i],
                                      metric = midf$Metric[i], n500 = TRUE)
cat(sprintf("A3 scope: audited later-replication flagged cells dropped: S1=%d S2=%d S4=%d\n",
            n_drop_s1, n_drop_s2, n_drop_s4))

## frozen lexicographic cell_index over (scenario, rep, method, metric)
keydf <- do.call(rbind, lapply(cells, function(z)
  data.frame(scenario = z$scenario, rep = z$rep, method = z$method,
             metric = z$metric, stringsAsFactors = FALSE)))
ord <- order(keydf$scenario, keydf$rep, keydf$method, keydf$metric)
cells <- cells[ord]
cat("tie-affected cells:", length(cells), "\n")
print(keydf[ord, ], row.names = FALSE)

## ---- recompute per scenario (production streams) --------------------------
dist_for <- function(xs, pr, met) {
  if (met == "Centers-Euclidean" || pr$type == "Point")
    as.matrix(stats::dist(pr$C)) else idist(pr$C, pr$R, met)
}
by_ctx <- split(seq_along(cells), vapply(cells, function(z)
  paste(z$scenario, z$rep), ""))
for (ctx in names(by_ctx)) {
  ids <- by_ctx[[ctx]]
  z1 <- cells[[ids[1]]]
  if (z1$scenario == "S1") {
    seed_r <- SEED_BASE + 1000L + z1$rep
    x <- gen_scenario1(seed_r); xs <- standardize(x)
    set.seed(seed_r + 500000L); rm <- run_methods_timed(xs)
    tie_seed <- seed_r; K <- K_S
  } else if (z1$scenario == "S2") {
    seed_r <- SEED_BASE + 2000L + z1$rep
    x <- gen_scenario2(seed_r); xs <- standardize(x)
    set.seed(seed_r + 500000L); rm <- run_methods_timed(x)  # raw_dr per production
    tie_seed <- seed_r; K <- K_S
  } else if (z1$scenario == "S4") {
    seed_r <- SEED_BASE + 7000L + z1$rep
    x <- gen_scenario4(seed_r); xs <- standardize(x)
    set.seed(seed_r + 500000L); rm <- run_methods_timed(xs)
    tie_seed <- seed_r + 300000L; K <- K_S
  } else {
    snap <- readRDS(file.path(DATA_PROC, "ghcn_seasonal_intervals.rds"))
    x <- interval_data(snap$centers, snap$radii); xs <- standardize(x)
    set.seed(SEED_BASE + 6000L); rm <- run_methods_timed(xs)
    tie_seed <- SEED_BASE + 6000L; K <- 10L
  }
  for (id in ids) {
    z <- cells[[id]]
    met <- z$metric
    Dh <- if (met == "Centers-Euclidean") as.matrix(stats::dist(xs$centers))
          else idist(xs$centers, xs$radii, met)
    Dl <- dist_for(xs, rm$projections[[z$method]], met)
    res <- audit_cell(Dh, Dl, K,
                      cell_key = list(tie_seed = tie_seed),
                      cell_index = id, n_500 = z$n500)
    add_row(z$scenario, z$rep, z$method, z$metric, id, res)
    cat(sprintf("cell %d/%d (%s rep %s %s %s) d50_pri=%.2e equi=%s elapsed %.1f min\n",
                id, length(cells), z$scenario, z$rep, z$method, z$metric,
                rows[[rr]]$d50_pri, rows[[rr]]$equivariant,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
}

aud <- do.call(rbind, rows)
save_result(aud, "h_tie_policy_audit",
            extra = list(n_cells = nrow(aud),
                         seed_rule = "priority seed 20260719+9000+cell_index; cell_index lexicographic (scenario, rep, method, metric)",
                         material_threshold = "any displayed-value change at 3 decimals or conclusion change",
                         prereg = "H_PREREGISTRATION.md sec 4 + A1.7 (post-A2)"))
cat("max |u50 - u500| over 500-cells:",
    suppressWarnings(max(aud$d50_500, na.rm = TRUE)), "\n")
cat("max |u50 - priority|:", max(aud$d50_pri), "\n")
cat("all equivariant:", all(aud$equivariant), "\n")
cat("h_tie_policy_audit done:", format(Sys.time()), "\n")
