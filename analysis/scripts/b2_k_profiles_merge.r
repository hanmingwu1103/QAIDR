# B2: merge Scenario IV K-grid rows (from b2_scenario4.r) into the K-profile
# summaries produced by b2_k_profiles.r, recomputing the aggregate and
# Kendall-stability outputs over all four scenarios. Cheap re-aggregation
# only; no embeddings are re-derived and no per-replication value changes.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))

K0 <- 10L; K_GRID <- c(5L, 10L, 20L, 40L); N_SUB <- 25L
idx <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")

long <- read.csv(file.path(OUTPUT_DIR, "b2_k_profiles_per_rep.csv"))
long <- long[long$scenario != "IV", ]        # idempotent
s4 <- file.path(OUTPUT_DIR, "scenario4_per_rep.csv")
stopifnot(file.exists(s4))
d4 <- read.csv(s4)
d4 <- d4[d4$K %in% K_GRID & d4$rep <= N_SUB, c("Method", "Metric", "K", idx, "rep")]
d4$scenario <- "IV"
long <- rbind(long[, names(d4)], d4)

summ <- aggregate(long[idx], long[c("scenario", "Method", "Metric", "K")],
                  function(v) c(mean = mean(v), sd = sd(v),
                                mcse = sd(v) / sqrt(length(v)),
                                q025 = unname(quantile(v, 0.025)),
                                q975 = unname(quantile(v, 0.975))))
taus <- list(); tt <- 0L
for (sc in unique(long$scenario)) for (met in unique(long$Metric))
  for (jj in idx) for (K in setdiff(K_GRID, K0)) for (r in seq_len(N_SUB)) {
    a <- long[long$scenario == sc & long$Metric == met & long$rep == r, ]
    v0 <- a[a$K == K0, ]; vK <- a[a$K == K, ]
    if (!nrow(v0) || !nrow(vK)) next
    mm <- intersect(v0$Method, vK$Method)
    if (length(mm) < 3) next
    tau <- suppressWarnings(cor(v0[[jj]][match(mm, v0$Method)],
                                vK[[jj]][match(mm, vK$Method)], method = "kendall"))
    tt <- tt + 1L
    taus[[tt]] <- data.frame(scenario = sc, Metric = met, index = jj, K = K,
                             rep = r, tau = tau)
  }
tau_df <- do.call(rbind, taus)
tau_summ <- aggregate(tau ~ scenario + Metric + index + K, tau_df,
                      function(v) c(mean = mean(v, na.rm = TRUE),
                                    q025 = unname(quantile(v, 0.025, na.rm = TRUE)),
                                    q975 = unname(quantile(v, 0.975, na.rm = TRUE))))
save_result(long, "b2_k_profiles_per_rep",
            inputs = s4,
            extra = list(K_grid = K_GRID, K0 = K0, n_sub = N_SUB, merged = "IV"))
save_result(do.call(data.frame, summ), "b2_k_profiles_summary")
save_result(do.call(data.frame, tau_summ), "b2_k_stability_kendall")
cat("Merged Scenario IV into K-profile summaries.\n")
