# B2 Part 2.7: Scenario II diagnostic warning-capture run (quick mode; NO
# manuscript numbers are produced). Runs 2 replications of the exact Scenario
# II code path with withCallingHandlers, records every warning with its stage
# and call, and writes a classified inventory.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
source(file.path(ANALYSIS_DIR, "scripts", "scenarios.r"))
save_sessioninfo("scenario2_warning_diagnostic")

warns <- list(); ww <- 0L
grab <- function(stage, expr) {
  withCallingHandlers(expr, warning = function(w) {
    ww <<- ww + 1L
    warns[[ww]] <<- data.frame(stage = stage, message = conditionMessage(w),
                               call = paste(deparse(conditionCall(w))[1], collapse = ""))
    invokeRestart("muffleWarning")
  })
}

for (r in 1:2) {
  seed_r <- SEED_BASE + 2000L + r
  x_raw <- grab("gen", gen_scenario2(seed_r))
  x_std <- grab("standardize", standardize(x_raw))
  set.seed(seed_r + 500000L)
  rm_raw <- grab("run_idr", run_methods_timed(x_raw))
  invisible(grab("evaluate", evaluate_cells(x_std, rm_raw$projections, K = 10,
                                            tie_seed = seed_r)))
  invisible(grab("assess_perm", assess_quality(x_std, rm_raw$projections, K = 10,
                                               perm_test = TRUE, n_perm = 19,
                                               seed = seed_r + 900000L,
                                               ties = "random", baseline = TRUE)))
}

inv <- if (ww > 0) do.call(rbind, warns) else
  data.frame(stage = character(), message = character(), call = character())
inv$class <- NA_character_
cls <- function(pat, label) inv$class[grepl(pat, inv$message)] <<- label
cls("failed creating initial embedding", "benign, documented: umap spectral initialization failed (disconnected/low-connectivity kNN graph at n=800); the umap package falls back to a random initial embedding; run is seeded so results are reproducible; affects Int-UMAP only")
cls("ties[.]method|rank", "expected: seeded random tie-breaking in rank()")
cls("no non-missing arguments|NaN|collapsing to unique", "numeric edge in dependency")
cls("umap|UMAP|n_neighbors|epochs", "benign: umap/RSDA configuration notice")
cls("recycl", "vector recycling in dependency (symbolicDA/RSDA internals)")
cls("deprecat", "benign: deprecation notice from dependency")
cls("index .* outside|out of range", "diagnostic: index range check")
inv$class[is.na(inv$class)] <- "UNCLASSIFIED - REVIEW"
tab <- aggregate(cbind(count = rep(1, nrow(inv))) ~ stage + class + message, inv, sum)
save_result(tab, "scenario2_warning_inventory",
            extra = list(n_reps = 2, n_warnings = ww,
                         unclassified = sum(inv$class == "UNCLASSIFIED - REVIEW")))
print(tab[, c("count", "stage", "class")])
cat("\nTotal warnings:", ww, " unclassified:",
    sum(inv$class == "UNCLASSIFIED - REVIEW"), "\n")
if (sum(inv$class == "UNCLASSIFIED - REVIEW") > 0) {
  cat("UNCLASSIFIED messages:\n")
  print(unique(inv$message[inv$class == "UNCLASSIFIED - REVIEW"]))
}
