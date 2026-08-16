OUT <- "analysis/output"
prof <- readRDS(file.path(OUT, "h_ushcn_lambda_profile.rds"))
comp <- readRDS(file.path(OUT, "h_ushcn_lambda_components.rds"))
cat("== components ==\n")
print(comp[, c("space", "n_pairs", "n_eligible", "mean_ratio",
               "frac_dminus_zero", "mean_dminus", "mean_dplus")], digits = 3)
cat("\n== Q_TC by method x lambda (subset) ==\n")
sub <- prof[prof$lambda %in% c(0, 0.25, 0.5, 0.75, 1) | is.na(prof$lambda), ]
w <- reshape(sub[, c("Method", "lambda", "Q_TC")], idvar = "Method",
             timevar = "lambda", direction = "wide")
print(w, digits = 3, row.names = FALSE)
cat("\n== Int-UMAP full profile ==\n")
iu <- prof[prof$Method == "Int-UMAP", c("lambda", "Q_TC", "B_TC", "Q_LC",
                                        "tie_affected", "tie_spread")]
print(iu, digits = 3, row.names = FALSE)
cat("\nties by lambda:\n")
print(aggregate(tie_affected ~ lambda, prof[!is.na(prof$lambda), ], sum),
      row.names = FALSE)
