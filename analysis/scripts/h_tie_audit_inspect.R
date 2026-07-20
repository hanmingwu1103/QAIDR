a <- readRDS("analysis/output/h_tie_policy_audit.rds")
cat("cells:", nrow(a), "\n")
cat("equivariant FALSE:", sum(!a$equivariant), "\n")
print(a[!a$equivariant, c("scenario", "rep", "Method", "Metric",
                          "cell_index", "d50_pri", "equivariant")],
      row.names = FALSE)
cat("\nd50_pri summary by scenario:\n")
print(aggregate(d50_pri ~ scenario, a, function(v)
  c(max = max(v), med = stats::median(v))))
cat("\nd50_500 max (where run):", suppressWarnings(max(a$d50_500, na.rm = TRUE)), "\n")
cat("cells with d50_pri > 1e-3:", sum(a$d50_pri > 1e-3), "\n")
print(a[a$d50_pri > 1e-3, c("scenario", "rep", "Method", "Metric", "d50_pri")],
      row.names = FALSE)
