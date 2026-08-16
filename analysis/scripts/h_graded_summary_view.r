OUT <- "analysis/output"
s <- readRDS(file.path(OUT, "h_graded_width_summary.rds"))
p <- readRDS(file.path(OUT, "h_graded_width_per_rep.rds"))
cat("== primary estimand: paired Qc_TC(Wasserstein) - Qc_TC(center) ==\n")
pri <- s[s$primary, c("level", "Method", "mean_qc_diff", "mcse")]
pri <- pri[order(pri$Method, pri$level), ]
print(pri, digits = 3, row.names = FALSE)
cat("\n== secondary: CR-Euclidean Q_TC contrast ==\n")
cr <- s[s$Metric == "CR-Euclidean" & s$index == "Q_TC",
        c("level", "Method", "mean_qc_diff", "mcse")]
cr <- cr[order(cr$Method, cr$level), ]
print(cr, digits = 3, row.names = FALSE)
cat("\n== raw Q_TC means (Wasserstein vs Centers) by level ==\n")
a <- aggregate(Q_TC ~ level + Metric, p[p$Metric %in% c("Wasserstein",
      "Centers-Euclidean", "CR-Euclidean"), ], mean)
print(reshape(a, idvar = "Metric", timevar = "level", direction = "wide"),
      digits = 3, row.names = FALSE)
cat("\nmeta: l3 agreement max delta recorded in per_rep meta json\n")
cat("tie-affected cells by level:\n")
print(aggregate(tie_affected ~ level, p, sum), row.names = FALSE)
