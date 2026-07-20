OUT <- "analysis/output"
p1 <- readRDS(file.path(OUT, "scenario1_per_rep.rds"))
cat("S1 flagged rows:", sum(p1$tie_affected), "of", nrow(p1),
    "reps:", paste(sort(unique(p1$rep[p1$tie_affected])), collapse = ","), "\n")
p2 <- readRDS(file.path(OUT, "scenario2_per_rep.rds"))
cat("S2 cols:", paste(names(p2), collapse = " "), "\n")
if ("tie_affected" %in% names(p2))
  cat("S2 flagged rows:", sum(p2$tie_affected),
      "reps:", paste(sort(unique(p2$rep[p2$tie_affected])), collapse = ","), "\n")
p4 <- readRDS(file.path(OUT, "scenario4_per_rep.rds"))
cat("S4 cols:", paste(names(p4), collapse = " "), "\n")
sp <- grep("_spread$", names(p4), value = TRUE)
if (length(sp)) {
  m <- as.matrix(p4[p4$K == 10 & p4$rep <= 25, sp])
  keep <- p4$K == 10 & p4$rep <= 25
  aff <- rowSums(m > 0) > 0
  cat("S4 K=10 rep<=25 cells with spread>0:", sum(aff), "of", sum(keep),
      "reps:", paste(sort(unique(p4$rep[keep][aff])), collapse = ","), "\n")
}
mi <- readRDS(file.path(OUT, "midscale_indices.rds"))
cat("USHCN daggered cells:", sum(mi$tie_affected), "of", nrow(mi), "\n")
print(mi[mi$tie_affected, c("Method", "Metric")], row.names = FALSE)
