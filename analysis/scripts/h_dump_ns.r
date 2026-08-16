# Dump deparsed bodies of every exported+internal QAIDR function for a
# functional diff between package builds. Usage:
#   Rscript h_dump_ns.r <outfile> [libpath]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) .libPaths(c(args[2], .libPaths()))
suppressMessages(library(QAIDR))
ns <- asNamespace("QAIDR")
objs <- sort(ls(ns, all.names = TRUE))
con <- file(args[1], "w", encoding = "UTF-8")
cat("QAIDR version:", as.character(utils::packageVersion("QAIDR")), "\n",
    "path:", find.package("QAIDR"), "\n\n", file = con)
for (o in objs) {
  v <- get(o, envir = ns)
  if (is.function(v)) {
    cat("#### ", o, " ####\n", file = con)
    writeLines(deparse(v, width.cutoff = 120), con)
    cat("\n", file = con)
  }
}
close(con)
cat("dumped", length(objs), "objects\n")
