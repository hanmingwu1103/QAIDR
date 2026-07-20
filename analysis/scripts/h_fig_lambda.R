# Prompt H: compact USHCN lambda-profile figure (Q_TC vs lambda, six
# methods; center-only baseline as dashed reference lines).
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.R"))
suppressMessages(library(ggplot2))

prof <- readRDS(file.path(OUTPUT_DIR, "h_ushcn_lambda_profile.rds"))
d <- prof[!is.na(prof$lambda), ]
base <- prof[is.na(prof$lambda), c("Method", "Q_TC")]
meth_lv <- c("C-PCA", "V-PCA", "MR-PCA", "SPCA", "IMDS", "Int-UMAP")
d$Method <- factor(d$Method, levels = meth_lv)
base$Method <- factor(base$Method, levels = meth_lv)

g <- ggplot(d, aes(lambda, Q_TC, color = Method, shape = Method)) +
  geom_hline(data = base, aes(yintercept = Q_TC, color = Method),
             linetype = "dashed", linewidth = 0.3, show.legend = FALSE) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.4) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = expression(lambda), y = expression(Q[TC]^I)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.margin = margin(t = -6))
out <- file.path(OUTPUT_DIR, "USHCN_lambda_profile.pdf")
ggsave(out, g, width = 6.2, height = 3.4)
cat("wrote", out, "\n")
