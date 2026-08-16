# Prompt H (H.8): consolidate the five Face K-profile figures into one
# faceted figure: rows = evaluation dissimilarity (4 interval + center
# baseline), columns = the six indices; lines = methods. Zero reference
# line on behavior panels (exact null mean); Q_LC panels carry the exact
# null mean K/(n-1) as a dashed curve.
SCRIPT_FILE <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]), winslash = "/", mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_FILE)
source(file.path(SCRIPT_DIR, "_common.r"))
suppressMessages(library(ggplot2))

kp <- readRDS(file.path(OUTPUT_DIR, "face_k_profiles.rds"))
pcix <- c("Q_TC", "B_TC", "Q_RE", "B_RE", "Q_LC", "B_LC")
long <- do.call(rbind, lapply(pcix, function(cl)
  data.frame(Method = kp$Method, Metric = kp$Metric, K = kp$K,
             index = cl, value = kp[[cl]], stringsAsFactors = FALSE)))
long$index <- factor(long$index, levels = pcix,
                     labels = c("Q^I[TC]", "B^I[TC]", "Q^I[RE]",
                                "B^I[RE]", "Q^I[LC]", "B^I[LC]"))
met_lv <- c("Int-Euclidean", "Hausdorff", "Ichino-Yaguchi", "Wasserstein",
            "Centers-Euclidean")
long$Metric <- factor(long$Metric, levels = met_lv,
                      labels = c("Int-Euclid.", "Hausdorff", "Ichino-Yag.",
                                 "Wasserstein", "Centers"))
meth_lv <- c("C-PCA", "V-PCA", "MR-PCA", "SPCA", "IMDS", "Int-UMAP")
long$Method <- factor(long$Method, levels = meth_lv)

n <- 27
null_df <- do.call(rbind, lapply(levels(long$Metric), function(mm) rbind(
  data.frame(Metric = mm, index = "Q^I[LC]", K = 1:25, y = (1:25) / (n - 1)),
  data.frame(Metric = mm, index = "B^I[TC]", K = 1:25, y = 0),
  data.frame(Metric = mm, index = "B^I[RE]", K = 1:25, y = 0),
  data.frame(Metric = mm, index = "B^I[LC]", K = 1:25, y = 0))))
null_df$index <- factor(null_df$index, levels = levels(long$index))
null_df$Metric <- factor(null_df$Metric, levels = levels(long$Metric))

g <- ggplot(long, aes(K, value, color = Method)) +
  geom_line(data = null_df, aes(K, y), inherit.aes = FALSE,
            linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_line(linewidth = 0.4) +
  facet_grid(Metric ~ index, scales = "free_y",
             labeller = labeller(index = label_parsed)) +
  labs(x = "Neighborhood size K", y = NULL) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(size = 7),
        legend.margin = margin(t = -5))
out <- file.path(OUTPUT_DIR, "RealData_QBK_faceted.pdf")
ggsave(out, g, width = 6.8, height = 7.6)
cat("wrote", out, "\n")
