#' @importFrom rlang .data
NULL

#' Plot 2D interval projections
#'
#' Creates a faceted plot of 2D interval projections from multiple DR methods,
#' displaying intervals as rectangles.
#'
#' @param projections An \code{idr_projections} object or a named list.
#' @param labels A factor or character vector of class labels.
#' @param obs_labels Optional character vector of observation labels for text.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' data(cars_mm)
#' x <- standardize(cars_mm)
#' proj <- run_idr(x)
#' plot_projections(proj, labels = cars_mm$labels)
#' }
plot_projections <- function(projections, labels = NULL, obs_labels = NULL) {
  proj_df <- data.frame()

  for (m in names(projections)) {
    out <- projections[[m]]
    n <- nrow(out$C)

    row_labels <- if (!is.null(obs_labels)) {
      obs_labels
    } else if (!is.null(rownames(out$C))) {
      rownames(out$C)
    } else {
      seq_len(n)
    }

    proj_df <- rbind(
      proj_df,
      data.frame(
        Method = m,
        x = out$C[, 1],
        y = out$C[, 2],
        w = out$R[, 1],
        h = out$R[, 2],
        Label = row_labels,
        stringsAsFactors = FALSE
      )
    )
  }

  if (!is.null(labels)) {
    proj_df$class <- as.factor(rep(labels, length(names(projections))))
    p <- ggplot2::ggplot(proj_df) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$x - .data$w,
          xmax = .data$x + .data$w,
          ymin = .data$y - .data$h,
          ymax = .data$y + .data$h,
          color = .data$class
        ),
        fill = "blue", alpha = 0.1
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$Label),
        size = 2.5, check_overlap = TRUE
      )
  } else {
    p <- ggplot2::ggplot(proj_df) +
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = .data$x - .data$w,
          xmax = .data$x + .data$w,
          ymin = .data$y - .data$h,
          ymax = .data$y + .data$h
        ),
        fill = "blue", color = "darkblue", alpha = 0.1
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$Label),
        size = 2.5, check_overlap = TRUE
      )
  }

  p <- p +
    ggplot2::facet_wrap(~ Method, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "", x = "Dim 1", y = "Dim 2")

  p
}


#' Plot quality/behavior index profiles over K
#'
#' Creates a 2 x 3 grid of line plots (Quality and Behavior for T&C, MRRE,
#' LCMC) for a given distance metric.
#'
#' @param profile_data A data frame as returned by \code{k_profiles()}.
#' @param metric Character string specifying which distance metric to plot.
#'   If \code{NULL}, returns a list of plots for all metrics.
#' @return A \code{ggplot} (or list of plots if \code{metric} is \code{NULL}).
#' @export
#' @examples
#' \dontrun{
#' data(cars_mm)
#' x <- standardize(cars_mm)
#' proj <- run_idr(x)
#' profiles <- k_profiles(x, proj, K_max = 10)
#' plot_k_profiles(profiles, metric = "Wasserstein")
#' }
plot_k_profiles <- function(profile_data, metric = NULL) {

  my_theme_legend <- ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "inside",
      legend.direction = "horizontal",
      legend.position.inside = c(0.5, 0.2),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.key.size = ggplot2::unit(0.5, "cm"),
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 10)
    )

  my_theme_no_legend <- ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 11, face = "bold")
    )

  make_panel <- function(sub_df, met_name) {
    p1 <- ggplot2::ggplot(sub_df,
                          ggplot2::aes(x = .data$K, y = .data$Q_TC,
                                       color = .data$Method)) +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$Method), linewidth = 1) +
      my_theme_no_legend +
      ggplot2::labs(title = "Quality (T&C)", y = "Q Score") +
      ggplot2::ylim(0, 1)

    p2 <- ggplot2::ggplot(sub_df,
                          ggplot2::aes(x = .data$K, y = .data$B_TC,
                                       color = .data$Method)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$Method), linewidth = 1) +
      my_theme_no_legend +
      ggplot2::labs(title = "Behavior (T&C)", y = "B Score")

    p3 <- ggplot2::ggplot(sub_df,
                          ggplot2::aes(x = .data$K, y = .data$Q_RE,
                                       color = .data$Method)) +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$Method), linewidth = 1) +
      my_theme_no_legend +
      ggplot2::labs(title = "Quality (MRRE)", y = "Q Score") +
      ggplot2::ylim(0, 1)

    p4 <- ggplot2::ggplot(sub_df,
                          ggplot2::aes(x = .data$K, y = .data$B_RE,
                                       color = .data$Method)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$Method), linewidth = 1) +
      my_theme_no_legend +
      ggplot2::labs(title = "Behavior (MRRE)", y = "B Score")

    p5 <- ggplot2::ggplot(sub_df,
                          ggplot2::aes(x = .data$K, y = .data$Q_LC,
                                       color = .data$Method)) +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$Method), linewidth = 1) +
      my_theme_no_legend +
      ggplot2::labs(title = "Quality (LCMC)", y = "Q Score") +
      ggplot2::ylim(0, 1)

    p6 <- ggplot2::ggplot(sub_df,
                          ggplot2::aes(x = .data$K, y = .data$B_LC,
                                       color = .data$Method)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "gray") +
      ggplot2::geom_line(ggplot2::aes(linetype = .data$Method), linewidth = 1) +
      my_theme_legend +
      ggplot2::labs(title = "Behavior (LCMC)", y = "B Score")

    gridExtra::grid.arrange(
      p1, p3, p5, p2, p4, p6, ncol = 3,
      top = grid::textGrob(
        paste("Evaluation Metric:", met_name),
        gp = grid::gpar(fontsize = 16, font = 2)
      )
    )
  }

  available_metrics <- unique(profile_data$Metric)

  if (!is.null(metric)) {
    sub_df <- profile_data[profile_data$Metric == metric, ]
    make_panel(sub_df, metric)
  } else {
    plots <- list()
    for (met in available_metrics) {
      sub_df <- profile_data[profile_data$Metric == met, ]
      plots[[met]] <- make_panel(sub_df, met)
    }
    invisible(plots)
  }
}
