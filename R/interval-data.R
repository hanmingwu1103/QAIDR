#' Create an interval_data object
#'
#' Constructs an \code{interval_data} object from center and radii matrices.
#'
#' @param centers A numeric matrix of interval midpoints (n x p).
#' @param radii A numeric matrix of interval half-widths (n x p).
#' @param labels An optional factor or character vector of class labels (length n).
#' @return An object of class \code{interval_data}.
#' @export
#' @examples
#' C <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
#' R <- matrix(c(0.1, 0.2, 0.3, 0.1, 0.2, 0.3), nrow = 3)
#' x <- interval_data(C, R)
interval_data <- function(centers, radii, labels = NULL) {
  centers <- as.matrix(centers)
  radii <- as.matrix(radii)
  stopifnot(
    is.numeric(centers), is.numeric(radii),
    nrow(centers) == nrow(radii),
    ncol(centers) == ncol(radii),
    all(radii >= 0)
  )
  if (!is.null(labels)) {
    labels <- as.factor(labels)
    stopifnot(length(labels) == nrow(centers))
  }
  structure(
    list(centers = centers, radii = radii, labels = labels),
    class = "interval_data"
  )
}


#' Create interval_data from a min-max data frame
#'
#' Reads a data frame with alternating min/max columns and extracts
#' centers and radii. This corresponds to the former \code{extract_CR_from_MM}.
#'
#' @param mm_df A data frame with min/max column pairs.
#' @param int_min_at Integer vector of column indices for interval minima.
#' @param y_at Integer or \code{NULL}. Column index for class labels.
#' @return An \code{interval_data} object.
#' @export
#' @examples
#' df <- data.frame(
#'   A.min = c(1, 3), A.max = c(2, 5),
#'   B.min = c(10, 20), B.max = c(12, 25),
#'   class = c("a", "b")
#' )
#' x <- interval_data_from_mm(df, int_min_at = c(1, 3), y_at = 5)
interval_data_from_mm <- function(mm_df, int_min_at, y_at = NULL) {
  mm_df <- as.data.frame(mm_df)
  C <- (mm_df[, int_min_at, drop = FALSE] +
          mm_df[, int_min_at + 1, drop = FALSE]) / 2
  R <- (mm_df[, int_min_at + 1, drop = FALSE] -
          mm_df[, int_min_at, drop = FALSE]) / 2

  # Clean column names (strip trailing ".min" or similar suffixes)
  cnames <- colnames(C)
  cnames <- sub("\\.[Mm]in$", "", cnames)
  cnames <- sub("\\.[Ll]ower$", "", cnames)
  colnames(C) <- cnames
  colnames(R) <- cnames

  if (!is.null(rownames(mm_df))) {
    rownames(C) <- rownames(mm_df)
    rownames(R) <- rownames(mm_df)
  }

  labels <- NULL
  if (!is.null(y_at)) {
    labels <- as.factor(mm_df[, y_at])
  }

  interval_data(as.matrix(C), as.matrix(R), labels = labels)
}


#' Create interval_data from a 3D min-max array
#'
#' Converts a 3D array (n x p x 2) where \code{[,,1]} is the lower bound
#' and \code{[,,2]} is the upper bound into an \code{interval_data} object.
#'
#' @param arr A 3D numeric array of dimensions (n, p, 2).
#' @param labels An optional factor or character vector of class labels.
#' @return An \code{interval_data} object.
#' @export
#' @examples
#' arr <- array(c(1, 2, 3, 4, 5, 6, 2, 4, 5, 6, 8, 9),
#'              dim = c(3, 2, 2))
#' x <- interval_data_from_array(arr)
interval_data_from_array <- function(arr, labels = NULL) {
  stopifnot(is.array(arr), length(dim(arr)) == 3, dim(arr)[3] == 2)
  centers <- (arr[, , 1] + arr[, , 2]) / 2
  radii   <- (arr[, , 2] - arr[, , 1]) / 2
  interval_data(centers, radii, labels = labels)
}


#' Standardize interval data
#'
#' Centers the midpoints and scales both midpoints and radii by the midpoint
#' standard deviation per variable. This is the standard symbolic
#' standardization procedure.
#'
#' @param x An \code{interval_data} object.
#' @return A new \code{interval_data} object with standardized values.
#' @export
#' @examples
#' C <- matrix(rnorm(30), 10, 3)
#' R <- matrix(runif(30, 0.1, 1), 10, 3)
#' x <- interval_data(C, R)
#' xs <- standardize(x)
standardize <- function(x) {
  UseMethod("standardize")
}

#' @export
standardize.interval_data <- function(x) {
  c_mean <- colMeans(x$centers, na.rm = TRUE)
  c_sd <- apply(x$centers, 2, sd, na.rm = TRUE)
  c_sd[c_sd == 0] <- 1

  centers_new <- scale(x$centers, center = c_mean, scale = c_sd)
  radii_new <- sweep(x$radii, 2, c_sd, "/")

  # Remove scale/center attributes from scale()
  attr(centers_new, "scaled:center") <- NULL
  attr(centers_new, "scaled:scale") <- NULL

  interval_data(centers_new, radii_new, labels = x$labels)
}


#' @export
print.interval_data <- function(x, ...) {
  n <- nrow(x$centers)
  p <- ncol(x$centers)
  cat(sprintf("interval_data: %d observations, %d variables\n", n, p))
  if (!is.null(x$labels)) {
    cat(sprintf("Labels: %s\n", paste(levels(x$labels), collapse = ", ")))
  }
  invisible(x)
}


#' @export
summary.interval_data <- function(object, ...) {
  n <- nrow(object$centers)
  p <- ncol(object$centers)
  cat(sprintf("interval_data: %d observations, %d variables\n", n, p))
  cat("\nCenter summary:\n")
  print(summary(object$centers))
  cat("\nRadii summary:\n")
  print(summary(object$radii))
  if (!is.null(object$labels)) {
    cat("\nLabel distribution:\n")
    print(table(object$labels))
  }
  invisible(object)
}


#' @export
as.data.frame.interval_data <- function(x, ...) {
  p <- ncol(x$centers)
  vnames <- colnames(x$centers)
  if (is.null(vnames)) vnames <- paste0("V", seq_len(p))

  mins <- x$centers - x$radii
  maxs <- x$centers + x$radii

  df <- data.frame(matrix(0, nrow(x$centers), p * 2))
  min_at <- seq(1, p * 2, 2)
  max_at <- seq(2, p * 2, 2)
  df[, min_at] <- mins
  df[, max_at] <- maxs
  colnames(df)[min_at] <- paste0(vnames, ".min")
  colnames(df)[max_at] <- paste0(vnames, ".max")

  if (!is.null(rownames(x$centers))) {
    rownames(df) <- rownames(x$centers)
  }
  if (!is.null(x$labels)) {
    df$class <- x$labels
  }
  df
}
