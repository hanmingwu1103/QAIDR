# Internal conversion helpers (not exported)

#' Convert interval_data to a 3D symbolicDA array
#'
#' @param x An \code{interval_data} object.
#' @return A 3D array (n x p x 2) with dimnames.
#' @noRd
.to_sym_array <- function(x) {
  if (inherits(x, "interval_data")) {
    centers <- x$centers
    radii <- x$radii
  } else {
    # Accept raw centers/radii
    centers <- x$centers
    radii <- x$radii
  }
  N <- nrow(centers)
  P <- ncol(centers)
  arr <- array(0, dim = c(N, P, 2))
  arr[, , 1] <- centers - radii
  arr[, , 2] <- centers + radii

  rnames <- rownames(centers)
  if (is.null(rnames)) rnames <- paste0("O", seq_len(N))
  cnames <- colnames(centers)
  if (is.null(cnames)) cnames <- paste0("V", seq_len(P))

  dimnames(arr) <- list(rnames, cnames, c("Min", "Max"))
  arr
}


#' Convert a 3D array to a min-max data frame
#'
#' @param data_MM_array A 3D array (n x p x 2).
#' @return A data frame with alternating .min/.max columns.
#' @noRd
.to_mm_df <- function(data_MM_array) {
  n <- dim(data_MM_array)[1]
  p <- dim(data_MM_array)[2]

  MM_df <- matrix(0, ncol = p * 2, nrow = n)
  min_at <- seq(1, p * 2, 2)
  max_at <- seq(1, p * 2, 2) + 1
  MM_df[, min_at] <- as.matrix(data_MM_array[, , 1])
  MM_df[, max_at] <- as.matrix(data_MM_array[, , 2])
  MM_df <- as.data.frame(MM_df)

  vnames <- dimnames(data_MM_array)[[2]]
  if (!is.null(vnames)) {
    colnames(MM_df)[min_at] <- paste0(vnames, ".min")
    colnames(MM_df)[max_at] <- paste0(vnames, ".max")
  }
  MM_df
}


#' Convert min-max data frame to RSDA symbolic object
#'
#' @param df A data frame with .min/.max columns.
#' @param class_name Optional class column name.
#' @return An RSDA symbolic data object.
#' @noRd
.to_rsda <- function(df, class_name = NULL) {
  if (!requireNamespace("RSDA", quietly = TRUE)) {
    stop("Package 'RSDA' is required. Install it with install.packages('RSDA').")
  }

  min_cols <- grep("\\.min$", names(df), value = TRUE)
  base_vars <- sub("\\.min$", "", min_cols)

  min_at <- grep("\\.min", names(df))
  max_at <- grep("\\.max", names(df))

  df_min <- df[, min_at, drop = FALSE]
  df_max <- df[, max_at, drop = FALSE]
  colnames(df_min) <- base_vars
  colnames(df_max) <- base_vars

  df_tmp <- rbind(df_min, df_max)

  if (!is.null(class_name)) {
    df_long <- data.frame(temp_ID = seq_len(nrow(df)), df_tmp, df[class_name])
    colnames(df_long)[ncol(df_long)] <- class_name
    at <- which(vapply(df_long, class, character(1)) == "character")
    if (length(at) != 0) {
      df_long[, at] <- lapply(df_long[, at, drop = FALSE], as.factor)
    }
  } else {
    df_long <- data.frame(temp_ID = seq_len(nrow(df)), df_tmp)
  }

  rsda_obj <- RSDA::classic.to.sym(
    df_long,
    concept = "temp_ID",
    default.categorical = RSDA::sym.set
  )
  rsda_obj
}
