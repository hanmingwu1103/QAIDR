# Internal formatting helpers (not exported)

#' Format a value with significance marker
#'
#' @param val Numeric value.
#' @param pval P-value.
#' @param alpha Significance level (default 0.05).
#' @param digits Number of decimal places (default 3).
#' @return A character string.
#' @noRd
.fmt_pval <- function(val, pval, alpha = 0.05, digits = 3) {
  s <- sprintf(paste0("%.", digits, "f"), val)
  if (pval < alpha) s <- paste0(s, "*")
  s
}
