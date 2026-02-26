#' Generate interval Swiss roll data
#'
#' Generates a synthetic interval-valued dataset embedded on a Swiss roll
#' manifold. Useful for benchmarking interval dimensionality reduction methods.
#'
#' @param n Integer number of observations (default 800).
#' @param seed Integer random seed (default 1).
#' @param t_min Minimum angle parameter (default \code{1.5 * pi}).
#' @param t_max Maximum angle parameter (default \code{4.5 * pi}).
#' @param y_max Maximum value for the linear dimension (default 70).
#' @param noise_sd Standard deviation of Gaussian noise on centers (default 0.02).
#' @param r_min Minimum interval radius (default 0.2).
#' @param r_max Maximum interval radius (default 0.5).
#' @return An \code{interval_data} object with 3 variables (x, y, z).
#' @export
#' @examples
#' dat <- gen_interval_swissroll(n = 100, seed = 42)
#' print(dat)
gen_interval_swissroll <- function(n = 800, seed = 1,
                                   t_min = 1.5 * pi, t_max = 4.5 * pi,
                                   y_max = 70,
                                   noise_sd = 0.02,
                                   r_min = 0.2, r_max = 0.5) {
  set.seed(seed)
  t <- sort(runif(n, t_min, t_max))
  x <- t * cos(t)
  z <- t * sin(t)
  y <- runif(n, 0, y_max)

  centers <- cbind(x, y, z) + matrix(rnorm(n * 3, sd = noise_sd), n, 3)
  radii <- matrix(runif(n * 3, r_min, r_max), n, 3)
  colnames(centers) <- c("x", "y", "z")
  colnames(radii) <- c("x", "y", "z")

  interval_data(centers, radii)
}
