#' Cars interval data
#'
#' Interval-valued measurements of 27 car models across 4 variables:
#' Price, Engine Capacity (EngCap), Top Speed, and Acceleration.
#' Each observation is an interval representing the range of values
#' within a car class.
#'
#' @format An \code{interval_data} object with 27 observations and 4 variables:
#'   \describe{
#'     \item{centers}{27 x 4 matrix of interval midpoints.}
#'     \item{radii}{27 x 4 matrix of interval half-widths.}
#'     \item{labels}{Factor with 4 levels: Berlina, Luxury, Sportive, Utilitarian.}
#'   }
#' @source Billard, L. and Diday, E. (2006).
#'   \emph{Symbolic Data Analysis: Conceptual Statistics and Data Mining}.
#'   Wiley.
"cars_mm"


#' Face anthropometry interval data
#'
#' Interval-valued anthropometric face measurements of 27 individuals
#' across 6 variables: AD, BC, AH, DH, EH, GH (distances between
#' facial landmarks).
#'
#' @format An \code{interval_data} object with 27 observations and 6 variables:
#'   \describe{
#'     \item{centers}{27 x 6 matrix of interval midpoints.}
#'     \item{radii}{27 x 6 matrix of interval half-widths.}
#'     \item{labels}{Factor with 9 levels (3-letter person identifiers).}
#'   }
#' @source Billard, L. and Diday, E. (2006).
#'   \emph{Symbolic Data Analysis: Conceptual Statistics and Data Mining}.
#'   Wiley.
"facedata_mm"
