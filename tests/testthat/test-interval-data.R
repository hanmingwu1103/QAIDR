test_that("interval_data constructor works", {
  C <- matrix(1:6, nrow = 3)
  R <- matrix(rep(0.1, 6), nrow = 3)
  x <- interval_data(C, R)

  expect_s3_class(x, "interval_data")
  expect_equal(nrow(x$centers), 3)
  expect_equal(ncol(x$centers), 2)
  expect_null(x$labels)
})

test_that("interval_data rejects mismatched dimensions", {
  C <- matrix(1:6, nrow = 3)
  R <- matrix(1:4, nrow = 2)
  expect_error(interval_data(C, R))
})

test_that("interval_data rejects negative radii", {
  C <- matrix(1:6, nrow = 3)
  R <- matrix(c(0.1, -0.1, 0.1, 0.1, 0.1, 0.1), nrow = 3)
  expect_error(interval_data(C, R))
})

test_that("interval_data with labels works", {
  C <- matrix(1:6, nrow = 3)
  R <- matrix(rep(0.1, 6), nrow = 3)
  x <- interval_data(C, R, labels = c("a", "b", "a"))

  expect_s3_class(x$labels, "factor")
  expect_equal(length(x$labels), 3)
})

test_that("interval_data_from_mm works", {
  df <- data.frame(
    A.min = c(1, 3), A.max = c(2, 5),
    B.min = c(10, 20), B.max = c(12, 25),
    class = c("a", "b")
  )
  x <- interval_data_from_mm(df, int_min_at = c(1, 3), y_at = 5)

  expect_s3_class(x, "interval_data")
  expect_equal(nrow(x$centers), 2)
  # Centers should be (min+max)/2
  expect_equal(x$centers[1, 1], 1.5)
  expect_equal(x$centers[1, 2], 11)
  # Radii should be (max-min)/2
  expect_equal(x$radii[1, 1], 0.5)
  expect_equal(x$radii[1, 2], 1)
  expect_equal(as.character(x$labels), c("a", "b"))
})

test_that("interval_data_from_array works", {
  arr <- array(0, dim = c(3, 2, 2))
  arr[, , 1] <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
  arr[, , 2] <- matrix(c(2, 4, 5, 12, 25, 35), nrow = 3)
  x <- interval_data_from_array(arr)

  expect_s3_class(x, "interval_data")
  expect_equal(x$centers[1, 1], 1.5)
  expect_equal(x$radii[1, 1], 0.5)
})

test_that("interval_data_from_dataSDA parses lower-upper character columns", {
  d <- data.frame(
    A = c("1,2", "3,5"),
    B = c("10.0,12.0", "20.0,25.0"),
    group = c("g1", "g2"),
    row.names = c("one", "two")
  )
  x <- interval_data_from_dataSDA(d, label_col = "group")

  expect_s3_class(x, "interval_data")
  expect_equal(x$centers, structure(matrix(c(1.5, 4, 11, 22.5), 2, 2),
                                    dimnames = list(c("one", "two"),
                                                    c("A", "B"))))
  expect_equal(x$radii, structure(matrix(c(0.5, 1, 1, 2.5), 2, 2),
                                  dimnames = list(c("one", "two"),
                                                  c("A", "B"))))
  expect_equal(as.character(x$labels), c("g1", "g2"))
})

test_that("interval_data_from_dataSDA converts the CRAN example data", {
  skip_if_not_installed("dataSDA")

  utils::data("cars.int", package = "dataSDA", envir = environment())
  cars <- interval_data_from_dataSDA(cars.int)
  expect_equal(dim(cars$centers), c(27L, 4L))
  expect_equal(rownames(cars$centers)[1], "Alfa145")
  expect_equal(as.character(cars$labels)[1:3],
               c("Utilitarian", "Berlina", "Sportive"))
  expect_equal(cars$centers[1, "Acceleration"], 9.75)

  utils::data("face.iGAP", package = "dataSDA", envir = environment())
  person <- sub("[[:digit:]]+$", "", rownames(face.iGAP))
  face <- interval_data_from_dataSDA(face.iGAP, labels = person)
  expect_equal(dim(face$centers), c(27L, 6L))
  expect_equal(face$centers["FRA1", "AD"], 156)
  expect_equal(face$radii["FRA1", "BC"], 1.505)
  expect_equal(levels(face$labels),
               c("FRA", "HUS", "INC", "ISA", "JPL", "KHA", "LOT",
                 "PHI", "ROM"))
})

test_that("interval_data_from_dataSDA validates ambiguous inputs", {
  d <- data.frame(A = c("1,2", "3,4"), group = c("a", "b"))
  expect_error(interval_data_from_dataSDA(d, label_col = "group",
                                          labels = c("a", "b")),
               "either")
  expect_error(interval_data_from_dataSDA(data.frame(x = 1:2)),
               "No dataSDA interval")
  expect_error(interval_data_from_dataSDA(
    data.frame(A = c("2,1", "3,4"))), "lower <= upper", fixed = TRUE)
})

test_that("standardize works correctly", {
  set.seed(42)
  C <- matrix(c(10, 20, 30, 100, 200, 300), nrow = 3)
  R <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
  x <- interval_data(C, R)
  xs <- standardize(x)

  expect_s3_class(xs, "interval_data")
  # After standardization, column means of centers should be ~0
  expect_equal(colMeans(xs$centers), c(0, 0), tolerance = 1e-10)
  # Column SDs of centers should be ~1
  expect_equal(apply(xs$centers, 2, sd), c(1, 1), tolerance = 1e-10)
})

test_that("as.data.frame.interval_data roundtrips", {
  C <- matrix(c(1.5, 4, 11, 22.5), nrow = 2)
  R <- matrix(c(0.5, 1, 1, 2.5), nrow = 2)
  x <- interval_data(C, R, labels = c("a", "b"))
  df <- as.data.frame(x)

  expect_equal(ncol(df), 5)  # 4 min/max cols + 1 class
  expect_true("class" %in% names(df))
  # min = center - radius
  expect_equal(df[1, 1], 1.0)
  # max = center + radius
  expect_equal(df[1, 2], 2.0)
})

test_that("print.interval_data works", {
  C <- matrix(1:6, nrow = 3)
  R <- matrix(rep(0.1, 6), nrow = 3)
  x <- interval_data(C, R)
  expect_output(print(x), "interval_data: 3 observations, 2 variables")
})
