test_that("plot_projections returns a ggplot", {
  proj <- structure(
    list(
      "A" = list(
        C = matrix(rnorm(10), 5, 2),
        R = matrix(0.1, 5, 2),
        type = "Interval"
      )
    ),
    class = "idr_projections"
  )

  p <- plot_projections(proj, labels = c("a", "a", "b", "b", "a"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_projections works without labels", {
  proj <- structure(
    list(
      "A" = list(
        C = matrix(rnorm(10), 5, 2),
        R = matrix(0.1, 5, 2),
        type = "Interval"
      )
    ),
    class = "idr_projections"
  )

  p <- plot_projections(proj)
  expect_s3_class(p, "ggplot")
})
