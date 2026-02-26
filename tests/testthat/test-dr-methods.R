test_that("umap_config_default returns a list", {
  cfg <- umap_config_default()
  expect_type(cfg, "list")
  expect_equal(cfg$n_components, 2)
  expect_equal(cfg$n_neighbors, 10)
})

test_that("umap_config_strong returns a list with extra params", {
  cfg <- umap_config_strong()
  expect_type(cfg, "list")
  expect_equal(cfg$n_epochs, 2000)
  expect_equal(cfg$negative_sample_rate, 10)
})

test_that("idr_projections print method works", {
  proj <- structure(
    list(
      "C-PCA" = list(C = matrix(0, 5, 2), R = matrix(0, 5, 2), type = "Point"),
      "V-PCA" = list(C = matrix(0, 5, 2), R = matrix(0.1, 5, 2), type = "Interval")
    ),
    class = "idr_projections"
  )
  expect_output(print(proj), "idr_projections: 2 methods")
  expect_output(print(proj), "C-PCA")
  expect_output(print(proj), "V-PCA")
})
