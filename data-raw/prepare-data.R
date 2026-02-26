## Prepare package datasets from CSV files

# --- Cars dataset ---
cars_df <- read.csv("Cars_MM.csv")
cars_mm <- QAIDR::interval_data_from_mm(
  cars_df,
  int_min_at = c(1, 3, 5, 7),
  y_at = 9
)
usethis::use_data(cars_mm, overwrite = TRUE)


# --- Face dataset ---
face_df <- read.csv("facedata_MM.csv")
facedata_mm <- QAIDR::interval_data_from_mm(
  face_df,
  int_min_at = seq(2, 12, 2),
  y_at = 1
)
# Shorten labels to 3-character prefix for class grouping
facedata_mm$labels <- as.factor(substr(as.character(facedata_mm$labels), 1, 3))
usethis::use_data(facedata_mm, overwrite = TRUE)
