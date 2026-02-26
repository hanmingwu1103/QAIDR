# Run from Rcode-QAIDR directory: Rscript QAIDR/data-raw/make_data.R

# --- Cars data (27 obs x 4 interval variables) ---
cars_df <- read.csv("Cars_MM.csv")
int_min_at <- c(1, 3, 5, 7)
y_at <- 9

C <- (cars_df[, int_min_at] + cars_df[, int_min_at + 1]) / 2
R <- (cars_df[, int_min_at + 1] - cars_df[, int_min_at]) / 2
cnames <- sub("\\.[Mm]in$", "", colnames(C))
colnames(C) <- cnames
colnames(R) <- cnames

cars_mm <- structure(
  list(centers = as.matrix(C), radii = as.matrix(R),
       labels = as.factor(cars_df[, y_at])),
  class = "interval_data"
)

# --- Face data (27 obs x 6 interval variables) ---
face_df <- read.csv("facedata_MM.csv")
int_min_at_f <- seq(2, 12, 2)
y_at_f <- 1

Cf <- (face_df[, int_min_at_f] + face_df[, int_min_at_f + 1]) / 2
Rf <- (face_df[, int_min_at_f + 1] - face_df[, int_min_at_f]) / 2
cnames_f <- sub("\\.[Mm]in$", "", colnames(Cf))
colnames(Cf) <- cnames_f
colnames(Rf) <- cnames_f
rownames(Cf) <- face_df[, y_at_f]
rownames(Rf) <- face_df[, y_at_f]

facedata_mm <- structure(
  list(centers = as.matrix(Cf), radii = as.matrix(Rf),
       labels = as.factor(substr(face_df[, y_at_f], 1, 3))),
  class = "interval_data"
)

# Save to package data directory
save(cars_mm, file = "QAIDR/data/cars_mm.rda", compress = "xz")
save(facedata_mm, file = "QAIDR/data/facedata_mm.rda", compress = "xz")

cat("cars_mm:", nrow(cars_mm$centers), "obs x", ncol(cars_mm$centers), "vars\n")
cat("facedata_mm:", nrow(facedata_mm$centers), "obs x", ncol(facedata_mm$centers), "vars\n")
cat("Datasets saved successfully.\n")
