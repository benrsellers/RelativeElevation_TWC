# =============================================================================
# Relative Elevation Model (REM) Workflow
# Inputs:  GeoTIFF DSM + Shapefile centerline
# Method:  Centerline sampling -> rolling min -> spline -> nearest-point ref surface
# Optimised for: high-res DSMs (5cm), variable channel width (1-15m)
# Packages: terra, sf, zoo
# =============================================================================

library(terra)
library(sf)
library(zoo)
library(ggplot2)

# --- 0. USER INPUTS -----------------------------------------------------------

dsm_path        <- "D:/CoveringGround/TWC_RemCreation/Agisoft_proc/Wanecka/outputs/Wanecka_DSM.tif"
centerline_path <- "D:/CoveringGround/TWC_RemCreation/REMs/Wanecka/centerline.shp"
output_rem_path <- "D:/CoveringGround/TWC_RemCreation/REMs/Wanecka/Waneka_REM2.tif"

sample_spacing  <- 0.5   # Distance between centerline sample points (m)
# At 5cm DSM, 0.5m = every 10 pixels — good balance

vert_tolerance  <- 1.0   # Max elevation difference from centerline point
# allowed when filtering bad samples (m)

roll_window_m   <- 5    # Rolling minimum window length (m)
# Should be >= widest channel section (~15m here)

spline_spar     <- 0.5   # Spline smoothing (0 = wiggly, 1 = very smooth)
# Increase if profile still looks noisy after rolling min

# =============================================================================
# STEP 1: Load data & align CRS
# =============================================================================

dsm        <- rast(dsm_path)
centerline <- st_read(centerline_path) |> st_geometry()

# Reproject centerline to DSM CRS if needed
dsm_crs <- st_crs(crs(dsm))
if (!isTRUE(st_crs(centerline) == dsm_crs)) {
  message("Reprojecting centerline to match DSM CRS...")
  centerline <- st_transform(centerline, crs = dsm_crs)
}

# Merge multipart to single linestring
# centerline <- st_union(centerline) |> st_cast("LINESTRING") |> st_line_merge()
stopifnot("Centerline must be a single linestring after merging" = length(centerline) == 1)

# Extract bare sfg for geometry operations
cl_line <- centerline[[1]]

cat("DSM loaded:", nrow(dsm), "x", ncol(dsm), "cells at", res(dsm)[1], "m resolution\n")
cat("Centerline length:", round(as.numeric(st_length(centerline)), 1), "m\n")

# =============================================================================
# STEP 2: Sample elevations densely along the centerline
# =============================================================================

cat("Sampling elevations along centerline...\n")

total_len <- as.numeric(st_length(st_sfc(cl_line)))
distances <- seq(0, total_len, by = sample_spacing)

# Interpolate points along centerline
cl_pts <- lapply(distances, function(d) {
  st_line_interpolate(st_sfc(cl_line), d, normalize = FALSE)
})

# Build sf point layer
thalweg_sf <- do.call(rbind, lapply(seq_along(cl_pts), function(i) {
  st_sf(
    dist_along = distances[i],
    geometry   = st_sfc(cl_pts[[i]][[1]], crs = st_crs(centerline))
  )
}))

# Extract DSM elevation at each centerline point
thalweg_sf$elevation <- extract(dsm, vect(thalweg_sf))[, 2]

# Remove NAs (centerline extends outside DSM extent)
thalweg_sf <- thalweg_sf[!is.na(thalweg_sf$elevation), ]

cat("Centerline sample points:", nrow(thalweg_sf), "\n")

# =============================================================================
# STEP 3: Filter elevation spikes using vertical tolerance
# =============================================================================
# Water surface artifacts, boulders, and debris cause elevation spikes.
# Flag and remove points that jump more than vert_tolerance from a local median.

cat("Filtering elevation artifacts...\n")

roll_n <- max(3, round(roll_window_m / sample_spacing))  # window in # of points

# Rolling median as a local reference
local_median <- rollapply(
  thalweg_sf$elevation,
  width   = roll_n,
  FUN     = median,
  fill    = NA,
  align   = "center",
  partial = TRUE,
  na.rm   = TRUE
)

# Keep points within vert_tolerance of local median
valid     <- abs(thalweg_sf$elevation - local_median) <= vert_tolerance
n_removed <- sum(!valid)
thalweg_sf <- thalweg_sf[valid, ]

cat("Removed", n_removed, "artifact points (", round(n_removed / length(valid) * 100, 1), "%)\n")

# =============================================================================
# STEP 4: Rolling minimum to pull profile to true channel bed
# =============================================================================
# Takes the minimum within a moving window — physically this finds the lowest
# point (thalweg) within each reach window, ignoring high bank/bar points.

cat("Applying rolling minimum...\n")

thalweg_sf$elev_rollmin <- rollapply(
  thalweg_sf$elevation,
  width   = roll_n,
  FUN     = min,
  fill    = NA,
  align   = "center",
  partial = TRUE,
  na.rm   = TRUE
)

# =============================================================================
# STEP 5: Smooth the profile with a spline
# =============================================================================
# Rolling min creates a stepped profile — the spline restores a smooth gradient.

cat("Smoothing thalweg profile...\n")

spline_fit <- smooth.spline(
  x    = thalweg_sf$dist_along,
  y    = thalweg_sf$elev_rollmin,
  spar = spline_spar
)

thalweg_sf$elev_smooth <- predict(spline_fit, thalweg_sf$dist_along)$y

# --- QC plot: inspect the longitudinal profile --------------------------------
# Run this block manually to assess smoothing quality before continuing.
# Increase spline_spar if the red line is too wiggly; decrease if over-smoothed.

plot(thalweg_sf$dist_along, thalweg_sf$elevation,
     type = "l", col = "grey70", lwd = 0.8,
     xlab = "Distance downstream (m)", ylab = "Elevation (m)",
     main = "Thalweg longitudinal profile")
lines(thalweg_sf$dist_along, thalweg_sf$elev_rollmin, col = "steelblue", lwd = 1.2)
lines(thalweg_sf$dist_along, thalweg_sf$elev_smooth,  col = "red",       lwd = 2)
legend("topright",
       legend = c("Raw centerline", "Rolling min", "Smoothed spline"),
       col    = c("grey70", "steelblue", "red"),
       lty    = 1, lwd = c(1, 1.2, 2))

# =============================================================================
# STEP 6: Build reference surface — nearest centerline point method
# =============================================================================
# Strategy: rasterize each thalweg station's index onto the DSM grid, then
# use focal/distance to propagate the nearest station index to every cell.
# Avoids converting the full DSM to points.

cat("Building reference surface via nearest centerline point...\n")

# Rasterize thalweg station index (1..n) onto DSM grid
thalweg_sf$station_id <- seq_len(nrow(thalweg_sf))
thalweg_vect          <- vect(thalweg_sf)

# Rasterize station IDs — each cell gets the ID of the thalweg point that
# falls in it (NA everywhere else)
station_rast <- rasterize(thalweg_vect, dsm, field = "station_id")

# Propagate nearest station ID to all NA cells using distance-based fill
# This is equivalent to a Voronoi/nearest-neighbour assignment in raster space
station_filled <- focal(station_rast,
                        w        = 3,
                        fun      = "modal",
                        na.policy = "only",
                        na.rm    = TRUE)

# Iteratively expand until no NAs remain (handles large gaps)
max_iter <- 200
iter     <- 0
while (global(station_filled, "isNA")$isNA > 0 && iter < max_iter) {
  station_filled <- focal(station_filled,
                          w         = 3,
                          fun       = "modal",
                          na.policy = "only",
                          na.rm     = TRUE)
  iter <- iter + 1
}
cat("Station ID propagation complete after", iter, "iterations\n")

# Build lookup vector: station_id -> smoothed elevation
elev_lookup <- thalweg_sf$elev_smooth  # index 1..n

# Map station IDs to smoothed elevations using terra classify
# Build a reclass matrix: from station_id -> to elev_smooth
reclass_mat <- cbind(
  from = thalweg_sf$station_id,
  to   = thalweg_sf$elev_smooth
)
ref_surface <- classify(station_filled, reclass_mat, others = NA)

# =============================================================================
# STEP 6b: Smooth the reference surface to remove focal expansion artifacts
# =============================================================================
# The iterative focal expansion creates blocky Voronoi-like boundaries between
# station zones, which appear as sharp edges in the REM. A Gaussian smooth
# blends these boundaries into a continuous surface.

cat("Smoothing reference surface...\n")

# Gaussian kernel — sigma controls the blend radius
# sigma ~ 2x your sample_spacing is a good start (e.g. 0.5m spacing -> sigma = 1.0)
# Increase sigma if edges are still visible; don't go too large or you lose
# the downstream elevation gradient
smooth_sigma  <- sample_spacing * 2   # metres — adjust as needed
smooth_cells  <- max(3, round(smooth_sigma / res(dsm)[1]))  # convert to pixels
if (smooth_cells %% 2 == 0) smooth_cells <- smooth_cells + 1  # must be odd

gauss_kernel  <- focalMat(dsm, smooth_sigma, type = "Gauss")
ref_surface   <- focal(ref_surface,
                       w         = gauss_kernel,
                       fun       = "sum",       # pre-weighted by focalMat
                       na.policy = "omit",
                       na.rm     = TRUE)

cat("Reference surface smoothed (sigma =", smooth_sigma, "m,", smooth_cells, "x", smooth_cells, "kernel)\n")

cat("Reference surface built\n")

# =============================================================================
# STEP 7: Compute the REM
# =============================================================================

cat("Computing REM...\n")

rem <- dsm - ref_surface

# =============================================================================
# STEP 8: Save outputs
# =============================================================================

writeRaster(rem, output_rem_path, overwrite = TRUE)
cat("REM saved to:", output_rem_path, "\n")

# --- Final plot ---------------------------------------------------------------

plot(rem,
     main  = "Relative Elevation Model",
     col   = hcl.colors(100, "RdYlBu"),
     range = c(-2, 5))  # adjust range to your channel relief
plot(vect(centerline), add = TRUE, col = "black", lwd = 1.5)

cat("\nDone.\n")

# =============================================================================
#Plotting histogram 
# =============================================================================

# Compute histogram without plotting
rem_valslim <- rem[rem <= 5 & rem > 0]
h <- hist(rem_valslim, breaks = 25, plot = FALSE)

# Convert to data frame
hist_df <- data.frame(
  xmin = h$breaks[-length(h$breaks)],
  xmax = h$breaks[-1],
  count = h$counts
)

# Midpoint of each bin (used for coloring)
hist_df$mid <- (hist_df$xmin + hist_df$xmax) / 2

ggplot(hist_df) +
  geom_rect(aes(
    xmin = xmin,
    xmax = xmax,
    ymin = 0,
    ymax = count,
    fill = mid
  ),
  color = NA
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "yellow",
    high = "red",
    midpoint = median(hist_df$mid),
    name = "Relative \nElevation (m)"
  ) +
  labs(
    title = "Wanecka Grasslands \n REM Value Distribution Within Restoration Area",
    x = "Relative \nElevation (m)",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
