# ==============================================================================
# Relative Elevation Model (REM) Generation Using rrrem
# ==============================================================================
#
# Description:
# This script generates a Relative Elevation Model (REM) from a digital
# elevation model (DEM) and a stream centerline using the rrrem package.
#
# A REM represents the vertical elevation difference between each raster
# cell and the nearest point along the stream centerline. REMs are commonly
# used to visualize floodplain topography, identify geomorphic surfaces,
# and support fluvial hazard and restoration analyses.
#
# Requirements:
#   - rrrem package: https://github.com/mikemahoney218/rrrem
#   - sf
#   - terra
#
# Inputs:
#   1. Stream centerline shapefile
#   2. DEM raster covering the study area
#
# Outputs:
#   1. GeoTIFF REM raster
#   2. Optional plot of the REM
#
# Author: Ben Sellers
# Date: 2026-06-15
# ==============================================================================

#------------------------------------------------------------------------------
# Install and load required packages
#------------------------------------------------------------------------------

# Install rrrem from GitHub (run once)
# install.packages("remotes")
# remotes::install_github("mikemahoney218/rrrem")

library(remotes)
library(rrrem)
library(sf)
library(terra)
library(ggplot2)

#------------------------------------------------------------------------------
# Define input and output file paths
#------------------------------------------------------------------------------

# Stream centerline shapefile
centerline_path <- "D:/CoveringGround/TWC_RemCreation/REMs/Wanecka/centerline.shp"

# Input DEM raster
dem_path <- "D:/CoveringGround/TWC_RemCreation/Agisoft_proc/Wanecka/outputs/Wanecka_DTM.tif"

# Output REM raster
output_rem_path <- "D:/CoveringGround/TWC_RemCreation/REMs/Wanecka/WaneckaREM.tif"

# Output Plots Path
output_plots_path <- "D:/CoveringGround/TWC_RemCreation/REMs/Wanecka/plots"

# Clip Results Polygon for Histogram
clip_path <- "D:/CoveringGround/TWC_RemCreation/REMs/Wanecka/cliphistorgram.shp"
#------------------------------------------------------------------------------
# Read input datasets
#------------------------------------------------------------------------------

cat("Reading centerline...\n")
centerline <- st_read(centerline_path)

cat("Reading DEM...\n")
dem <- rast(dem_path)

#------------------------------------------------------------------------------
# Ensure matching coordinate reference systems
#------------------------------------------------------------------------------

# rrrem requires the centerline and DEM to share the same CRS.
# Transform the centerline if necessary.

dem_crs <- st_crs(crs(dem))

if (st_crs(centerline) != dem_crs) {
  cat("Transforming centerline to match DEM CRS...\n")
  centerline <- st_transform(centerline, dem_crs)
}

#------------------------------------------------------------------------------
# Generate Relative Elevation Model
#------------------------------------------------------------------------------

# point number = how many points to be sampled along the centerline
#
# Example:
#   - UTM CRS → meters
#   - State Plane feet → feet

point_number <- 500

cat("Generating REM...\n")

rem <- make_rem(
  dem        = dem,
  centerline = centerline,
  n_points = point_number
)

#------------------------------------------------------------------------------
# Save REM to disk
#------------------------------------------------------------------------------

cat("Writing REM raster...\n")

writeRaster(
  rem,
  output_rem_path,
  overwrite = TRUE
)

#------------------------------------------------------------------------------
# Visualize results
#------------------------------------------------------------------------------

cat("Plotting REM...\n")

plot(
  rem,
  main = "Relative Elevation Model (REM)"
)

cat("REM creation complete.\n")
cat("Output written to:\n")
cat(output_rem_path, "\n")

# --- Final plot ---------------------------------------------------------------

plot(rem,
     main  = "Relative Elevation Model",
     col   = hcl.colors(100, "RdYlBu"),
     range = c(-2, 5))  # adjust range to your channel relief
plot(vect(centerline), add = TRUE, col = "black", lwd = 1.5)

cat("\nDone.\n")


#------------------------------------------------------------------------------
# Create a 15 m centerline buffer
#------------------------------------------------------------------------------

# The buffer distance is in the units of the DEM CRS.
# For UTM coordinate systems, this will typically be meters.

buffer_distance <- 10

cat("Creating centerline buffer...\n")

centerline_buffer <- st_buffer(
  centerline,
  dist = buffer_distance
)

cat("Clipping REM to centerline buffer...\n")

# Convert sf polygon to terra vector
buffer_vect <- vect(clip_path)

# Crop to buffer extent for efficiency
rem_crop <- crop(rem, buffer_vect)

# Mask outside the buffer
rem_clip <- mask(rem_crop, buffer_vect)

# =============================================================================
#Plotting histogram 
# =============================================================================

# Compute histogram without plotting
# Extract raster values
rem_vals <- values(rem_clip)[,1]

# Remove NAs and filter
rem_vals <- rem_vals[
  !is.na(rem_vals) &
    rem_vals > 0 &
    rem_vals <= 5
]


h <- hist(rem_vals, breaks = 25, plot = FALSE)

# Convert to data frame
hist_df <- data.frame(
  xmin = h$breaks[-length(h$breaks)],
  xmax = h$breaks[-1],
  count = h$counts
)

# Midpoint of each bin (used for coloring)
hist_df$mid <- (hist_df$xmin + hist_df$xmax) / 2

hist_df$pct <- 100 * hist_df$count / sum(hist_df$count)

ggplot(hist_df) +
  geom_rect(aes(
    xmin = xmin,
    xmax = xmax,
    ymin = 0,
    ymax = count,
    fill = mid
  ),
  color = "black"
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "gray",
    high = "white",
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
  ) +
  theme(
    panel.background = element_rect(fill = "gray90")
  )

#Each bar as a percentage of total pixels
p <- ggplot(hist_df) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = 0,
      ymax = pct,
      fill = mid
    ),
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "gray",
    high = "white",
    midpoint = median(hist_df$mid),
    name = "Relative \nElevation (m)"
  ) +
  labs(
    title = "Wanecka Grasslands\nREM Value Distribution Within Restoration Area",
    x = "Relative \nElevation (m)",
    y = "Percent of Pixels (%)"
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "gray90", color = NA),
    plot.background = element_rect(fill = "gray90", color = NA),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
plot(p)
# Write plots and final CSV
write.csv(
  hist_df,
  file.path(output_plots_path, "REM_histogram_percentages.csv"),
  row.names = FALSE
)

ggsave(
  filename = file.path(output_plots_path, "REM_histogram_percentages.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
