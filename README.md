[README.md](https://github.com/user-attachments/files/28888845/README.md)

![REM header](waneckagithubheader.png)


# Relative Elevation Model (REM) Generation

An R script that generates a Relative Elevation Model (REM) from a digital elevation model (DEM) and a stream centerline, using the [`rrrem`](https://github.com/mikemahoney218/rrrem) package. The script also clips the REM to a buffer/polygon of interest and produces summary histograms of relative elevation values.

## What is a REM?

A Relative Elevation Model represents the vertical elevation difference between each raster cell and the nearest point along a stream centerline. REMs are widely used in fluvial geomorphology to visualize floodplain topography, delineate geomorphic surfaces (e.g., active channel, floodplain, terraces), and support river restoration and flood hazard analyses.

## Features

- Reads a stream centerline shapefile and a DEM raster
- Automatically reprojects the centerline to match the DEM's CRS if needed
- Generates a REM using `rrrem::make_rem()`
- Exports the REM as a GeoTIFF
- Plots the REM with the centerline overlaid
- Clips the REM to a buffer/polygon for focused analysis
- Computes and visualizes histograms of relative elevation values (raw counts and percent of pixels)
- Exports histogram data as CSV and the plot as a PNG

## Requirements

- R (≥ 4.0 recommended)
- Packages:
  - [`rrrem`](https://github.com/mikemahoney218/rrrem) (installed from GitHub via `remotes`)
  - `sf`
  - `terra`
  - `ggplot2`
  - `remotes`

Install dependencies:

```r
install.packages(c("sf", "terra", "ggplot2", "remotes"))
remotes::install_github("mikemahoney218/rrrem")
```

## Inputs

| Input | Description |
|---|---|
| Stream centerline shapefile | Vector line representing the stream/channel centerline |
| DEM raster (GeoTIFF) | Digital elevation model covering the study area |
| Clip polygon shapefile | Polygon used to clip the REM for the histogram analysis |

## Outputs

| Output | Description |
|---|---|
| `*REM.tif` | REM raster (GeoTIFF), full extent |
| REM plot | On-screen plot of the REM with centerline overlay |
| `REM_histogram_percentages.csv` | Histogram bin data (counts and percentages) |
| `REM_histogram_percentages.png` | Histogram plot of relative elevation distribution |

## Usage

1. Open the script in R or RStudio.
2. Update the file paths in the **Define input and output file paths** section:
   ```r
   centerline_path     <- "path/to/centerline.shp"
   dem_path             <- "path/to/dem.tif"
   output_rem_path      <- "path/to/output_REM.tif"
   output_plots_path    <- "path/to/plots_folder"
   clip_path            <- "path/to/clip_polygon.shp"
   ```
3. Adjust parameters as needed:
   - `point_number` — number of points sampled along the centerline (default: 500)
   - `buffer_distance` — buffer distance around the centerline, in DEM CRS units (default: 10)
   - Histogram value range and bin count (`breaks`, value filtering bounds)
4. Run the script. It will:
   - Generate and save the REM raster
   - Display REM plots
   - Clip the REM to the buffer/polygon
   - Generate and save histogram outputs

## Notes

- The centerline and DEM must share the same coordinate reference system (CRS); the script handles reprojection automatically if they differ.
- Buffer and point-sampling distances are in the units of the DEM's CRS (e.g., meters for UTM, feet for State Plane in feet).
- The plotting color range (`range = c(-2, 5)`) should be adjusted to match the relief of your study site.

## Author

Ben Sellers

## License

Add a license of your choice (e.g., MIT) before publishing.
