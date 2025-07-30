#' Suggest Optimal Number of Principal Components
#'
#' This function analyzes a PCA object to determine the recommended number of 
#' principal components (PCs) to retain, using three common criteria:
#' - Cumulative explained variance (threshold-based)
#' - Kaiser criterion (eigenvalue > 1)
#' - Elbow method (first minimum of successive variance drops)
#' 
#' Optionally, a diagnostic plot (scree and cumulative variance) is shown.
#'
#' @param pca_obj A PCA object returned by [prcomp()] or similar. Must contain `$sdev`.
#' @param variance_cutoff Numeric; cumulative variance threshold for selecting PCs 
#'                        (default is 0.8 = 80% explained variance).
#' @param plot Logical; if `TRUE`, a scree plot and cumulative variance plot are displayed.
#'
#' @return A list containing:
#' \describe{
#'   \item{n_pcs}{Number of PCs recommended based on the variance threshold.}
#'   \item{explained_variance}{Number of PCs needed to reach the variance cutoff.}
#'   \item{kaiser}{Number of PCs with eigenvalue > 1 (Kaiser criterion).}
#'   \item{elbow}{Position of elbow point (first minimal drop in explained variance).}
#'   \item{info_table}{A summary table showing the values for each criterion.}
#' }
#'
#' @details
#' - **Cumulative Variance Threshold**: Retain the smallest number of components such that
#'   the cumulative proportion of explained variance meets or exceeds `variance_cutoff`.
#'
#' - **Kaiser Criterion**: Retain all components with eigenvalues greater than 1. Assumes
#'   data has been standardized (mean-centered and scaled). See Kaiser (1960).
#'
#' - **Elbow Method**: Finds the index where the decrease in explained variance is smallest,
#'   i.e., where the "knee" or "elbow" appears in the scree plot.
#'
#' @references
#' - Jolliffe, I. T. (2002). *Principal Component Analysis*. Springer Series in Statistics.
#' - Kaiser, H. F. (1960). The application of electronic computers to factor analysis.
#'   *Educational and Psychological Measurement*, 20(1), 141–151.
#' - Cattell, R. B. (1966). The scree test for the number of factors. *Multivariate Behavioral Research*, 1(2), 245–276.
#'
#' @examples
#' pca <- prcomp(USArrests, scale. = TRUE)
#' suggest_n_pcs(pca, variance_cutoff = 0.9)
#'
#' @export
suggest_n_pcs <- function(pca_obj, variance_cutoff = 0.8, plot = TRUE) {
  # Extract standard deviations of the principal components
  std_dev <- pca_obj$sdev
  
  # Compute proportion of variance explained by each PC
  var_explained <- std_dev^2 / sum(std_dev^2)
  
  # Compute cumulative explained variance
  cum_var_explained <- cumsum(var_explained)
  
  # Criterion 1: Number of components needed to reach variance_cutoff
  n_var <- which(cum_var_explained >= variance_cutoff)[1]
  
  # Criterion 2: Kaiser criterion – eigenvalue > 1
  eigenvalues <- std_dev^2
  n_kaiser <- sum(eigenvalues > 1)
  
  # Criterion 3: Elbow method – where decrease in explained variance flattens
  diffs <- diff(var_explained)
  elbow <- which.min(diffs)[1]
  
  # Assemble criterion comparison table
  info_table <- data.frame(
    Criterion = c("Variance Cutoff", "Kaiser Criterion", "Elbow Method"),
    Num_Components = c(n_var, n_kaiser, elbow),
    Cumulative_Explained_Variance = c(
      round(cum_var_explained[n_var], 3),
      round(cum_var_explained[n_kaiser], 3),
      round(cum_var_explained[elbow], 3)
    )
  )
  
  # Final recommendation is based on variance_cutoff only (can be modified as needed)
  n_final <- n_var
  
  # Print summary
  cat("📊 Summary of PCA Component Selection Criteria:\n")
  print(info_table, row.names = FALSE)
  cat("\n✅ Recommended number of PCs (based on variance_cutoff =", variance_cutoff, "):", n_final, "\n")
  
  # Optional plots
  if (plot) {
    par(mfrow = c(1, 2))
    
    # Scree plot: individual variance explained
    plot(var_explained, type = "b", pch = 19, col = "steelblue",
         xlab = "Component", ylab = "Explained Variance",
         main = "Scree Plot")
    abline(h = 1, col = "red", lty = 2)     # Kaiser line
    abline(v = elbow, col = "darkgreen", lty = 3)  # Elbow marker
    legend("topright", legend = c("Kaiser (λ > 1)", "Elbow"),
           col = c("red", "darkgreen"), lty = c(2, 3), bty = "n")
    
    # Cumulative variance plot
    plot(cum_var_explained, type = "b", pch = 19, col = "darkorange",
         xlab = "Component", ylab = "Cumulative Variance",
         main = "Cumulative Explained Variance")
    abline(h = variance_cutoff, col = "red", lty = 2)
    abline(v = n_var, col = "blue", lty = 3)
    legend("bottomright", legend = c("Cutoff", "Selected Components"),
           col = c("red", "blue"), lty = c(2, 3), bty = "n")
    
    par(mfrow = c(1, 1))  # reset plotting layout
  }
  
  # Return results silently for use in pipelines
  invisible(list(
    n_pcs = n_final,
    explained_variance = n_var,
    kaiser = n_kaiser,
    elbow = elbow,
    info_table = info_table
  ))
}



#' Split Z Coordinates into Vertical Slices and Count Points per Slice
#'
#' This function takes a vector of Z-coordinates (heights) and bins them into
#' 1-meter horizontal slices. It returns the count of points in each slice, ensuring that
#' all slices from 0 to `maxZ` are represented, even if some slices have zero points.
#'
#' @param Z A numeric vector of Z coordinates (e.g., heights of LiDAR points in meters).
#' @param maxZ Integer; the maximum height to consider (defines the highest slice boundary).
#'
#' @return A named list containing point counts per 1-meter height slice. The names are
#' formatted for clarity (e.g., `"ground_0_1m"`, `"pulses_1_2m"`, …).
#'
#' @details
#' - This is a foundational step in computing vertical vegetation structure such as
#'   Leaf Area Density (LAD) profiles.
#' - The slicing assumes a 1-meter vertical resolution and bins by floor(Z).
#' - Empty slices (no points) are included with count 0 to preserve structure for later matrix assembly.
#'
#' @examples
#' z_vals <- runif(1000, 0, 20)
#' pointsByZSlice(z_vals, maxZ = 20)
#'
#' @export
pointsByZSlice <- function(Z, maxZ) {
  # Floor Z-values to get integer bin index (0-based)
  heightSlices <- as.integer(Z)
  
  # Create data.table for potential grouping (not used further here)
  zSlice <- data.table::data.table(Z = Z, heightSlices = heightSlices)
  
  # Count number of points per height slice using base aggregate
  sliceCount <- stats::aggregate(list(V1 = Z), list(heightSlices = heightSlices), length)
  sliceCount$V1 <- as.numeric(sliceCount$V1)  # Ensure numeric (not integer or factor)
  
  # Ensure all expected slice bins [0, maxZ] exist (fill with 0 if missing)
  colRange <- 0:maxZ
  missing <- colRange[!colRange %in% sliceCount$heightSlices]
  if (length(missing) > 0) {
    fill <- data.frame(heightSlices = missing, V1 = as.numeric(0))
    sliceCount <- rbind(sliceCount, fill)
  }
  
  # Order slices from bottom to top
  sliceCount <- sliceCount[order(sliceCount$heightSlices), ]
  
  # Create readable column names for each slice
  colNames <- as.character(sliceCount$heightSlices)
  colNames[1] <- "ground_0_1m"  # Name for the lowest bin
  colNames[-1] <- paste0("pulses_", sliceCount$heightSlices[-1], "_", sliceCount$heightSlices[-1] + 1, "m")
  
  # Create named list of metrics
  metrics <- list()
  metrics[colNames] <- sliceCount$V1
  
  return(metrics)
}


#' Recommend DEM Interpolation Method Based on Ground Point Quality
#'
#' This function analyzes a LAS object or LAS file and recommends an appropriate interpolation
#' method (`tin()`, `knnidw()`, or `kriging()`) for `lidR::rasterize_terrain()` based on
#' ground point density, ratio, and nearest-neighbor distance.
#'
#' @param las A LAS object or character path to a .las/.laz file.
#' @param res Numeric. Raster resolution (in meters) for ground point density estimation. Default is 1.
#' @param verbose Logical. If TRUE, prints diagnostic information. Default is TRUE.
#'
#' @return A character string with the recommended interpolation function (e.g., `"tin()"`)
#'
#' @details
#' This function implements a rule-based scoring system to select an appropriate terrain
#' interpolation algorithm for `lidR::rasterize_terrain()`. The recommendation is based on:
#'
#' \itemize{
#'   \item Ground point ratio (percentage of points classified as ground)
#'   \item Mean ground point density (pts/m²)
#'   \item Mean nearest-neighbor distance between ground points (meters)
#' }
#'
#' Depending on these indicators, one of the following interpolation algorithms is suggested:
#' \describe{
#'   \item{\code{"tin()"}}{Recommended when ground point distribution is dense and regular.}
#'   \item{\code{"knnidw(k = 6, p = 2)"}}{Used under intermediate conditions with moderate density.}
#'   \item{\code{"kriging(k = 10)"}}{Recommended for sparse or clustered ground points.}
#' }
#'
#' This approach follows best practices from airborne LiDAR filtering literature, including:
#' \itemize{
#'   \item Zhang et al. (2016): Cloth Simulation Filtering (CSF) – \doi{10.3390/rs8060501}
#'   \item Ma et al. (2025): Partitioned Cloth Simulation Filtering (PCSF) – \doi{10.3390/forests16071179}
#'   \item Chen et al. (2024): Adaptive DEM filtering with roughness-based interpolation – \url{https://www.sciencedirect.com/science/article/pii/S0924271624002636}
#' }
#'
#' These studies suggest that interpolation performance depends strongly on the spatial characteristics
#' of ground point clouds, especially in forested terrain. The chosen metrics are commonly used to
#' quantify LiDAR completeness and ground visibility.
#'
#' @seealso \code{\link[lidR]{rasterize_terrain}}, \code{\link[lidR]{filter_ground}}, \code{\link[lidR]{grid_density}}
#'
#' @examples
#' \dontrun{
#'   las <- readLAS("data/las/forest_tile.las")
#'   method <- recommend_dem_interpolation(las, res = 1)
#'   dem <- rasterize_terrain(las, res = 1, algorithm = eval(parse(text = method)))
#' }
#'
#' @importFrom lidR readLAS filter_ground grid_density
#' @importFrom RANN nn2
#' @export
recommend_dem_interpolation <- function(las, res = 1, verbose = TRUE) {
  if (inherits(las, "character")) las <- readLAS(las)
  if (is.empty(las)) stop("LAS file is empty or invalid")
  
  ground <- filter_ground(las)
  
  cls_tab <- table(las@data$Classification)
  n_total <- sum(cls_tab)
  n_ground <- if ("2" %in% names(cls_tab)) cls_tab["2"] else 0
  ground_pct <- 100 * as.numeric(n_ground) / n_total
  
  density_map <- grid_density(ground, res = res)
  density_vals <- values(density_map)
  density_vals <- density_vals[!is.na(density_vals)]
  mean_density <- if (length(density_vals) > 0) mean(density_vals) else 0
  
  if (nrow(ground@data) >= 2) {
    xy <- ground@data[, c("X", "Y")]
    nn_dist <- RANN::nn2(xy, k = 2)$nn.dists[, 2]
    mean_nn <- mean(nn_dist)
  } else {
    mean_nn <- Inf
  }
  
  score <- 0
  if (ground_pct > 30) score <- score + 1 else if (ground_pct < 10) score <- score - 1
  if (mean_density > 1) score <- score + 1 else if (mean_density < 0.3) score <- score - 1
  if (mean_nn < 1.5) score <- score + 1 else if (mean_nn > 3) score <- score - 1
  
  method <- if (score >= 2) {
    "tin()"
  } else if (score <= -1) {
    "kriging(k = 10)"
  } else {
    "knnidw(k = 6, p = 2)"
  }
  
  if (verbose) {
    message(sprintf("📊 Ground point ratio:     %.1f%%", ground_pct))
    message(sprintf("📊 Mean ground density:   %.2f pts/m²", mean_density))
    message(sprintf("📊 Mean NN distance:      %.2f m", mean_nn))
    message(sprintf("✅ Recommended method:    %s", method))
  }
  
  return(method)
}


#' Merge LAS/LAZ tiles into a single file using `lidR::catalog_retile`
#'
#' This function merges LAS/LAZ tiles from a directory into a single file.
#' Internally, it loads the directory as a `LAScatalog`, sets chunking to cover the entire extent,
#' and writes a single merged `.laz` or `.las` file. Uses parallel processing if desired.
#'
#' @param tile_dir Path to directory containing LAS/LAZ tiles (or a single LAS/LAZ file).
#' @param output_file Full path to the merged output file (default: `"merged.laz"`).
#' @param chunk_size Optional internal chunking for processing (default: `10000` m).
#' @param workers Number of parallel workers (default: `4`).
#'
#' @return Character string path to the created merged `.laz` file.
#'
#' @examples
#' \dontrun{
#' merge_las_tiles("tiles/", "merged.laz", workers = 6)
#' }
#'
#' @export
merge_las_tiles <- function(tile_dir,
                            output_file = "merged.laz",
                            chunk_size = 10000,
                            workers = 4) {
  if (!dir.exists(tile_dir) && !file.exists(tile_dir)) {
    stop("Input tile directory or file does not exist.")
  }
  
  library(lidR)
  library(future)
  
  set_lidr_threads(workers)
  future::plan(multisession, workers = workers)
  
  ctg <- readLAScatalog(tile_dir)
  opt_chunk_size(ctg) <- chunk_size
  opt_chunk_buffer(ctg) <- 0
  opt_output_files(ctg) <- sub("\\.la[sz]$", "", output_file)
  
  message("Merging tiles into: ", output_file)
  catalog_retile(ctg)  # This writes the file
  
  # Return final path with correct extension
  merged_path <- paste0(sub("\\.la[sz]$", "", output_file), ".las")
  return(merged_path)
}

#' Retile a LAS/LAZ file into regular tiles
#'
#' This function splits a large LAS/LAZ file into regular square tiles using `lidR::catalog_retile`.
#' It supports parallel processing and optional compression.
#'
#' @param input_file Path to the input LAS/LAZ file.
#' @param out_dir Directory to write the resulting tiles (default: `"tiles/"`).
#' @param chunk_size Numeric. Tile size in meters (default: `250`).
#' @param output_ext File extension of output tiles: either `"laz"` or `"las"` (default: `"laz"`).
#' @param buffer Buffer size between tiles, in meters (default: `0`, i.e. no overlap).
#' @param workers Number of parallel workers for processing (default: `4`).
#'
#' @return A `LAScatalog` object referencing the tiled files.
#'
#' @examples
#' \dontrun{
#' retile_las("data/input.las", out_dir = "tiles/", chunk_size = 200)
#' }
#'
#' @export
retile_las <- function(input_file,
                       out_dir = "tiles/",
                       chunk_size = 250,
                       output_ext = "laz",
                       buffer = 0,
                       workers = 4) {
  if (!file.exists(input_file)) stop("Input file not found.")
  if (!output_ext %in% c("laz", "las")) stop("Invalid extension: use 'laz' or 'las'.")
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  library(lidR)
  library(future)
  
  # Enable parallel processing
  set_lidr_threads(workers)
  future::plan(multisession, workers = workers)
  
  # Read LAS file as a catalog
  ctg <- readLAScatalog(input_file)
  opt_laz_compression(ctg) <- (output_ext == "laz")
  opt_chunk_size(ctg) <- chunk_size
  opt_chunk_buffer(ctg) <- buffer
  opt_output_files(ctg) <- file.path(out_dir, paste0("tile_{XLEFT}_{YBOTTOM}.", output_ext))
  
  message("Tiling ", input_file, " into ", out_dir, " with chunk size ", chunk_size, " m")
  tiled_ctg <- catalog_retile(ctg)
  
  return(tiled_ctg)
}



#' Preprocess ALS/TLS Point Cloud into Voxel Slice Counts
#'
#' This function filters a normalized LAS point cloud to a maximum height (`zmax`), 
#' splits the points into horizontal Z-slices (via `pointsByZSlice()`), 
#' and computes slice-wise point counts per XY pixel using `lidR::pixel_metrics()`.
#' The result is a flat data frame combining X/Y coordinates and vertical slice counts.
#'
#' @param normlas A normalized LAS object (e.g., output from `normalize_height()`) containing ground-aligned Z values.
#' @param res_xy Numeric. The horizontal resolution (in meters) of the XY voxel grid. Default is 2 m.
#' @param res_z Numeric. The vertical resolution (in meters) used for binning points into Z slices. Default is 2 m.
#'              This is passed indirectly to `pointsByZSlice()` via `zmax`.
#' @param zmax Numeric. The maximum height (in meters) to include for voxelization. Points above this value are excluded.
#'
#' @return A data frame where each row corresponds to an XY voxel column and contains:
#' \describe{
#'   \item{X, Y}{The center coordinates of the pixel (voxel column base).}
#'   \item{ground_0_1m, pulses_1_2m, ...}{Point counts per vertical slice from `pointsByZSlice()`.}
#' }
#'
#' @details
#' This function prepares voxel-based vertical profiles from a normalized LAS point cloud,
#' which is a common preprocessing step in vegetation structure analysis, such as for:
#' - Estimating Leaf Area Density (LAD)
#' - Building 3D vegetation models (e.g., for ENVI-met)
#' - Computing light extinction or aerodynamic roughness from LiDAR data
#'
#' The function performs the following steps:
#'
#' 1. **Z-Filtering**: Points are restricted to the height interval `[0, zmax]` to exclude
#'    noise (e.g., below ground level) and irrelevant outliers (e.g., birds, clouds).
#'
#' 2. **Safety Check for Empty Point Cloud**: If filtering removes all points, the function
#'    returns `NULL` to avoid errors in later processing stages.
#'
#' 3. **Dynamic Vertical Binning Limit**: The actual maximum height (`maxZ`) is computed as
#'    the minimum between the highest Z-value and the user-defined `zmax`. This ensures the
#'    binning range reflects both the data and physical modeling limits.
#'
#' 4. **Per-Pixel Vertical Slice Metrics**: Using `lidR::pixel_metrics()`, the function applies
#'    a custom-defined metric — `pointsByZSlice(Z, maxZ)` — to each XY cell. This splits the 
#'    vertical column above each pixel into 1-meter height bins (Z-slices) and counts the
#'    number of points in each slice. Empty bins are filled with 0 to ensure uniform output.
#'
#' 5. **Raster Geometry to XY Coordinates**: The function extracts the centroid (X, Y) of each
#'    pixel cell using `terra::xyFromCell()` so that slice metrics can be mapped spatially.
#'
#' 6. **Output Formatting**: The final result is a flat data frame where each row represents
#'    a voxel column. It includes the X/Y coordinate and point counts for each Z-slice,
#'    formatted with descriptive column names like `"ground_0_1m"`, `"pulses_2_3m"`, etc.
#'
#' This regularized output is designed to be compatible with downstream modeling frameworks
#' (e.g., Beer–Lambert LAD computation, ENVI-met's 3DPLANT input, or machine learning models).
#' It bridges the gap between unstructured 3D point clouds and gridded model inputs.
#'
#' @note The function assumes the input LAS object is already normalized to ground level 
#' (i.e., Z = 0 corresponds to terrain surface). Use `normalize_height()` beforehand if needed.
#'
#' @seealso 
#' [lidR::pixel_metrics()], 
#' [lidR::voxel_metrics()], 
#' [pointsByZSlice()], 
#' [normalize_height()], 
#' [terra::xyFromCell()]
#'
#' @examples
#' \dontrun{
#' las <- lidR::readLAS("path/to/normalized.laz")
#' voxel_df <- preprocess_voxels(las, res_xy = 2, res_z = 2, zmax = 40)
#' head(voxel_df)
#' }
#'
#' @export
preprocess_voxels <- function(normlas, res_xy = 2, res_z = 2, zmax = 40) {
  # Assign LAS to working variable for clarity
  las <- normlas
  
  # Step 1: Filter points to a valid height range [0, zmax]
  las <- filter_poi(las, Z >= 0 & Z <= zmax)
  
  # Step 2: Check if the LAS object is empty after filtering
  if (lidR::is.empty(las)) return(NULL)
  
  # Step 3: Define maximum vertical bin based on actual max Z (floored), capped at zmax
  maxZ <- min(floor(max(las@data$Z)), zmax)
  
  # Step 4: Define metric function for pixel-wise Z-slice counting
  func <- formula(paste0("~pointsByZSlice(Z, ", maxZ, ")"))
  
  # Step 5: Compute per-pixel vertical structure using custom Z-slice function
  voxels <- pixel_metrics(las, func, res = res_xy)
  
  # Step 6: Extract the X/Y centroid coordinates of each pixel cell
  xy <- terra::xyFromCell(voxels, seq_len(ncell(voxels)))
  
  # Step 7: Extract Z-slice point counts from each pixel cell
  vals <- terra::values(voxels)
  
  # Step 8: Combine XY coordinates with slice values into a data frame
  df <- cbind(xy, vals)
  colnames(df)[1:2] <- c("X", "Y")  # Rename coordinate columns
  
  return(df)
}


#' Convert Vertical Pulse Counts to LAD Using the Beer–Lambert Law
#'
#' This function applies the Beer–Lambert law to vertical LiDAR pulse count profiles
#' (from voxelized ALS/TLS data) to estimate Leaf Area Density (LAD) per vertical slice.
#' It normalizes the pulse density, corrects edge cases, applies the extinction formula,
#' and scales or clips LAD values to stay within biologically plausible ranges.
#'
#' @param df A data frame containing columns with pulse counts per Z-slice. These columns must be named with the `"pulses_"` prefix, as produced by `pointsByZSlice()` and `preprocess_voxels()`.
#' @param grainsize Numeric. Vertical resolution of each slice (in meters). Default is 2 m.
#' @param k Numeric. Extinction coefficient (typically between 0.3 and 0.5); default is 0.3.
#' @param scale_factor Numeric. Scaling factor applied after Beer–Lambert transformation to adjust LAD magnitude (default: 1.2).
#' @param lad_max Numeric or `NULL`. Maximum LAD value allowed (default: 2.5). Use `NULL` to disable.
#' @param lad_min Numeric or `NULL`. Minimum LAD value allowed (default: 0.05). Use `NULL` to disable.
#' @param keep_pulses Logical. If `TRUE`, the original `"pulses_"` columns are retained in the output. Default is `FALSE`.
#'
#' @return A data frame with the same structure as input, but with new `"lad_"` columns
#' for each original `"pulses_"` column, containing the estimated LAD values. 
#' Original pulse columns are optionally removed.
#'
#' @details
#' The Beer–Lambert law models **light attenuation** through a medium such as foliage. In the context of LiDAR, 
#' this relationship is inverted to estimate **leaf area density** from the relative decrease in returned pulses:
#'
#' \deqn{
#' LAD = -\frac{\log(1 - p)}{k \cdot \Delta z}
#' }
#'
#' where:
#' - \( p \) is the normalized proportion of pulses in a given voxel column and slice
#' - \( k \) is the extinction coefficient (typically 0.3–0.5)
#' - \( \Delta z \) is the vertical resolution of the slice (grainsize)
#'
#' To avoid mathematical issues:
#' - Values of \( p \geq 1 \) are clipped to 0.9999
#' - Values of \( p \leq 0 \) are clipped to 1e-5
#'
#' A **scaling factor** can be applied to tune the LAD magnitude (empirical correction).
#' LAD values are optionally clipped to a maximum and minimum for ecological realism or 
#' to avoid over-saturation artifacts in further modeling (e.g., radiative transfer, ENVI-met input).
#'
#' This approach assumes:
#' - Pulse counts per slice are proportional to occlusion probability
#' - Normalization by column maximum represents local beam extinction well enough
#' - The LiDAR pulse distribution is representative of foliage density (most valid in leaf-on conditions)
#'
#' @note
#' - For ALS data, this method provides an **approximate LAD profile** suitable for modeling but
#'   not physically exact. For more accurate LAD estimation, full waveform or TLS data is preferred.
#' - Normalization by max pulse count assumes that the densest slice corresponds to near-complete attenuation.
#'   This introduces uncertainty if canopy gaps dominate the scene.
#'
#' @seealso
#' [preprocess_voxels()], [pointsByZSlice()], [lidR::voxel_metrics()], [terra::rast()]
#'
#' @examples
#' \dontrun{
#' df_voxels <- preprocess_voxels(las)
#' df_lad <- convert_to_LAD_beer(df_voxels, grainsize = 2, k = 0.3)
#' head(df_lad)
#' }
#'
#' @export
convert_to_LAD_beer <- function(df, grainsize = 2, k = 0.3, scale_factor = 1.2,
                                lad_max = 2.5, lad_min = 0.000, keep_pulses = FALSE) {

  df_lad <- as.data.frame(df)
  # Find all columns containing vertical pulse count data
  pulse_cols <- grep("^pulses_", names(df_lad), value = TRUE)
  
  # Iterate over each pulse column (i.e., per Z-slice)
  for (col in pulse_cols) {
    # Construct name for LAD output column
    lad_col <- paste0("lad_", col)
    
    # Normalize pulse density per column (relative to max)
    p_rel <- df_lad[[col]] / max(df_lad[[col]], na.rm = TRUE)
    
    # Avoid values of 0 or 1 that would break log(1 - p)
    p_rel[p_rel >= 1] <- 0.9999
    p_rel[p_rel <= 0] <- 1e-5
    
    # Apply Beer–Lambert transformation to estimate LAD
    lad_vals <- -log(1 - p_rel) / (k * grainsize)
    
    # Apply empirical scale factor to adjust LAD magnitude
    lad_vals <- lad_vals * scale_factor
    
    # Enforce maximum and minimum LAD values (if set)
    if (!is.null(lad_max)) {
      lad_vals <- pmin(lad_vals, lad_max)
    }
    if (!is.null(lad_min)) {
      lad_vals <- pmax(lad_vals, lad_min)
    }
    
    # Store LAD values in new column
    df_lad[[lad_col]] <- lad_vals
    
    # Optionally remove original pulse column
    if (!keep_pulses) {
      df_lad[[col]] <- NULL
    }
  }
  
  return(df_lad)
}



#' Compute Canopy Traits from Long-Format LAD Profiles
#'
#' This function calculates key vegetation structure metrics from long-format 
#' Leaf Area Density (LAD) profiles, grouped by `ENVIMET_ID`. These traits 
#' include Leaf Area Index (LAI), maximum LAD, crown height, total height, 
#' leaf thickness, and aerodynamic roughness length, all of which are relevant 
#' for canopy modeling (e.g., in ENVI-met).
#'
#' @param lad_df A long-format data frame containing LAD values per voxel slice.
#' Must include at least `ENVIMET_ID`, `lad`, `z`, and optionally `Height_CHM`, `LeafThickness`, and `cluster`.
#' @param res_z Vertical voxel resolution (in meters). Default is 2 m.
#' @param LAD_cutoff Minimum LAD threshold for crown detection and trait computation (default: 0.05).
#' @param leaf_thickness_df Optional data frame mapping `cluster` to `LeafThickness` values.
#' @param roughness_fun A function to derive aerodynamic roughness length from total height.
#' Defaults to `height * 0.1`, but can be replaced with empirical models.
#' @param plantclass_prefix Character string used to label the plant class, typically `"Tree"` or `"Plant"`.
#'
#' @return A data frame with one row per `ENVIMET_ID` and the following columns:
#' \describe{
#'   \item{ENVIMET_ID}{Unique identifier for the plant object (voxel column or point).}
#'   \item{LAI}{Leaf Area Index (sum of LAD across vertical profile, clipped at 95th percentile).}
#'   \item{MaxLAD}{Maximum LAD value (95th percentile capped at 3).}
#'   \item{CrownHeight}{Height (in res_z units) of topmost voxel with LAD above cutoff.}
#'   \item{Height}{Estimated total plant height, based on CHM or max Z.}
#'   \item{LeafThickness}{Estimated or provided leaf thickness (from lookup table or column).}
#'   \item{LADcutoff}{The LAD threshold used for crown height filtering.}
#'   \item{RoughnessLength}{Estimated aerodynamic roughness length (via `roughness_fun`).}
#'   \item{plantclass}{String ID used for matching plant profile in simulations.}
#' }
#'
#' @details
#' This function summarizes vertical LAD profiles (from ALS or TLS data) into trait values 
#' that are essential for microclimate or ecological modeling. The key components:
#'
#' - **LAI**: Computed as the sum of LAD values per profile, excluding extreme values above 
#'   the 95th percentile to avoid waveform or outlier artifacts.
#'
#' - **MaxLAD**: The 95th percentile LAD value, capped at 3 to avoid unphysical spikes.
#'
#' - **CrownHeight**: The highest voxel (Z/res_z) with LAD above the cutoff threshold, interpreted 
#'   as the vertical extent of the crown (not the total plant height).
#'
#' - **Height**: Either taken from an external `Height_CHM` column (if provided), or estimated 
#'   as the maximum Z value in the profile × `res_z`.
#'
#' - **LeafThickness**: Optionally joined from a lookup table (`leaf_thickness_df`) using `cluster`.
#'
#' - **RoughnessLength**: Estimated from height using a customizable function (default: 10% of height).
#'
#' This function ensures compatibility with ENVI-met’s 3DPLANT system, where LAI, MaxLAD, and 
#' CrownHeight directly map to physical vegetation parameters in `.pld` files.
#'
#' @note
#' - CrownHeight is returned in voxel units (Z index × res_z); it is not a biologically precise
#'   crown base height but an upper limit used for model placement.
#' - The input `lad_df` must be long-format, with rows representing individual vertical slices per plant.
#'
#' @seealso
#' [convert_to_LAD_beer()], [export_lad_to_envimet3d()], [normalize_height()], [rasterize_canopy()]
#'
#' @examples
#' \dontrun{
#' traits_df <- compute_traits_from_lad(lad_long_df, res_z = 2)
#' head(traits_df)
#' }
#'
#' @export
compute_traits_from_lad <- function(lad_df, res_z = 2, LAD_cutoff = 0.05,
                                    leaf_thickness_df = NULL,
                                    roughness_fun = function(height) height * 0.1,
                                    plantclass_prefix = "Tree") {
  
  # Optional: merge leaf thickness values if available
  if (!is.null(leaf_thickness_df)) {
    lad_df <- lad_df %>%
      left_join(leaf_thickness_df %>% select(cluster, LeafThickness), by = "cluster")
  }
  
  # Split by ENVIMET_ID (each plant or tree column)
  groups <- split(lad_df, lad_df$ENVIMET_ID)
  
  # Trait calculation for each plant object
  result_list <- lapply(groups, function(group_df) {
    
    # LAD filtering: keep values under 95th percentile and above cutoff
    lad_95 <- quantile(group_df$lad, 0.95, na.rm = TRUE)
    lad_filtered <- group_df$lad[group_df$lad <= lad_95 & group_df$lad > LAD_cutoff]
    
    # LAI: sum of LAD (unit is m²/m² if LAD is in m²/m³ × height slice)
    lai <- sum(lad_filtered, na.rm = TRUE)
    
    # MaxLAD: 95th percentile capped at a max of 3
    max_lad <- min(lad_95, 3)
    
    # Crown height: highest voxel slice above LAD threshold
    if (any(group_df$lad > LAD_cutoff, na.rm = TRUE)) {
      crown_z <- max(group_df$z[group_df$lad > LAD_cutoff]/res_z, na.rm = TRUE)
      crown_height <- crown_z
    } else {
      crown_height <- NA
    }
    
    # Total plant height: prefer CHM column, fall back to max Z
    height_value <- if ("Height_CHM" %in% names(group_df) && !all(is.na(group_df$Height_CHM))) {
      unique(group_df$Height_CHM)[1]
    } else {
      max(group_df$z, na.rm = TRUE) * res_z
    }
    
    # Leaf thickness: taken from joined lookup or left NA
    lt <- if ("LeafThickness" %in% names(group_df)) {
      unique(group_df$LeafThickness)[1]
    } else {
      NA
    }
    
    # Return one row of traits
    data.frame(
      ENVIMET_ID = unique(group_df$ENVIMET_ID),
      LAI = lai,
      MaxLAD = max_lad,
      CrownHeight = crown_height,
      Height = height_value,
      LeafThickness = lt,
      LADcutoff = LAD_cutoff,
      RoughnessLength = roughness_fun(height_value),
      plantclass = paste0(plantclass_prefix, "_", unique(group_df$ENVIMET_ID))
    )
  })
  
  # Combine all into one data frame
  traits_df <- do.call(rbind, result_list)
  
  # Remove rows with missing LAI (empty plant objects)
  traits_df <- traits_df[!is.na(traits_df$LAI), ]
  
  return(traits_df)
}




#' Export Clustered LAD Profiles to ENVI-met PLANT3D (.pld) XML Format
#'
#' This function exports vertical LAD profiles (e.g., from TLS or ALS clustering) to
#' ENVI-met-compatible 3DPLANT XML (.pld) files. Each profile is mapped to a species definition,
#' associated with a seasonal LAI and blossom cycle, and written in sparse matrix format.
#'
#' @param lad_df A long-format data frame with LAD values per voxel slice and per plant.
#' Must contain columns `ENVIMET_ID`, `z`, and `lad`. Optionally includes `species_class`.
#' @param file_out Output file path for the `.pld` XML file. Default is `"tls_envimet_tree.pld"`.
#' @param res_z Vertical resolution of LAD slices (in meters). Default is 1 m.
#' @param trait_df Optional data frame of additional plant traits (unused here but compatible for future use).
#'
#' @return Writes an XML file to disk in ENVI-met PLANT3D format and prints confirmation to console.
#'
#' @details
#' This function is used to convert voxelized and clustered LAD profiles derived from TLS/ALS into
#' XML `.pld` files for use in ENVI-met's 3DPLANT vegetation modeling system.
#'
#' For each unique `ENVIMET_ID` (typically corresponding to a tree column), a separate `<PLANT3D>` block
#' is generated. Each plant object is assigned:
#' - A default or class-based species name (e.g., `"Fagus sylvatica"`)
#' - A height based on the number of vertical LAD layers
#' - LAD values in sparse 3D format
#' - A monthly season and blossom profile (12 months)
#'
#' The function uses a predefined lookup table (`custom_profiles`) for assigning seasonal LAI and blossom 
#' curves to known species. Unmapped species default to generic "broadleaf" or "conifer" profiles.
#'
#' Species-specific traits such as `Albedo`, `Transmittance`, `RootDiameter`, and `LeafType` are hard-coded
#' per class ID but can be extended as needed.
#'
#' @note
#' - `species_class` must be numeric and match predefined mappings. Unknown classes fall back to a default.
#' - LAD values are averaged across horizontal slices (X/Y = 1) and rounded to 5 decimals.
#' - The LAD format is `"sparematrix-3D"` with only z-direction layers, suitable for columnar tree shapes.
#'
#' @seealso
#' [compute_traits_from_lad()], [convert_to_LAD_beer()], [ENVI-met documentation on PLANT3D profiles]
#'
#' @examples
#' \dontrun{
#' export_lad_to_envimet_p3d(lad_df = my_lad_data, file_out = "trees.pld", res_z = 1)
#' }
#'
#' @export
export_lad_to_envimet_p3d <- function(lad_df, file_out = "tls_envimet_tree.pld", res_z = 1, trait_df = NULL) {
  # --- Define seasonal LAI/blossom profiles per species ---
  season_profiles <- list(
    broadleaf = list(
      Season = c(0.3, 0.3, 0.3, 0.4, 0.7, 1, 1, 1, 0.8, 0.6, 0.3, 0.3),
      Blossom = c(0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    conifer = list(
      Season = rep(0.8, 12),
      Blossom = rep(0, 12)
    )
  )
  
  # --- Assign species-specific profiles if known ---
  custom_profiles <- list(
    "Fagus sylvatica" = season_profiles$broadleaf,
    "Alnus glutinosa" = season_profiles$broadleaf,
    "Fraxinus excelsior" = season_profiles$broadleaf,
    "Pseudotsuga menziesii" = season_profiles$conifer,
    "Larix decidua" = season_profiles$broadleaf,
    "Quercus robur" = season_profiles$broadleaf,
    "Picea abies" = list(Season = rep(0.75, 12), Blossom = rep(0, 12))
  )
  
  # --- Clean input and encode Z-layer index ---
  lad_df <- lad_df[!is.na(lad_df$lad), ]
  z_map <- setNames(seq_along(sort(unique(lad_df$z))), sort(unique(lad_df$z)))
  lad_df$k <- z_map[as.character(lad_df$z)]
  lad_df$lad_value <- round(lad_df$lad, 5)
  
  # --- Initialize XML document ---
  tree_ids <- unique(lad_df$ENVIMET_ID)
  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  root <- newXMLNode("ENVI-MET_Datafile")
  header_node <- newXMLNode("Header")
  addChildren(header_node, newXMLNode("filetype", "DATA"))
  addChildren(header_node, newXMLNode("version", "1"))
  addChildren(header_node, newXMLNode("revisiondate", now))
  addChildren(header_node, newXMLNode("remark", "Clustered TLS-based trees"))
  addChildren(header_node, newXMLNode("fileInfo", "Clustered LAD Trees"))
  addChildren(header_node, newXMLNode("checksum", "32767"))
  addChildren(header_node, newXMLNode("encryptionlevel", "1699612"))
  addChildren(root, header_node)
  
  # --- Iterate over all unique tree profiles ---
  for (id in tree_ids) {
    tree_df <- lad_df[lad_df$ENVIMET_ID == id, ]
    
    # Average LAD across XY slices per height layer
    profile <- tree_df %>%
      group_by(z = k) %>%
      summarise(lad_value = mean(lad_value, na.rm = TRUE), .groups = "drop")
    
    # Define voxel profile parameters
    zlayers <- max(profile$z)
    dataI <- 1
    dataJ <- 1
    Height <- zlayers * res_z
    
    # --- Default plant parameters ---
    name <- "Fagus sylvatica"; albedo <- 0.18; trans <- 0.30; root_d <- 4.5; leaf_type <- 1
    
    # --- Override species traits if species_class is available ---
    if (!all(is.na(tree_df$species_class))) {
      class_val <- na.omit(unique(tree_df$species_class))[1]
      if (class_val == 2)      { name <- "Alnus glutinosa";       albedo <- 0.18; trans <- 0.35; root_d <- 3.5; leaf_type <- 1 }
      else if (class_val == 3) { name <- "Fraxinus excelsior";    albedo <- 0.19; trans <- 0.38; root_d <- 4.0; leaf_type <- 1 }
      else if (class_val == 4) { name <- "Fagus sylvatica";       albedo <- 0.18; trans <- 0.30; root_d <- 4.5; leaf_type <- 1 }
      else if (class_val == 5) { name <- "Pseudotsuga menziesii"; albedo <- 0.20; trans <- 0.18; root_d <- 4.2; leaf_type <- 2 }
      else if (class_val == 6) { name <- "Larix decidua";         albedo <- 0.23; trans <- 0.25; root_d <- 4.0; leaf_type <- 2 }
      else if (class_val == 7) { name <- "Quercus robur";         albedo <- 0.20; trans <- 0.35; root_d <- 5.0; leaf_type <- 1 }
      else if (class_val == 11){ name <- "Picea abies";           albedo <- 0.22; trans <- 0.15; root_d <- 3.0; leaf_type <- 2 }
    }
    
    # --- Assign seasonal profile ---
    profile_key <- name
    season_profile <- custom_profiles[[profile_key]]
    if (is.null(season_profile)) {
      season_profile <- season_profiles[[ifelse(leaf_type == 1, "broadleaf", "conifer")]]
    }
    
    # --- Build PLANT3D block ---
    plant_node <- newXMLNode("PLANT3D")
    addChildren(plant_node, newXMLNode("ID", id))
    addChildren(plant_node, newXMLNode("Description", "Clustered TLS Tree"))
    addChildren(plant_node, newXMLNode("AlternativeName", name))
    addChildren(plant_node, newXMLNode("Planttype", "0"))
    addChildren(plant_node, newXMLNode("Leaftype", as.character(leaf_type)))
    addChildren(plant_node, newXMLNode("Albedo", sprintf("%.5f", albedo)))
    addChildren(plant_node, newXMLNode("Eps", "0.96000"))
    addChildren(plant_node, newXMLNode("Transmittance", sprintf("%.5f", trans)))
    addChildren(plant_node, newXMLNode("Height", sprintf("%.5f", Height)))
    addChildren(plant_node, newXMLNode("Width", "1.00000"))
    addChildren(plant_node, newXMLNode("Depth", "1.00000"))
    addChildren(plant_node, newXMLNode("RootDiameter", sprintf("%.5f", root_d)))
    addChildren(plant_node, newXMLNode("cellsize", sprintf("%.5f", res_z)))
    addChildren(plant_node, newXMLNode("xy_cells", dataI))
    addChildren(plant_node, newXMLNode("z_cells", zlayers))
    addChildren(plant_node, newXMLNode("scalefactor", "1.00000"))
    
    # --- LAD Profile block (sparse 3D format) ---
    lad_lines <- apply(profile, 1, function(r) {
      sprintf("%d,%d,%d,%.5f", 1, 1, r[1], r[2])
    })
    lad_node <- newXMLNode("LAD-Profile",
                           attrs = c(type = "sparematrix-3D",
                                     dataI = dataI,
                                     dataJ = dataJ,
                                     zlayers = zlayers,
                                     defaultValue = "0.00000"),
                           .children = paste(lad_lines, collapse = "\n"))
    addChildren(plant_node, lad_node)
    
    # --- Add seasonality and flowering profiles ---
    addChildren(plant_node, newXMLNode("Season-Profile",
                                       paste(sprintf("%.5f", season_profile$Season), collapse = ",")))
    addChildren(plant_node, newXMLNode("Blossom-Profile",
                                       paste(sprintf("%.5f", season_profile$Blossom), collapse = ",")))
    
    # --- Disable L-System generation ---
    addChildren(plant_node, newXMLNode("L-SystemBased", "0"))
    addChildren(root, plant_node)
  }
  
  # --- Write XML to file ---
  saveXML(root, file = file_out, indent = TRUE, encoding = "UTF-8")
  message("✔ ENVI-met PLANT3D (.pld) written to: ", normalizePath(file_out))
}

#' Generate Compact ENVI-met-Compatible Base36 String IDs
#'
#' Converts an integer index to a left-padded base-36 string (digits + uppercase A–Z),
#' prefixed with `"S"`, suitable for use as compact `ENVIMET_ID`s (e.g., `"S0001A"`).
#' Useful when assigning unique but readable identifiers for synthetic plant elements
#' in 3D simulation domains.
#'
#' @param n An integer (scalar) to convert to base-36.
#' @param width Integer width of the resulting code (default: 5). Strings are zero-padded on the left.
#'
#' @return A character string in base-36 representation, prefixed with `"S"` and left-padded
#' to match the desired `width`. Returns e.g. `"S0000A"`, `"S0001Z"`, `"S00010"`, etc.
#'
#' @details
#' Base-36 encoding uses the characters `0–9` and `A–Z` to represent numbers in a compact alphanumeric form.
#' This is commonly used in modeling frameworks like ENVI-met where short string IDs are needed to:
#' - Uniquely label plant objects (`ENVIMET_ID`)
#' - Avoid numeric-only names (which may conflict with XML or database formats)
#' - Allow for high capacity in short formats (36⁵ = ~60 million unique IDs for `width = 5`)
#'
#' The function:
#' 1. Converts the number `n` to base-36 using digit/modulo arithmetic
#' 2. Left-pads the result to fixed width with zeros (`"0"`)
#' 3. Adds a prefix `"S"` to ensure the string starts with a non-numeric character
#'
#' @note
#' This function assumes positive integers (`n > 0`). No validation is done for negative or non-integer input.
#' It is up to the user to avoid duplicate IDs.
#'
#' @seealso
#' [export_lad_to_envimet_p3d()], [sprintf()], [strtoi()] for inverse conversion
#'
#' @examples
#' int_to_base36(1)      # "S00001"
#' int_to_base36(35)     # "S0000Z"
#' int_to_base36(36)     # "S00010"
#' int_to_base36(12345)  # "S009IX"
#'
#' @export
int_to_base36 <- function(n, width = 5) {
  # Base-36 character set: 0–9, A–Z
  chars <- c(0:9, LETTERS)
  base <- length(chars)
  
  # Convert to base-36 via division/remainder
  result <- character()
  while (n > 0) {
    result <- c(chars[(n %% base) + 1], result)
    n <- n %/% base
  }
  
  # Collapse to single string
  result <- paste(result, collapse = "")
  
  # Pad with zeros on the left to match width
  padded <- sprintf(paste0("%0", width, "s"), result)
  
  # Replace spaces with "0" and add "S" prefix
  paste0("S", substr(gsub(" ", "0", padded), 1, width))
}


#' Suggest Optimal Number of Principal Components
#'
#' This function evaluates a PCA object and recommends the number of components to retain
#' based on three common criteria: cumulative explained variance, Kaiser criterion, and the
#' elbow method. Optionally, a diagnostic scree plot and cumulative variance plot are shown.
#'
#' @param pca_obj A PCA object as returned by [prcomp()], which must contain the `sdev` vector.
#' @param variance_cutoff The cumulative explained variance threshold to use for the primary selection (default: 0.8).
#' @param plot Logical. If `TRUE`, produces a scree plot and a cumulative variance plot for visual inspection.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{n_pcs}{Recommended number of principal components based on variance cutoff.}
#'   \item{explained_variance}{Number of components needed to reach the cutoff.}
#'   \item{kaiser}{Number of components with eigenvalue > 1 (Kaiser criterion).}
#'   \item{elbow}{Index of the "elbow" point in the scree plot.}
#'   \item{info_table}{A summary data frame showing results for all three criteria.}
#' }
#'
#' @details
#' The number of principal components can be selected using multiple heuristics:
#'
#' - **Cumulative Explained Variance**: Retain the smallest number of components such that the
#' cumulative variance exceeds `variance_cutoff`. This is often considered the primary criterion.
#'
#' - **Kaiser Criterion**: Retain all components with eigenvalue > 1. Assumes input data has been scaled.
#'
#' - **Elbow Method**: Identifies the point where the marginal gain in explained variance drops off,
#' i.e., the inflection point of the scree plot. Here approximated by the first minimum in the
#' differences of explained variance.
#'
#' The function outputs all three estimates but uses only the **variance cutoff** for final selection.
#'
#' @note
#' The scree plot assumes components are sorted by variance (as returned by [prcomp()]).
#' The elbow method used here is a simplified heuristic and may not capture complex knees.
#'
#' @seealso
#' [prcomp()], [factoextra::fviz_eig()], [psych::principal()]
#'
#' @examples
#' pca <- prcomp(USArrests, scale. = TRUE)
#' suggest_n_pcs(pca, variance_cutoff = 0.9)
#'
#' @export
suggest_n_pcs <- function(pca_obj, variance_cutoff = 0.8, plot = TRUE) {
  # Extract standard deviations and calculate explained variance
  std_dev <- pca_obj$sdev
  var_explained <- std_dev^2 / sum(std_dev^2)
  cum_var_explained <- cumsum(var_explained)
  
  # Criterion 1: variance cutoff
  n_var <- which(cum_var_explained >= variance_cutoff)[1]
  
  # Criterion 2: Kaiser criterion (eigenvalue > 1)
  eigenvalues <- std_dev^2
  n_kaiser <- sum(eigenvalues > 1)
  
  # Criterion 3: elbow (minimum change in explained variance)
  diffs <- diff(var_explained)
  elbow <- which.min(diffs)[1]
  
  # Summarize selection criteria
  info_table <- data.frame(
    Criterion = c("Variance Cutoff", "Kaiser Criterion", "Elbow Method"),
    Components = c(n_var, n_kaiser, elbow),
    Cumulative_Variance = c(
      round(cum_var_explained[n_var], 3),
      round(cum_var_explained[n_kaiser], 3),
      round(cum_var_explained[elbow], 3)
    )
  )
  
  # Final decision based on variance criterion only
  n_final <- n_var
  
  # Output summary
  cat("📊 PCA Component Selection Summary:\n")
  print(info_table, row.names = FALSE)
  cat("\n✅ Recommended number of PCs (variance_cutoff =", variance_cutoff, "):", n_final, "\n")
  
  # Optional visualization
  if (plot) {
    par(mfrow = c(1, 2))
    
    # Scree plot
    plot(var_explained, type = "b", pch = 19, col = "steelblue",
         xlab = "Component", ylab = "Explained Variance",
         main = "Scree Plot")
    abline(h = 1, col = "red", lty = 2)
    abline(v = elbow, col = "darkgreen", lty = 3)
    legend("topright", legend = c("Kaiser (λ > 1)", "Elbow"),
           col = c("red", "darkgreen"), lty = c(2, 3), bty = "n")
    
    # Cumulative variance plot
    plot(cum_var_explained, type = "b", pch = 19, col = "darkorange",
         xlab = "Component", ylab = "Cumulative Variance",
         main = "Cumulative Explained Variance")
    abline(h = variance_cutoff, col = "red", lty = 2)
    abline(v = n_var, col = "blue", lty = 3)
    legend("bottomright", legend = c("Cutoff", "Selected Components"),
           col = c("red", "blue"), lty = c(2, 3), bty = "n")
    
    par(mfrow = c(1, 1))
  }
  
  # Return results invisibly
  invisible(list(
    n_pcs = n_final,
    explained_variance = n_var,
    kaiser = n_kaiser,
    elbow = elbow,
    info_table = info_table
  ))
}

