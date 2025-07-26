#------------------------------------------------------------------------------
# Script Type: Processing script
# Script Name: 30_filter_classified_species_map.R
# Author: Chris Reudenbach, creuden@gmail.com
#
# Description:
#   - Loads a classified raster map of tree species
#   - Applies spatial smoothing using OTB ClassificationMapRegularization
#   - Aggregates to 1 m resolution using median filtering
#   - Performs contextual correction: replaces isolated Douglas-fir pixels
#     with Beech or Oak if those are dominant in the local neighborhood
#
# Input:
#   - RDS and GeoTIFF classification of tree species (0.2 m resolution)
#
# Output:
#   - Cleaned and aggregated species raster (1 m resolution)
#
# Dependencies:
#   - OTB 9.1+ with PATH and OTB_APPLICATION_PATH correctly set
#
# Copyright: Chris Reudenbach 2021, GPL (>= 3)
# Git: https://github.com/gisma-courses/microclimate.git
#
# Commentary:
#   This workflow separates two crucial but distinct steps in map cleaning:
#   1. **Noise reduction (smoothing):** Applied via OTB's ClassificationMapRegularization.
#      - It performs a fast majority-based smoothing using a local moving window (e.g., 3x3).
#      - This step removes small speckles or misclassified pixels in homogeneous areas.
#      - It is computationally efficient due to OTB's C++-based implementation.
#
#   2. **Semantic filtering:** Performed in R via a contextual reclassification function.
#      - Specifically targets ecologically unlikely or isolated Douglas-fir pixels.
#      - These are replaced with surrounding Beech or Oak pixels if they dominate locally.
#      - Allows flexible, rule-based filtering that OTB cannot natively perform.
#
#   ➤ Both steps are technically possible in R using terra::focal(), but:
#     - Smoothing with `focal()` is **much slower** on large rasters (single-threaded).
#     - OTB is highly recommended for performance.
#
#   ➤ The R-based semantic filtering step is **required** if logical replacement
#     rules (like Douglas-fir substitution) are needed. This goes beyond statistical smoothing.
#------------------------------------------------------------------------------


# === Libraries ===
library(terra)           # raster handling
library(RColorBrewer)    # color palettes
library(link2GI)         # OTB integration
library(rprojroot)
library(tools)           # file name tools
library(mapview)         # interactive maps
library(dplyr)           # data manipulation

# === Environment and paths ===
root_folder <- find_rstudio_root_file()

# Set up OTB environment
otb <- link2GI::linkOTB(searchLocation = "~/apps/OTB-9.1.0-Linux/")
Sys.setenv(OTB_APPLICATION_PATH = file.path(dirname(as.character(otb$pathOTB)), "lib/otb/applications"))
Sys.setenv(PATH = paste(otb$pathOTB, Sys.getenv("PATH"), sep = ":"))

# === Parameters ===
target_res <- 1                # desired resolution in meters
min_tree_height <- 2           # (not used yet)
fn <- "5-25_MOF_rgb"           # image stem
epsg <- 25832                  # UTM32N
sapflow_ext <- raster::extent(477500, 478218, 5631730, 5632500)  # area of interest

# === Class ID legend ===
ts <- data.frame(
  ID = 1:12,
  value = c(
    "agriculture",
    "alder",
    "ash",
    "beech",
    "douglas_fir",
    "larch",
    "oak",
    "pastures",
    "roads",
    "settlements",
    "spruce",
    "water"
  )
)

#------------------------------------------------------------------------------
# FUNCTION: Replace isolated Douglas-fir with Beech or Oak if dominant around
#------------------------------------------------------------------------------
replace_douglas_in_buche_eiche <- function(rast_input,
                                           window_size = 5,
                                           douglas_value = 5,
                                           target_values = c(4, 7),
                                           target_res = 1.0) {
  if (window_size %% 2 == 0)
    stop("window_size must be odd")
  
  # Focal window matrix (square)
  w <- matrix(1, nrow = window_size, ncol = window_size)
  
  # Run OTB ClassificationMapRegularization to compute local mode
  cmr <- parseOTBFunction("ClassificationMapRegularization", otb)
  cmr$io.in <- sprintf("data/aerial/%s_%sm.tif",
                       tools::file_path_sans_ext(basename("data/aerial/aggregate.tif")),
                       target_res)
  cmr$io.out <- sprintf("data/aerial/%s_%sm.tif",
                        tools::file_path_sans_ext(basename("data/aerial/aggregate_mode.tif")),
                        window_size)
  cmr$ip.radius <- as.character((window_size - 1) / 2)  # for 5x5 window: radius = (5 - 1)/2 = 2
  cmr$progress <- "true"
  
  runOTB(cmr, gili = otb$pathOTB, quiet = FALSE)

  # Identify Douglas-fir pixels and surrounding Beech/Oak dominance
  r_mode = rast(cmr$io.out)
  rast_input = rast(cmr$io.in)
  is_douglas <- rast_input == douglas_value
  is_oak_beech_mode <- r_mode %in% target_values
  replace_mask <- is_douglas & is_oak_beech_mode
  
  # Replace Douglas-fir where Beech or Oak dominate
  r_new <- rast_input
  r_new[replace_mask] <- r_mode[replace_mask]
  
  # Construct output path 
  outname <- paste0("data/aerial/",
                    "agg_cleand_",
                    as.character(target_res),
                    "m.tif")
  writeRaster(r_new, outname,overwrite = TRUE)
  
  return(r_new)
}

#------------------------------------------------------------------------------
# STEP 1: Read tree species classification from RDS
#------------------------------------------------------------------------------
sapflow_species <- readRDS("data/aerial/sfprediction_ffs_5-25_MOF_rgb.rds")

# Write to GeoTIFF for further processing
raster::writeRaster(
  sapflow_species,
  "data/aerial/prediction_ffs.tif",
  progress = "text",
  overwrite = TRUE
)

# Crop to sapflow test area
sapflow_species <- raster::crop(sapflow_species, sapflow_ext)
raster::writeRaster(
  sapflow_species,
  "data/aerial/prediction_ffs_cut.tif",
  progress = "text",
  overwrite = TRUE
)

#------------------------------------------------------------------------------
# STEP 2: Run OTB ClassificationMapRegularization (majority filter)
#------------------------------------------------------------------------------
cmr <- parseOTBFunction("ClassificationMapRegularization", otb)
cmr$io.in <- "data/aerial/prediction_ffs.tif"
cmr$io.out <- "data/aerial/majority_out.tif"
cmr$progress <- "true"
cmr$ip.radius <- "1"

filter_treespecies <- runOTB(cmr,
                             gili = otb$pathOTB,
                             quiet = FALSE,
                             retRaster = TRUE)

#------------------------------------------------------------------------------
# STEP 3: Aggregate to 1 m resolution using median
#------------------------------------------------------------------------------
r <- rast("data/aerial/majority_out.tif")
cur_res <- res(r)[1]
fact <- round(target_res / cur_res)

if (target_res <= cur_res)
  stop("Zielauflösung ist kleiner als aktuelle.")

r_agg <- aggregate(r,
                   fact = fact,
                   fun = median,
                   na.rm = TRUE)

# Build automatic filename
outfile <- sprintf("data/aerial/%s_%sm.tif",
                   tools::file_path_sans_ext(basename("data/aerial/aggregate.tif")),
                   target_res)

# Save aggregated raster
writeRaster(r_agg, outfile, overwrite = TRUE)

#------------------------------------------------------------------------------
# STEP 4: Clean Douglas-fir patches contextually
#------------------------------------------------------------------------------
species_cleaned <- replace_douglas_in_buche_eiche(window_size = 9)

#------------------------------------------------------------------------------
# STEP 5: Visualize intermediate steps (interactive)
#------------------------------------------------------------------------------


library(mapview)
library(leafsync)
library(htmlwidgets)
library(terra)
library(RColorBrewer)

# Define common parameters
palette <- brewer.pal(12, "Paired")
zoom_center <- list(lng = 8.68443, lat = 50.84089, zoom = 18)

# -- Map 1: Raw species classification (0.2 m)
m1 <- mapview(
  terra::crop(sapflow_species, sapflow_ext),
  col.regions = palette,
  at = ts$ID,
  layer.name = "Species 0.2m"
)

# -- Map 2: OTB 3×3 modal smoothing
m2 <- mapview(
  terra::crop(filter_treespecies, sapflow_ext),
  col.regions = palette,
  at = ts$ID,
  fgb = TRUE,
  layer.name = "3x3 modal_filt"
)

# -- Map 3: Aggregated to 1 m resolution
m3 <- mapview(
  terra::crop(r_agg, sapflow_ext),
  col.regions = palette,
  at = ts$ID,
  layer.name = "Aggregated 1m"
)

# -- Map 4: Douglas-fir replaced by contextual rules
m4 <- mapview(
  terra::crop(species_cleaned, sapflow_ext),
  col.regions = palette,
  at = ts$ID,
  fgb = TRUE,
  layer.name = "Douglas out 1m"
)

# Convert to leaflet and apply zoom center
lm1 <- m1@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom)
lm2 <- m2@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom)
lm3 <- m3@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom)
lm4 <- m4@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom)

# Synchronize maps side-by-side
sync(lm1, lm2, lm3, lm4)

