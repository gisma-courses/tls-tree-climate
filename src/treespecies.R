#------------------------------------------------------------------------------
# Script Name: 30_filter_classified_species_map.R
# Author: Chris Reudenbach
# Description: Cleans classified tree species raster using OTB and contextual rules
#------------------------------------------------------------------------------

# === Libraries ===
library(terra)
library(RColorBrewer)
library(link2GI)
library(rprojroot)
library(mapview)
library(dplyr)
library(tools)

# === Environment and paths ===
root_folder <- find_rstudio_root_file()
otb <- link2GI::linkOTB(searchLocation = "~/apps/OTB-9.1.0-Linux/")
Sys.setenv(OTB_APPLICATION_PATH = file.path(dirname(as.character(otb$pathOTB)), "lib/otb/applications"))
Sys.setenv(PATH = paste(otb$pathOTB, Sys.getenv("PATH"), sep = ":"))

# === Parameters ===
target_res <- res_xy  # define outside or replace with numeric (e.g. 1)
epsg <- 25832
sapflow_ext <- raster::extent(477500, 478218, 5631730, 5632500)

# === Static file names ===
pred_rds                <- "data/aerial/sfprediction_ffs_5-25_MOF_rgb.rds"
pred_tif                <- "data/aerial/prediction_ffs.tif"
majority_tif            <- "data/aerial/majority_out.tif"
aggregate_tif           <- "data/aerial/aggregate.tif"
aggregate_mode_tif      <- "data/aerial/aggregate_mode.tif"
treespecies_cleaned_tif <- "data/aerial/treespecies_cleaned.tif"

# === Class ID legend ===
ts <- data.frame(
  ID = 1:12,
  value = c(
    "agriculture", "alder", "ash", "beech", "douglas_fir",
    "larch", "oak", "pastures", "roads", "settlements",
    "spruce", "water"
  )
)

#------------------------------------------------------------------------------
# FUNCTION: Replace Douglas-fir if surrounded by Beech/Oak
#------------------------------------------------------------------------------
replace_douglas_in_buche_eiche <- function(rast_input,
                                           mode_raster_path,
                                           output_path,
                                           window_size = 5,
                                           douglas_value = 5,
                                           target_values = c(4, 7)) {
  if (window_size %% 2 == 0)
    stop("window_size must be odd")
  
  # OTB: majority mode filter
  cmr <- parseOTBFunction("ClassificationMapRegularization", otb)
  cmr$io.in <- aggregate_tif
  cmr$io.out <- mode_raster_path
  cmr$ip.radius <- as.character((window_size - 1) / 2)
  cmr$progress <- "true"
  runOTB(cmr, gili = otb$pathOTB, quiet = FALSE)
  
  r_mode <- rast(mode_raster_path)
  is_douglas <- rast_input == douglas_value
  is_oak_beech_mode <- r_mode %in% target_values
  r_new <- rast_input
  r_new[is_douglas & is_oak_beech_mode] <- r_mode[is_douglas & is_oak_beech_mode]
  
  writeRaster(r_new, output_path, overwrite = TRUE)
  return(r_new)
}

#------------------------------------------------------------------------------
# STEP 1: Load and export predicted species raster
#------------------------------------------------------------------------------
sapflow_species <- readRDS(pred_rds)
raster::writeRaster(sapflow_species, pred_tif, progress = "text", overwrite = TRUE)

sapflow_species <- raster::crop(sapflow_species, sapflow_ext)
raster::writeRaster(sapflow_species, "data/aerial/prediction_ffs_cut.tif", overwrite = TRUE)

#------------------------------------------------------------------------------
# STEP 2: OTB modal filter (3x3 smoothing)
#------------------------------------------------------------------------------
cmr <- parseOTBFunction("ClassificationMapRegularization", otb)
cmr$io.in <- pred_tif
cmr$io.out <- majority_tif
cmr$progress <- "true"
cmr$ip.radius <- "1"
filter_treespecies <- runOTB(cmr, gili = otb$pathOTB, quiet = FALSE, retRaster = TRUE)

#------------------------------------------------------------------------------
# STEP 3: Aggregate to 1 m resolution (median)
#------------------------------------------------------------------------------
r <- rast(majority_tif)
cur_res <- res(r)[1]
fact <- round(target_res / cur_res)
if (target_res <= cur_res)
  stop("Zielauflösung ist kleiner als aktuelle.")

r_agg <- aggregate(r, fact = fact, fun = median, na.rm = TRUE)
writeRaster(r_agg, aggregate_tif, overwrite = TRUE)

#------------------------------------------------------------------------------
# STEP 4: Contextual cleanup: Douglas-fir correction
#------------------------------------------------------------------------------
species_cleaned <- replace_douglas_in_buche_eiche(
  rast_input = r_agg,
  mode_raster_path = aggregate_mode_tif,
  output_path = treespecies_cleaned_tif,
  window_size = 9
)

#------------------------------------------------------------------------------
# STEP 5: Interactive comparison
#------------------------------------------------------------------------------
if (visualize)
{palette <- brewer.pal(12, "Paired")
zoom_center <- list(lng = 8.68443, lat = 50.84089, zoom = 18)

m1 <- mapview(crop(rast(pred_tif), sapflow_ext), col.regions = palette, at = ts$ID, layer.name = "Species 0.2m")
m2 <- mapview(crop(rast(majority_tif), sapflow_ext), col.regions = palette, at = ts$ID, layer.name = "3x3 modal_filt")
m3 <- mapview(crop(rast(aggregate_tif), sapflow_ext), col.regions = palette, at = ts$ID, layer.name = "Aggregated 1m")
m4 <- mapview(crop(rast(treespecies_cleaned_tif), sapflow_ext), col.regions = palette, at = ts$ID, layer.name = "Douglas out 1m")

out <- sync(
  m1@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom),
  m2@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom),
  m3@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom),
  m4@map %>% leaflet::setView(zoom_center$lng, zoom_center$lat, zoom_center$zoom)
)

out
}