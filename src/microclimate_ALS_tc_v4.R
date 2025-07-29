#' --- ENVI-met 3DPLANT column generator from voxelized ALS data ---
#' Full pipeline from ALS voxelization to XML export for ENVI-met 3DPLANT.
#' Includes normalization, topographic metrics, LAD calculation, clustering, and XML writing.

# === Load Libraries ===
library(lidR)
library(terra)
library(dplyr)
library(tidyr)
library(sf)
library(here)
library(XML)
library(stats)
library(tibble)
library(rprojroot)
library(tools)
library(RANN)
library(clusternomics)
library(e1071)
library(entropy)
library(NbClust)
library(matrixStats)

# === Configuration ===
ts <- data.frame(
  ID = 1:12,
  value = c("agriculture", "alder", "ash", "beech", "douglas_fir", "larch",
            "oak", "pastures", "roads", "settlements", "spruce", "water")
)
valid_species <- c("alder", "ash", "beech", "douglas_fir", "larch", "oak", "spruce")
valid_ids <- ts$ID[ts$value %in% valid_species]

visualize <- FALSE
las_file <- here("data/ALS/tiles/")
res_xy <- 2
res_z <- 2
k <- 0.3
scale_factor <- 1.2
crs_code <- 25832
output_gpkg <- "data/envimet/envimet_p3dtree_points.gpkg"
xml_output_file <- "data/envimet/als_envimet_trees.pld"
species_raster <- rast("data/aerial/treespecies_cleaned.tif")

dir.create("data/output", showWarnings = FALSE, recursive = TRUE)
source("src/new_utils.R")

# === ⬛ STAGE: LAS Merging and Normalization ===
las_fn <- merge_las_tiles(las_file, "data/ALS/merged_output.laz", chunk_size = 10000, workers = 6)
las <- readLAS(las_fn)
crs(las) <- "EPSG:25832"

# === ⬛ STAGE: Ground Classification and DEM/CHM Generation ===
# ⚠ MEMORY: large point cloud + rasters
chm_pre <- rasterize_canopy(las, res = res_xy, algorithm = pitfree(c(0,1,3,6,9,12,16)))
rugosity <- terra::focal(chm_pre, w = matrix(1, 3, 3), fun = sd, na.rm = TRUE)
mean_rug <- global(rugosity, fun = "mean", na.rm = TRUE)[[1]]

csf_params <- if (mean_rug > 1) {
  message("Detected complex/dense canopy – using fine CSF settings")
  csf(cloth_resolution = 0.5, rigidness = 2, class_threshold = 0.4, iterations = 800)
} else {
  message("Detected open canopy – using coarse CSF settings")
  csf(cloth_resolution = 1.5, rigidness = 4, class_threshold = 0.6, iterations = 300)
}

las <- classify_ground(las, csf_params)
dem_algo <- eval(parse(text = recommend_dem_interpolation(las, res_xy)))
dem <- rasterize_terrain(las, res = res_xy, algorithm = dem_algo)
las_norm <- normalize_height(las, algorithm = knnidw(k = 6L, p = 2))

# ✅ Free original LAS
las <- NULL; gc()

# === ⬛ STAGE: CHM + DSM + Topographic Metrics ===
pit_algo <- pitfree(c(0, 1, 3, 6, 9, 12, 16))
chm <- rasterize_canopy(las_norm, res = res_xy, algorithm = pit_algo)
dsm <- rasterize_canopy(las_norm, res = res_xy, algorithm = p2r())
slope <- terrain(dem, "slope", unit = "radians")
aspect <- terrain(dem, "aspect", unit = "degrees")
TPI <- terra::focal(terrain(dsm, "TPI"), w = matrix(1,3,3), fun = mean)
TPI[TPI < 0] <- -1; TPI[TPI > 0] <- 1

topo <- c(dem, dsm, chm, slope, aspect, TPI)
names(topo) <- c("dem", "dsm", "chm", "slope", "aspect", "TPI")
writeRaster(topo, "data/ALS/topo_stack.tif", overwrite = TRUE)
plot(topo)

# === ⬛ STAGE: Voxelization + LAD ===
# ⚠ MEMORY: voxel and LAD matrices can be huge
voxels <- preprocess_voxels(las_norm, res_xy, res_z)
lad_df <- convert_to_LAD_beer(voxels, grainsize = res_z, k = k, scale_factor = scale_factor)

# ✅ Free voxelized LAS
las_norm <- NULL; voxels <- NULL; gc()

# === ⬛ STAGE: Add topographic & species data ===
topo_stack <- rast("data/ALS/topo_stack.tif")
dem <- topo_stack[["dem"]]
dsm <- topo_stack[["dsm"]]

# Erzeuge ein sf-Objekt aus lad_df
sf_lad <- st_as_sf(as.data.frame(lad_df), coords = c("X", "Y"), crs = crs_code)
geom_only <- st_geometry(sf_lad)
lad_df$elev <- exactextractr::exact_extract( dem, st_buffer(geom_only, dist = 0.1), "mean")
lad_df$slope  <- exactextractr::exact_extract( terrain(dem, "slope"), st_buffer(geom_only, dist = 0.1), "mean")
lad_df$aspect <- exactextractr::exact_extract( terrain(dem, "aspect"), st_buffer(geom_only, dist = 0.1), "mean")
lad_df$TPI    <- exactextractr::exact_extract( terrain(dem, "TPI"), st_buffer(geom_only, dist = 0.1), "mean")
lad_df$CHM    <- exactextractr::exact_extract( dsm-dem, st_buffer(geom_only, dist = 0.1), "mean")
lad_df$species_class <- exactextractr::exact_extract( species_raster, st_buffer(geom_only, dist = 0.1), "mean")

# ✅ Remove sf and stack
sf_lad <- NULL; topo_stack <- NULL;geom_only=NULL; gc()


# === ⬛ STAGE: Compute structural indices ===
# Only keep numeric LAD columns

lad_matrix <- select(lad_df, starts_with("lad_")) |> as.matrix()
rownames(lad_matrix) <- paste(lad_df$X, lad_df$Y, sep = "_")

# More efficient summary statistics
lad_df$LAD_mean        <- matrixStats::rowMeans2(lad_matrix, na.rm = TRUE)
lad_df$LAD_max         <- matrixStats::rowMaxs(lad_matrix, na.rm = TRUE)
lad_df$LAD_height_max  <- max.col(lad_matrix, ties.method = "first") * res_z

# Skewness, kurtosis, entropy — still using apply (no native rowSkewness etc.)
lad_df$LAD_skewness    <- apply(lad_matrix, 1, skewness)
gc()
lad_df$LAD_kurtosis    <- apply(lad_matrix, 1, kurtosis)
gc()
lad_df$LAD_entropy     <- apply(lad_matrix, 1, entropy)
gc()

# === ⬛ STAGE: Ecological indices ===
n_layers <- ncol(lad_matrix)
top_third <- seq(ceiling(2 * n_layers / 3), n_layers)
lad_df$Gap_Fraction <- rowMeans(lad_matrix == 0, na.rm = TRUE)
lad_df$Canopy_Cover_Top <- rowMeans(lad_matrix[, top_third] > 0, na.rm = TRUE)
lad_df$LAD_CV <- apply(lad_matrix, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))

# Vertical evenness
lad_sums <- rowSums(lad_matrix, na.rm = TRUE)
lad_df$Vertical_Evenness <- sapply(seq_len(nrow(lad_matrix)), function(i) {
  p <- lad_matrix[i, ] / lad_sums[i]
  p <- p[p > 0]
  -sum(p * log(p)) / log(length(p))
})

# ✅ Remove LAD matrix if done
gc()

# === ⬛ STAGE: Combine data for clustering ===
lad_matrix <- lad_df %>%
  select(starts_with("lad_pulses_"), starts_with("LAD_"), elev, slope, aspect, TPI, CHM) %>%
  as.matrix()
rownames(lad_matrix) <- paste(lad_df$X, lad_df$Y, sep = "_")

# === ⬛ STAGE: Filter for valid species and complete rows ===
valid_idx <- lad_df$species_class %in% valid_ids
lad_df <- lad_df[valid_idx, ]
lad_matrix <- lad_matrix[valid_idx, ]

pulse_cols <- grep("^lad_pulses_", colnames(lad_matrix), value = TRUE)
lad_pulses <- lad_matrix[, pulse_cols]

# Filter valid rows for clustering
valid_rows <- apply(lad_pulses, 1, function(x) all(is.finite(x) & !is.na(x)))
lad_pulses <- lad_pulses[valid_rows, ]
lad_df <- lad_df[valid_rows, ]

# === ⬛ STAGE: PCA sampling + NbClust ===
set.seed(42)
sample_idx <- sample(seq_len(nrow(lad_pulses)), floor(nrow(lad_pulses) * 0.01))
sample_data <- lad_pulses[sample_idx, ]

cat("Sample size before NA filtering:", nrow(sample_data), "\n")
keep_cols <- which(colMeans(is.na(sample_data)) <= 0.2)
sample_data <- sample_data[, keep_cols]
sample_data <- sample_data[, apply(sample_data, 2, var, na.rm = TRUE) > 1e-10]
sample_data <- na.omit(sample_data)
cat("Sample size after cleaning:", nrow(sample_data), "\n")

# ⚠ MEMORY: PCA and NbClust can explode RAM use
pca_res <- prcomp(sample_data, scale. = TRUE)
pc_info <- suggest_n_pcs(pca_res, variance_cutoff = 0.8)
sample_data_pca <- pca_res$x[, 1:pc_info$n_pcs]

nb <- NbClust(sample_data_pca, distance = "euclidean", min.nc = 2, max.nc = 30, method = "kmeans")
optimal_k <- as.integer(nb$Best.nc[1])

        # --- Bereinige ungültige Zeilen für Clustering ---
        valid_rows <- apply(lad_pulses, 1, function(x) all(is.finite(x) & !is.na(x)))
        lad_pulses <- lad_pulses[valid_rows, ]
        lad_df <- lad_df[valid_rows, ]
        
        # --- Clustering der LAD-Pulse-Daten ---
        km_arma  <- ClusterR::KMeans_arma(
          lad_pulses,
          clusters = optimal_k,
          n_iter = 100,
          seed_mode = "random_subset"
        )
     # --- Cluster-Ergebnisse in lad_df schreiben ---
        lad_df$cluster <- as.integer(predict_KMeans(lad_pulses, km_arma))
        
        #lad_df$cluster
        
        # --- Mittelung der LAD-Profile pro Cluster zu synthetischen Profilen ---
        cluster_profiles <- lad_df %>%
          group_by(cluster) %>%
          summarise(across(starts_with("lad_"), mean, na.rm = TRUE)) %>%
          arrange(cluster)
        
        cluster_profiles_lad <- convert_to_LAD_beer(
          df = cluster_profiles,
          grainsize = res_z,     # z. B. 2 m
          k = 0.5,               # Extinktionskoeffizient (typisch 0.3–0.5)
          scale_factor = 1.2,    # optional, empirisch
          lad_max = 3.0,         # realistische Obergrenze für LAD (z. B. 3 m²/m³)
          lad_min = 0.05,        # untere Schranke, optional
          keep_pulses = FALSE    # Originaldaten entfernen
        )
        
        
        # --- Long-Format für LAD-Profile ---
        # --- Long-Format für LAD-Profile ---
        lad_profiles_long <- cluster_profiles_lad %>%
          pivot_longer(
            cols = starts_with("lad_pulses_"),
            names_to = "layer",
            values_to = "lad"
          ) %>%
          mutate(
            z = as.integer(gsub("lad_pulses_|_.*", "", layer)) * res_z,
            ENVIMET_ID = cluster
          ) %>%
          select(ENVIMET_ID, cluster, z, lad)
        
        cluster_ids <- sort(unique(lad_profiles_long$cluster))
        cluster_mapping <- data.frame(
          cluster = cluster_ids,
          ENVIMET_ID = sapply(cluster_ids, int_to_base36)
        )
        
        # --- Mapping-Tabellen vorbereiten ---
        species_mapping <- tibble::tibble(
          species_class = c(1, 2, 3, 4, 5, 6, 7),
          species_name = c("Acer pseudoplatanus", "Betula pendula", "Fagus sylvatica",
                           "Picea abies", "Pinus sylvestris", "Quercus robur", "Tilia cordata")
        )
        
        leaf_thickness_lookup <- tibble::tibble(
          species_name = c("Fagus sylvatica", "Quercus robur", "Acer pseudoplatanus",
                           "Pinus sylvestris", "Picea abies", "Betula pendula", "Tilia cordata"),
          LeafThickness = c(0.0025, 0.0025, 0.0022, 0.0015, 0.0016, 0.0021, 0.0023)
        )
        
        # --- Cluster-IDs generieren (z. B. S00001 ...) ---
        cluster_mapping <- lad_profiles_long %>%
          distinct(cluster) %>%
          mutate(ENVIMET_ID = paste0("S", formatC(cluster, width = 5, flag = "0")))
        
        # --- Voting-basierte Artzuweisung (häufigste Art pro Cluster) ---
        cluster_species <- lad_df %>%
          group_by(cluster, species_class) %>%
          summarise(n = n(), .groups = "drop") %>%
          group_by(cluster) %>%
          slice_max(n, with_ties = FALSE) %>%
          ungroup()
        
        # --- lad_profiles_long neu aufbauen und sauber joinen ---
        lad_profiles_long <- lad_profiles_long %>%
          select(cluster, z, lad) %>%
          left_join(cluster_species, by = "cluster") %>%
          left_join(species_mapping, by = "species_class") %>%
          left_join(leaf_thickness_lookup, by = "species_name") %>%
          left_join(cluster_mapping, by = "cluster")
        
        # --- Optional: Validierung ---
        stopifnot(!any(is.na(lad_profiles_long$species_class)))
        stopifnot(!any(is.na(lad_profiles_long$species_name)))
        stopifnot(!any(is.na(lad_profiles_long$LeafThickness)))
        stopifnot(!any(is.na(lad_profiles_long$ENVIMET_ID)))
        
        cluster_heights <- lad_df %>%
          group_by(cluster) %>%
          summarise(Height_CHM = mean(CHM, na.rm = TRUE), .groups = "drop")
        lad_profiles_long <- lad_profiles_long %>%
          left_join(cluster_heights, by = "cluster")
        
        # --- Traits berechnen ---
        plant_traits <- compute_traits_from_lad(
          lad_df = lad_profiles_long,
          res_z = res_z
        )
        
        
        plant_traits <-plant_traits %>%
          filter(!is.na(LAI), !is.na(MaxLAD))
        
        valid_ids <- plant_traits$ENVIMET_ID
        lad_profiles_long <- lad_profiles_long %>%
          filter(ENVIMET_ID %in% valid_ids)
        
        # --- ENVI-met 3DPLANT-Datei exportieren ---
        export_lad_to_envimet_p3d(
          lad_df = lad_profiles_long,
          file_out = "data/envimet/envimet_pseudo3Dtree.pld",
          res_z = res_z,
          trait_df = plant_traits
        )
        
        # --- Export GPKG with positions and ENVIMET_IDs ---
        lad_df <- lad_df %>% left_join(cluster_mapping, by = "cluster")
        point_df <- lad_df[!duplicated(paste0(lad_df$X, "_", lad_df$Y)), c("X", "Y", "ENVIMET_ID", "species_class", "cluster")]
        sf_points <- st_as_sf(point_df, coords = c("X", "Y"), crs = crs_code)
        st_write(sf_points, output_gpkg, delete_layer = TRUE)
        
        # --- Plot spatial distribution of clusters (optional diagnostics) ---
        if (visualize) {
          library(ggplot2)
          ggplot(lad_df, aes(x = X, y = Y, color = as.factor(cluster))) +
            geom_point(size = 0.8) +
            coord_equal() +
            scale_color_viridis_d(name = "Cluster") +
            labs(title = "Spatial distribution of LAD clusters") +
            theme_minimal()
        }
