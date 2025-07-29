        #' --- ENVI-met 3DPLANT column generator from voxelized ALS data ---
        #'
        #' @author Chris Reudenbach
        #'
        #' This script performs the full pipeline from voxelization of ALS data to the export
        #' of ENVI-met compatible 3DPLANT profiles. It includes LAD computation, clustering
        #' of similar profiles, and export of:
        #' - Point geometries with shared ENVIMET_IDs
        #' - XML-based 3D plant profile database (.pld)
        
        # Load required libraries
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
        library(NbClust)  # newly added for automatic cluster estimation
        
        # Parameters
        # Lookup-Tabelle der Klassen
        ts <- data.frame(
          ID = 1:12,
          value = c(
            "agriculture", "alder", "ash", "beech", "douglas_fir",
            "larch", "oak", "pastures", "roads", "settlements",
            "spruce", "water"
          )
        )
        # Definiere gültige Baumarten
        valid_species <- c("alder", "ash", "beech", "douglas_fir", "larch", "oak", "spruce")
        
        # Erstelle Vektor gültiger Klassen-IDs
        valid_ids <- ts$ID[ts$value %in% valid_species]
        
        visualize=FALSE
        las_file <- here("data/ALS/tiles/")
        res_xy <- 2
        res_z  <- 2
        k      <- 0.3
        scale_factor <- 1.2
        crs_code <- 25832
        output_gpkg <- "data/envimet/envimet_p3dtree_points.gpkg"
        xml_output_file <- "data/envimet/als_envimet_trees.pld"
        # Tree species raster (classified&cleaned)
        # source("src/treespecies.R")
        species_raster <- rast(treespecies_fn) # Tree species raster (classified)
        
        
        # Create output directory if necessary
        dir.create("data/output", showWarnings = FALSE, recursive = TRUE)
        
        source("src/new_utils.R")
       
        ############################################
        
        #--------------------

        # Read and normalize LAS
        las_fn=merge_las_tiles(
          tile_dir = las_file,
          output_file = "data/ALS/merged_output.laz",
          chunk_size = 10000,
          workers = 6
        )
        las=lidR::readLAS(las_fn)
        crs(las) <- "EPSG:25832"
        
        # Crop LAS file to extent (optional xhunking for github)
        #las_cropped <- clip_rectangle(las,sapflow_ext@xmin, sapflow_ext@xmax, sapflow_ext@ymin, sapflow_ext@ymax)
        # Write to file
        #writeLAS(las_cropped, "data/ALS/output_cropped.laz")
        #las <- readLAS("data/ALS/output_cropped.laz")
        
        
        #------------------------------------------------------------------------------
        # Description: Applies CSF ground classification with adaptive parameters
        #              based on CHM rugosity from a single LAS file
        #------------------------------------------------------------------------------
        # === Parameters ===
       ## output_path <- "data/ALS/csf.las"   # Pfad zur Ausgabedatei
        
        # === Step 1: Preliminary CHM for Rugosity Estimation ===
        # (non-normalized, quick canopy height approximation)
        chm_pre <- rasterize_canopy(las, res = res_xy, algorithm = pitfree(thresholds = c(0, 1, 3, 6, 9, 12, 16)))
        
        # === Step 2: Estimate Rugosity to Guide CSF Parameters ===
        rugosity <- terra::focal(chm_pre, w = matrix(1, 3, 3), fun = sd, na.rm = TRUE)
        mean_rug <- global(rugosity, fun = "mean", na.rm = TRUE)[[1]]
        message(sprintf("Mean CHM rugosity: %.3f", mean_rug))
        
        # === Step 3: Choose Adaptive CSF Parameters ===
        csf_params <- if (mean_rug > 1) {
          message("Detected complex/dense canopy – using fine CSF settings")
          csf(cloth_resolution = 0.5, rigidness = 2, class_threshold = 0.4, iterations = 800)
        } else {
          message("Detected open canopy – using coarse CSF settings")
          csf(cloth_resolution = 1.5, rigidness = 4, class_threshold = 0.6, iterations = 300)
        }
        
        # === Step 4: Classify Ground Using Selected CSF ===
        las <- classify_ground(las, csf_params)
        
        # === Step 5: Create DEM from Classified Ground Points ===
        dem_algo <- eval(parse(text = recommend_dem_interpolation(las, res = res_xy)))
        dem <- rasterize_terrain(las, res = res_xy, algorithm = dem_algo)
        
        # === Step 6: Normalize Heights Using DEM ===
        # This sets Z = canopy height above ground
        las_norm <- normalize_height(las, algorithm = knnidw(k = 6L, p = 2))
        
        # === Step 7: CHM from Normalized Point Cloud ===
        pit_algo <- pitfree(thresholds = c(0, 1, 3, 6, 9, 12, 16))
        chm <- rasterize_canopy(las_norm, res = res_xy, algorithm = pit_algo, pkg = "terra")
        
        # === Step 8: DSM from First Returns (Optional) ===
        dsm <- rasterize_canopy(las, res = res_xy, algorithm = p2r(), pkg = "terra")
        
        # === Step 9: Terrain Derivatives from DEM ===
        slope  <- terrain(dem, "slope", unit = "radians")
        aspect <- terrain(dem, "aspect", unit = "degrees")
        
        # === Step 10: TPI from DSM (Smoothed and Binarized) ===
        TPI <- terrain(dsm, "TPI")
        TPI <- terra::focal(TPI, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
        TPI[TPI < 0] <- -1
        TPI[TPI > 0] <- 1
        
        # --- Stack ---
        topo <- c(dem, dsm, chm, slope, aspect, TPI)
        names(topo) <- c("dem", "dsm", "chm", "slope", "aspect", "TPI")
        
        # --- Save ---
        out_file <- "data/ALS/topo_stack.tif"
        terra::writeRaster(topo, out_file, overwrite = TRUE)
        
        # --- Optional Plot ---
        plot(topo)
 
        #############################################
        
        # --- Apply original wide-pulse voxel processing and LAD conversion ---
        voxels <- preprocess_voxels(las_norm, res_xy = res_xy, res_z = res_z)
        lad_df <- convert_to_LAD_beer(voxels, grainsize = res_z, k = k, scale_factor = scale_factor)
        
        # --- Add topographic indices prior to clustering ---
        dem <- rast("data/ALS/topo_stack.tif")[["dem"]]
        dsm <- rast("data/ALS/topo_stack.tif")[["dsm"]]
        
        # Erzeuge ein sf-Objekt aus lad_df
        sf_lad <- st_as_sf(as.data.frame(lad_df), coords = c("X", "Y"), crs = crs_code)
        
        # Füge topografische Daten direkt über sf hinzu
        lad_df$elev  <- terra::extract(dem, sf_lad)[,2]
        lad_df$slope <- terra::extract(terrain(dem, "slope"), sf_lad)[,2]
        lad_df$aspect<- terra::extract(terrain(dem, "aspect"), sf_lad)[,2]
        lad_df$TPI   <- terra::extract(terrain(dsm, "TPI"), sf_lad)[,2]
        lad_df$CHM   <- terra::extract(dsm - dem, sf_lad)[,2]
        lad_df$species_class <- terra::extract(species_raster, sf_lad)[,2]
        
        
        mapviewlad_df = as.data.frame(lad_df)
        # --- Compute structural indices prior to clustering ---
        lad_matrix <- lad_df %>%
          dplyr::select(starts_with("pulses_")) %>%
          as.matrix()
        rownames(lad_matrix) <- paste(lad_df$X, lad_df$Y, sep = "_")
        
        lad_df$LAD_mean       <- rowMeans(lad_matrix, na.rm = TRUE)
        lad_df$LAD_max        <- apply(lad_matrix, 1, max, na.rm = TRUE)
        lad_df$LAD_height_max <- apply(lad_matrix, 1, which.max) * res_z
        lad_df$LAD_skewness   <- apply(lad_matrix, 1, skewness)
        lad_df$LAD_kurtosis   <- apply(lad_matrix, 1, kurtosis)
        lad_df$LAD_entropy    <- apply(lad_matrix, 1, entropy)
        
        
        # --- Additional ecological indices ---
        # Define vertical voxel indices
        n_layers <- ncol(lad_matrix)
        top_third <- seq(ceiling(2 * n_layers / 3), n_layers)
        
        # Gap fraction: Anteil leerer Voxel in Säule
        lad_df$Gap_Fraction <- rowMeans(lad_matrix == 0, na.rm = TRUE)
        
        # Canopy cover in upper third
        lad_df$Canopy_Cover_Top <- rowMeans(lad_matrix[, top_third] > 0, na.rm = TRUE)
        
        # LAD coefficient of variation
        lad_df$LAD_CV <- apply(lad_matrix, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
        
        # Optional: vertical evenness via normalized entropy
        norm_LAD <- lad_matrix / rowSums(lad_matrix, na.rm = TRUE)
        lad_df$Vertical_Evenness <- apply(norm_LAD, 1, function(p) {
          p <- p[p > 0]
          -sum(p * log(p)) / log(length(p))
        })
        # --- Combine LAD profile, structural indices, eco indices and topo features for clustering ---
        lad_matrix <- lad_df %>%
          dplyr::select(starts_with("lad_"), starts_with("LAD_"), Gap_Fraction, Canopy_Cover_Top, LAD_CV, Vertical_Evenness, elev, slope, aspect, TPI, CHM) %>%
          as.matrix()
        
        # --- Combine LAD profile, structural indices, and topo features for clustering ---
        lad_matrix <- lad_df %>%
          dplyr::select(starts_with("pulses_"), starts_with("LAD_"), elev, slope, aspect, TPI, CHM) %>%
          as.matrix()
        rownames(lad_matrix) <- paste(lad_df$X, lad_df$Y, sep = "_")
        
        
        
        
        
        # --- Filter nach gültigen Baumarten-IDs ---
        valid_idx <- lad_df$species_class %in% valid_ids
        lad_df <- lad_df[valid_idx, ]
        lad_matrix <- lad_matrix[valid_idx, ]
        
        # --- Nur Pulse-Spalten extrahieren ---
        pulse_cols <- grep("^pulses_", colnames(lad_matrix), value = TRUE)
        lad_pulses <- lad_matrix[, pulse_cols]
        
        # --- Sampling für Clusterzahlbestimmung ---
        set.seed(42)
        sample_idx <- sample(1:nrow(lad_pulses), floor(nrow(lad_pulses) * 0.01))
        sample_data <- lad_pulses[sample_idx, ]
        
        cat("Sample size before NA filtering:", nrow(sample_data), "rows and", ncol(sample_data), "columns\n")
        
        # Bereinige Sample-Daten
        na_threshold <- 0.2
        keep_cols <- which(colMeans(is.na(sample_data)) <= na_threshold)
        sample_data <- sample_data[, keep_cols]
        
        variances <- apply(sample_data, 2, var, na.rm = TRUE)
        sample_data <- sample_data[, variances > 1e-10]
        
        sample_data <- na.omit(sample_data)
        cat("Sample size after cleaning:", nrow(sample_data), "rows and", ncol(sample_data), "columns\n")
        
        # --- PCA aufbereiten ---
        pca_res <- prcomp(sample_data, scale. = TRUE)
        
        # --- Kriterien auswerten und empfohlene Anzahl PCs bestimmen ---
        pc_info <- suggest_n_pcs(pca_res, variance_cutoff = 0.8)
        
        # --- Reduziertes PCA-Dataset erzeugen ---
        sample_data_pca <- pca_res$x[, 1:pc_info$n_pcs]
        
        # --- Clusteranzahl bestimmen ---
        nb <- NbClust::NbClust(
          sample_data_pca,
          distance = "euclidean",
          min.nc = 2,
          max.nc = 30,
          method = "kmeans"
        )
        optimal_k <- as.integer(nb$Best.nc[1])
          
        # --- Bereinige ungültige Zeilen für Clustering ---
        valid_rows <- apply(lad_pulses, 1, function(x) all(is.finite(x) & !is.na(x)))
        lad_pulses <- lad_pulses[valid_rows, ]
        lad_df <- lad_df[valid_rows, ]
        
        # --- Clustering der LAD-Pulse-Daten ---
        km_arma <- KMeans_arma(
          lad_pulses,
          clusters = optimal_k,
          n_iter = 100,
          seed_mode = "random_subset"
        )
        
        # --- Cluster-Ergebnisse in lad_df schreiben ---
        lad_df$cluster <- as.integer(predict_KMeans(lad_pulses, km_arma))
        
        lad_df$cluster
        
        # --- Mittelung der LAD-Profile pro Cluster zu synthetischen Profilen ---
        cluster_profiles <- lad_df %>%
          group_by(cluster) %>%
          summarise(across(starts_with("pulses_"), mean, na.rm = TRUE)) %>%
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
        point_df <- lad_df[!duplicated(paste0(lad_df$X, "_", lad_df$Y)), c("X", "Y", "ENVIMET_ID.x", "species_class", "cluster")]
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
