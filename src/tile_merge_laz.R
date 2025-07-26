#' Retile a LAS/LAZ file into regular tiles
#'
#' This function splits a LAS/LAZ file into regular square tiles using `lidR::catalog_retile`.
#'
#' @param input_file Path to input LAS/LAZ file
#' @param out_dir Output directory for tiled files
#' @param chunk_size Size of each tile (in meters), default = 200
#' @param output_ext File extension, either "laz" or "las" (default: "laz")
#' @param buffer Buffer size between tiles (default: 0 for non-overlapping tiles)
#' @param workers Number of parallel workers to use (default: 4)
#' @return Path to resulting `LAScatalog` object
#' @export
#'
#' @examples
#' retile_las("data/input.las", out_dir = "tiles/", chunk_size = 100)

retile_las <- function(input_file,
                       out_dir = "tiles/",
                       chunk_size = 250,
                       output_ext = "laz",
                       buffer = 0,
                       workers = 4) {
  if (!file.exists(input_file)) stop("Input file not found.")
  if (!output_ext %in% c("laz", "las")) stop("Invalid extension: use 'laz' or 'las'.")
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load packages
  library(lidR)
  library(future)
  
  # Enable multithreading
  set_lidr_threads(workers)
  future::plan(multisession, workers = workers)
  
  # Read and prepare catalog
  ctg <- readLAScatalog(input_file)
  opt_laz_compression(ctg) <- TRUE
  opt_chunk_size(ctg)   <- chunk_size
  opt_chunk_buffer(ctg) <- buffer
  opt_output_files(ctg) <- file.path(out_dir, paste0("tile_{XLEFT}_{YBOTTOM}.", output_ext))
  
  # Retile into square chunks
  message("Tiling ", input_file, " into ", out_dir)
  tiled_ctg <- catalog_retile(ctg)
  
  return(tiled_ctg)
}


#' Merge LAS/LAZ tiles into a single file using lidR::catalog_retile
#'
#' This function merges tiled LAS/LAZ files from a LAScatalog by retiling the whole extent
#' into one chunk (no overlap) and writing a single merged output file.
#'
#' @param tile_dir Directory containing LAS/LAZ tiles (can be a single file too)
#' @param output_file Path to output LAS/LAZ file (default: "merged.laz")
#' @param chunk_size Optional chunk size to control memory (default: 10000)
#' @param workers Number of parallel workers (default: 4)
#' @return Path to merged file as character string
#' @export
#'
#' @examples
#' merge_las_tiles("tiles/", "merged.laz")

merge_las_tiles <- function(tile_dir,
                            output_file = "merged.laz",
                            chunk_size = 10000,
                            workers = 4) {
  if (!dir.exists(tile_dir) && !file.exists(tile_dir)) {
    stop("Input tile directory or file does not exist.")
  }
  
  # Load packages
  library(lidR)
  library(future)
  
  # Set parallel options
  set_lidr_threads(workers)
  future::plan(multisession, workers = workers)
  
  # Read catalog (even if just one file)
  ctg <- readLAScatalog(tile_dir)
  
  # Configure catalog for merging
  opt_chunk_size(ctg)   <- chunk_size
  opt_chunk_buffer(ctg) <- 0
  opt_output_files(ctg) <- sub("\\.la[sz]$", "", output_file)
  
  message("Merging tiles into: ", output_file)
  
  merged_ctg <- catalog_retile(ctg)
  
  return(paste0(sub("\\.la[sz]$", "", output_file), ".laz"))
  
  

  
  
  
  retile_las(
    input_file = "data/ALS/output_cropped",
    out_dir = "data/ALS/tiles",
    chunk_size = 260,
    output_ext = "laz",
    buffer = 0,
    workers = 6
  )
  
  merge_las_tiles(
    tile_dir = "data/ALS/tiles",
    output_file = "data/ALS/merged_output.laz",
    chunk_size = 10000,
    workers = 6
  )
}