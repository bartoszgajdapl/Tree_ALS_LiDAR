#### Environment setup and required libraries ####
library(lidR)
library(sf)
library(terra)
library(dplyr)
library(mapview)
library(fs)
library(stringr)

#### Load LiDAR data catalog ####
ctg <- readLAScatalog("pointclouds_t/") # Load directory with LiDAR data
opt_chunk_size(ctg) <- 0 # Set chunk (tile) size to 0 for auto (let lidR decide)
opt_chunk_buffer(ctg) <- 20 # Set 20m buffer for better stitching of adjacent tiles
opt_filter(ctg) <- "-keep_class 2 4 5" # Filter: keep ground (2), medium (4) and high vegetation (5) classes
opt_output_files(ctg) <- "output/chm/{*}_CHM" # Set output path template for CHM rasters
opt_laz_compression(ctg) <- TRUE # Save in LAZ format to save space
plot(ctg) # Quick preview of available tiles

# Initialize dataframe for tracking original vs. new tile IDs
tile_id_map <- data.frame(old = character(), new = character(), stringsAsFactors = FALSE)

#### Prepare the "output" folder ####

# Check if the "output" folder exists; 
# If not – create it. If it exists:
#     → check if it is empty. 
#     → if not empty, ask in console whether to clean it up
# Display a message with user's choice and the consequences

if (dir_exists("output")) {
  files_inside <- dir_ls("output", recurse = TRUE)
  
  if (length(files_inside) > 0) {
    cat("Folder 'output/' contains data. Do you want to clean it? (y/n): ")
    answer <- tolower(readline())
    
    if (answer == "y") {
      cat("Cleaning 'output/' directory...\n")
      dir_delete("output")
      dir_create("output")
      dir_create("output/centroids")
      dir_create("output/crowns")
      dir_create("output/chm")
      dir_create("output/pointcloud_segmented")
      cat("Folder 'output/' cleaned and ready.\n")
    } else {
      cat("Folder was NOT cleaned. Results may be overwritten or duplicated!\n")
    }
  } else {
    cat("Folder 'output/' exists and is empty — OK.\n")
  }
} else {
  cat("Creating a new 'output/' folder and subfolders...\n")
  dir_create("output")
  dir_create("output/centroids")
  dir_create("output/crowns")
  dir_create("output/chm")
  dir_create("output/pointcloud_segmented")
  cat("Folders created!\n")
}

#### Define function for processing LiDAR tiles in the catalog ####

process_tile <- function(cluster, ...) {
  las <- readLAS(cluster)
  if (is.null(las) || npoints(las) == 0) return(NULL)
  
  # Get the filename from the cluster
  tile_filename <- basename(cluster@files[1])
  
  tile_counter <<- tile_counter + 1
  tile_id_final <- sprintf("tile_%06d", tile_counter)
  tile_id_map <<- rbind(tile_id_map, data.frame(file = tile_filename, tile_id = tile_id_final, stringsAsFactors = FALSE))
  cat("\nProcessing:", tile_filename, "as", tile_id_final, "\n")
  
  # Height normalization
  las_norm <- normalize_height(las, algorithm = knnidw(), add_lasattribute = TRUE)
  
  # CHM generation
  CHM <- rasterize_canopy(
    las_norm, # point cloud after normalization
    res = 1, # CHM raster resolution (recommended: 0.75–1)
    algorithm = pitfree( 
      thresholds = c(0, 2, 5, 10, 15, 20, 30),
      max_edge = c(0, 4),
      subcircle = 0.3,
      highest = TRUE
    )
  )
  if (all(is.na(values(CHM)))) return(NULL)
  
  # Smooth CHM
  kernel <- matrix(1, 5, 5)
  CHM_smooth <- terra::focal(CHM, w = kernel, fun = mean, na.rm = TRUE)
  
  # Tree top detection and segmentation
  vwf_fun <- function(x) { 0.15 * x + 5 } # Variable window function
  ttops <- locate_trees(CHM_smooth, lmf(vwf_fun, hmin = 6)) # Detect local maxima (tree tops) higher than 6m
  las_segmented <- segment_trees(las_norm, dalponte2016(CHM_smooth, ttops)) # Segment the point cloud
  
  # Assign global treeID
  tree_ids <- sort(unique(na.omit(las_segmented@data$treeID)))
  n_trees <- length(tree_ids)
  tile_num <- as.numeric(str_extract(tile_id_final, "\\d+"))
  tree_prefix <- tile_num * 1e6
  new_ids <- tree_prefix + seq_len(n_trees)
  
  mapping <- data.frame(treeID = tree_ids, treeID_global = new_ids)
  
  las_segmented@data <- las_segmented@data %>%
    left_join(mapping, by = "treeID") %>%
    mutate(treeID = treeID_global) %>%
    select(-treeID_global)
  
  # Crown outlines (delineate_crowns)
  crowns_sf <- suppressWarnings(delineate_crowns(las_segmented, func = .stdtreemetrics, type = "concave", concavity = 2))
  
  # Convert to sf and filter
  crowns_sf <- st_as_sf(crowns_sf)
  crowns_sf <- crowns_sf %>% filter(Z >= 6) %>% 
    select(treeID, Z, npoints, convhull_area, geometry)
  
  # Centroids
  centroids_geom <- st_centroid(st_geometry(crowns_sf))
  centroids_sf <- st_sf(
    treeID = crowns_sf$treeID,
    Z = crowns_sf$Z,
    geometry = centroids_geom,
    crs = st_crs(crowns_sf)
  )
  
  # Tree tops
  ttops_sf <- st_zm(st_as_sf(ttops))
  
  # Save outputs
  st_write(centroids_sf, paste0("output/centroids/", tile_id_final, "_centroids.geojson"),
           delete_dsn = TRUE, quiet = TRUE)
  st_write(crowns_sf, paste0("output/crowns/", tile_id_final, "_crowns.geojson"),
           delete_dsn = TRUE, quiet = TRUE)
  
  # Save .laz file with treeIDs
  writeLAS(las_segmented, paste0("output/pointcloud_segmented/", tile_id_final, ".laz"))
  
  cat("Saved:", tile_id_final, "\n")
  return(CHM_smooth)
}

#### Run processing function for files in the catalog ####

# The process_tile function will be called for each LAS/LAZ file
tile_counter <- 0  # Reset before starting

# Run parallel processing for all LiDAR tiles in the catalog
catalog_apply(ctg, function(lascluster, ...) {
  process_tile(lascluster, .filename = attr(lascluster, "filename"))
})

# Save to CSV: mapping of tiles with their random ID + timestamp and new assigned ID
tile_id_map_unique <- tile_id_map[!duplicated(tile_id_map$tile_id), ]
write.csv(tile_id_map_unique, "output/tile_id_mapping.csv", row.names = FALSE)

#### Merge all results into a single file ####

# Merge all centroids and crowns into a single vector layer
# Also, merge all CHM rasters into a single mosaic
centroid_files <- list.files("output/centroids", pattern = "\\.geojson$", full.names = TRUE)
centroids_all <- do.call(rbind, lapply(centroid_files, st_read, quiet = TRUE))
crown_files <- list.files("output/crowns", pattern = "\\.geojson$", full.names = TRUE)
crowns_all <- do.call(rbind, lapply(crown_files, st_read, quiet = TRUE))

# Save merged files
st_write(centroids_all, "output/centroids/centroids_all.geojson", delete_dsn = TRUE)
st_write(crowns_all, "output/crowns/crowns_all.geojson", delete_dsn = TRUE)

# Merge and save the CHM rasters mosaic for the whole area
chm_files <- list.files("output/chm", pattern = "\\.tif$", full.names = TRUE)
chm_rasters <- lapply(chm_files, terra::rast)
CHM_merged <- do.call(terra::mosaic, chm_rasters)
writeRaster(CHM_merged, "output/chm/CHM_merged.tif", overwrite = TRUE)

#### Handle duplicates in buffer zone ####

tile_bounds <- st_as_sf(ctg) # Convert catalog extents to sf object
plot(tile_bounds$geometry) # Preview tile boundaries

# Function to expand tile bbox by buffer size
expand_bbox <- function(poly_sf, buffer) {
  bb <- st_bbox(poly_sf)
  bb_exp <- bb
  bb_exp["xmin"] <- bb["xmin"] - buffer
  bb_exp["ymin"] <- bb["ymin"] - buffer
  bb_exp["xmax"] <- bb["xmax"] + buffer
  bb_exp["ymax"] <- bb["ymax"] + buffer
  st_as_sfc(bb_exp)
}

tile_bounds <- st_as_sf(ctg)  # Tiles as sf
buffer_size <- 20  # Buffer in meters

# Iterate each tile to create a buffered version of its extent
buffered_list <- lapply(1:nrow(tile_bounds), function(i) {
  geom <- expand_bbox(tile_bounds[i, ], buffer_size)
  st_sf(tile_id = basename(tile_bounds$filename[i]), geometry = geom, crs = st_crs(tile_bounds))
})

# Combine list of buffered tiles into one sf object
tile_bounds_buffered <- do.call(rbind, buffered_list)

# Preview tile extents and buffers
mapview::mapview(tile_bounds_buffered, zcol = "tile_id") +
  mapview::mapview(tile_bounds, col.regions = "black")

# Create all combinations of tile pairs
tile_combos <- combn(1:nrow(tile_bounds_buffered), 2)

# Calculate intersections only where buffers actually overlap
overlap_areas <- lapply(1:ncol(tile_combos), function(i) {
  a <- tile_bounds_buffered[tile_combos[1, i], ]
  b <- tile_bounds_buffered[tile_combos[2, i], ]
  
  if (st_intersects(a, b, sparse = FALSE)[1, 1]) {
    st_intersection(a, b)
  } else {
    NULL
  }
})

# Remove empty results and merge to one layer
overlap_areas <- do.call(rbind, overlap_areas[!sapply(overlap_areas, is.null)])

# Save the detected buffer overlap area to file
st_write(overlap_areas, "output/buffer_overlap_zones.geojson", delete_dsn = TRUE)

#### Define function for removing duplicated trees and centroids in buffer overlaps ####

remove_duplicates_from_overlap <- function(crowns_all, centroids_all, overlap_areas, threshold = 0.7) {
  cat("Removing duplicate trees in buffer overlap zones...\n")
  
  # Union all buffer areas into a single polygon
  overlap_union <- suppressWarnings(st_union(overlap_areas))
  
  # Intersect crowns with buffer area
  crowns_in_overlap <- suppressWarnings(st_intersection(crowns_all, overlap_union))
  if (nrow(crowns_in_overlap) == 0) {
    warning("No crowns in buffer area — skipping duplicate removal.")
    return(list(crowns = crowns_all, centroids = centroids_all))
  }
  
  # Assign IDs
  crowns_in_overlap$id <- seq_len(nrow(crowns_in_overlap))
  
  # Check intersections
  intersections <- st_intersects(crowns_in_overlap)
  to_remove <- c()
  pb <- txtProgressBar(min = 1, max = length(intersections), style = 3)
  
  for (i in seq_along(intersections)) {
    overlaps <- intersections[[i]]
    for (j in overlaps) {
      if (i >= j) next
      
      a <- crowns_in_overlap[i, ]
      b <- crowns_in_overlap[j, ]
      inter <- suppressWarnings(st_intersection(a, b))
      
      # Safe check
      if (!is.null(inter) && nrow(inter) > 0) {
        area_inter <- tryCatch(st_area(inter), error = function(e) NA)
        area_a <- tryCatch(st_area(a), error = function(e) NA)
        area_b <- tryCatch(st_area(b), error = function(e) NA)
        
        if (is.na(area_inter) || is.na(area_a) || is.na(area_b)) next
        
        min_area <- min(area_a, area_b)
        ratio <- as.numeric(area_inter / min_area)
        
        if (!is.na(ratio) && ratio > threshold) {
          id_to_drop <- if (area_a < area_b) a$id else b$id
          to_remove <- c(to_remove, id_to_drop)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  crowns_to_drop <- crowns_in_overlap$treeID[crowns_in_overlap$id %in% to_remove]
  crowns_clean <- crowns_all[!crowns_all$treeID %in% crowns_to_drop, ]
  centroids_clean <- centroids_all[!centroids_all$treeID %in% crowns_to_drop, ]
  
  # Export results
  st_write(crowns_clean, "output/crowns/crowns_clean.geojson", delete_dsn = TRUE, quiet = TRUE)
  st_write(centroids_clean, "output/centroids/centroids_clean.geojson", delete_dsn = TRUE, quiet = TRUE)
  st_write(crowns_all[crowns_all$treeID %in% crowns_to_drop, ], "output/crowns/crowns_removed.geojson", delete_dsn = TRUE, quiet = TRUE)
  st_write(centroids_all[centroids_all$treeID %in% crowns_to_drop, ], "output/centroids/centroids_removed.geojson", delete_dsn = TRUE, quiet = TRUE)
  
  cat("\nDuplicate trees in the buffer zone removed and saved!\n")
  cat("Removed", length(crowns_to_drop), "duplicate trees.\n")
  cat("\nScript finished, all files have been processed and saved in the output location!\n")
  return(list(crowns = crowns_clean, centroids = centroids_clean))
}

#### Run the function to remove duplicated trees detected in buffer overlap ####
result <- remove_duplicates_from_overlap(
  crowns_all = crowns_all,
  centroids_all = centroids_all,
  overlap_areas = overlap_areas,
  threshold = 0.7 # Percentage of overlap (here: 70%) that triggers removal
)

#### OPTIONAL – visualization, statistics ####
### Display interactive map ###
# Load generated crowns and crown centroids into variables
crowns_clean <- st_read("output/crowns/crowns_clean.geojson", quiet = TRUE)
centroids_clean <- st_read("output/centroids/centroids_clean.geojson", quiet = TRUE)

# Show interactive map
mapview(crowns_clean, col.regions = "forestgreen", alpha.regions = 0.4, layer.name = "Tree crowns") +
  mapview(centroids_clean, col.regions = "red", label = centroids_clean$treeID, layer.name = "Tree centroids")

### Display point cloud ###
pointcloud_segm <- readLAS("output/pointcloud_segmented/tile_000001.laz", filter = "-keep_class 4 5") # Load a segmented LAZ tile (skip ground)
las_trees <- filter_poi(pointcloud_segm, !is.na(treeID)) # Filter points that have assigned treeID (are trees)
tree_ids <- unique(na.omit(pointcloud_segm@data$treeID)) # Get list of treeIDs in the tile
tree_ids <- sort(tree_ids) # Sort treeID list
tree_ids # Print list of treeIDs

plot(las_trees, color = "treeID") # Show whole tile with detected trees (color = treeID)
plot(filter_poi(pointcloud_segm, treeID == 1000071), color = "Z", legend=TRUE, axis=TRUE) # Show one tree by its treeID

### Additional interesting statistics ###

# Tallest trees (top 3)
top3_height <- crowns_clean %>%
  arrange(desc(Z)) %>%
  slice(1:3) %>%
  select(treeID, Z)
print(top3_height)

# Largest crown area (top 3)
top3_crownarea <- crowns_clean %>%
  arrange(desc(convhull_area)) %>%
  slice(1:3) %>%
  select(treeID, convhull_area, Z)
print(top3_crownarea)