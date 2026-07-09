# =============================================================================
# LANDSCAPE METRICS  (Stage 2 of the pipeline)
# =============================================================================
# Calculates habitat amount and configuration metrics in 150 / 500 / 1000 m
# buffers around each survey point, for each fragmentation definition
# (hyp1 = Definition A, hyp2 = B, hyp3 = C). Produces the three files that feed
# the occupancy analysis (script 01):
#     hyp{1,2,3}_landscape_metrics_2m_res.csv
#
# INPUTS (produced by Stage 1, the ArcGIS Pro notebook):
#   dissolved_habitat_hyp{h}_{year}.shp   suitable-habitat patches, disturbances
#                                         removed, for definition h and year
#   patch_metrics_hyp{h}_{year}.csv        survey points (with coordinates and
#                                         point metadata) for that definition/year
# These inputs are derived from the licensed AVI / Human Footprint layers and are
# NOT redistributed; see the repository README for how to regenerate them.
#
# OUTPUT: outputs/landscape_metrics/2m_resolution/hyp{h}_landscape_metrics_2m_res.csv
#   (copy the three hyp*_landscape_metrics_2m_res.csv files into the occupancy
#    analysis's data/ folder to run script 01)
#
# Habitat patches are rasterised at 2 m so narrow features (e.g. 2 m seismic
# lines) are preserved rather than merged with adjacent habitat.
# =============================================================================

library(sf)
library(landscapemetrics)
library(tidyverse)
library(raster)
library(fasterize)

# ---- Global parameters ------------------------------------------------------
RESOLUTION   <- 2                    # metres (fine enough to keep 2 m features)
BUFFER_SIZES <- c(150, 500, 1000)    # metres (7 ha, 78 ha, 314 ha)

# =============================================================================
# METRICS FOR ONE POINT AT ONE BUFFER SIZE
# =============================================================================
# Clips the habitat polygons to a circular buffer, rasterises them at 2 m, and
# computes the class-level metrics for the habitat class. prop_habitat is taken
# from the polygon areas over the true circular buffer area. Where a buffer
# contains no habitat, configuration metrics that are undefined without habitat
# (aggregation, cohesion, contagion) are returned as NA.

calculate_buffer_metrics <- function(habitat_polygons, point, buffer_radius,
                                     resolution = RESOLUTION) {
  gc(verbose = FALSE)
  point_buffer      <- st_buffer(point, buffer_radius)
  habitat_in_buffer <- suppressWarnings(st_intersection(habitat_polygons, point_buffer))

  buffer_area  <- pi * buffer_radius^2
  habitat_area <- if (nrow(habitat_in_buffer) == 0) 0 else
    as.numeric(sum(st_area(habitat_in_buffer)))
  prop_habitat <- habitat_area / buffer_area

  empty <- list(prop_habitat = prop_habitat, edge_density = 0, aggregation_index = NA,
                cohesion = NA, n_patches = 0, mean_patch_area = 0, contagion = NA)
  if (nrow(habitat_in_buffer) == 0) return(empty)

  # Rasterise habitat = 1, everything else = 0, then clip to the circular buffer.
  habitat_in_buffer$value <- 1
  bbox <- st_bbox(point_buffer)
  r <- raster(xmn = bbox["xmin"], xmx = bbox["xmax"],
              ymn = bbox["ymin"], ymx = bbox["ymax"],
              res = resolution, crs = st_crs(habitat_polygons))
  r <- fasterize(habitat_in_buffer, r, field = "value", background = 0)
  r <- mask(r, point_buffer)                       # outside the circle -> NA
  if (cellStats(r, "sum") == 0) return(empty)

  # Pull the value for the HABITAT class (class == 1) from a class-level metric.
  c_val <- function(f, cls = 1) {
    out <- tryCatch(f(r), error = function(e) NULL)
    if (is.null(out) || !"class" %in% names(out)) return(NA_real_)
    v <- out$value[out$class == cls]
    if (length(v) == 0) NA_real_ else v[1]
  }

  ed <- c_val(lsm_c_ed)
  list(
    prop_habitat      = prop_habitat,
    edge_density      = if (is.na(ed)) 0 else ed,   # one class => no edge => 0
    aggregation_index = c_val(lsm_c_ai),
    cohesion          = c_val(lsm_c_cohesion),
    n_patches         = c_val(lsm_c_np),
    mean_patch_area   = c_val(lsm_c_area_mn),
    contagion         = tryCatch(lsm_l_contag(r)$value, error = function(e) NA_real_)
  )
}

# =============================================================================
# PROCESS ALL POINTS x BUFFERS FOR ONE DEFINITION (hypothesis)
# =============================================================================
# Loops over years, reads that year's habitat patches and survey points, and
# computes metrics at every buffer size for every point. Writes one combined CSV
# per definition plus per-year and interim CSVs for safety on long runs.

process_landscape_metrics <- function(hyp_num, input_dir, output_dir,
                                      test_year = NULL, buffer_sizes = BUFFER_SIZES) {

  res_folder <- file.path(output_dir, paste0(RESOLUTION, "m_resolution"))
  dir.create(res_folder, recursive = TRUE, showWarnings = FALSE)

  habitat_pattern <- paste0("^dissolved_habitat_hyp", hyp_num, "_\\d{4}\\.shp$")
  habitat_files   <- list.files(input_dir, pattern = habitat_pattern, full.names = TRUE)

  if (!is.null(test_year)) {
    habitat_files <- habitat_files[grep(test_year, habitat_files)]
    if (length(habitat_files) == 0)
      stop(sprintf("No habitat files for hypothesis %d and year %s", hyp_num, test_year))
  }

  all_results <- data.frame()

  for (habitat_file in habitat_files) {
    gc(verbose = FALSE)
    year <- as.numeric(stringr::str_extract(basename(habitat_file), "\\d{4}"))
    if (is.na(year)) next

    points_csv <- file.path(input_dir, paste0("patch_metrics_hyp", hyp_num, "_", year, ".csv"))
    if (!file.exists(points_csv)) {
      cat(sprintf("No points file for hypothesis %d, year %d\n", hyp_num, year)); next
    }

    cat(sprintf("\nProcessing hypothesis %d, year %d\n", hyp_num, year))
    habitat     <- st_read(habitat_file, quiet = TRUE)
    points_data <- read.csv(points_csv)

    # Find the coordinate columns (naming varies by export).
    if (all(c("Easting", "Northing") %in% names(points_data))) {
      coords <- c("Easting", "Northing")
    } else if (all(c("POINT_X", "POINT_Y") %in% names(points_data))) {
      coords <- c("POINT_X", "POINT_Y")
    } else if (all(c("X", "Y") %in% names(points_data))) {
      coords <- c("X", "Y")
    } else {
      cat("  Could not find coordinate fields, skipping year", year, "\n"); next
    }

    points_sf <- st_as_sf(points_data, coords = coords, crs = st_crs(habitat))
    cat(sprintf("  Processing %d points...\n", nrow(points_sf)))

    for (i in 1:nrow(points_sf)) {
      if (i %% 25 == 0 || i == nrow(points_sf))
        cat(sprintf("  point %d of %d\n", i, nrow(points_sf)))

      point    <- points_sf[i, ]
      min_dist <- st_distance(point, habitat) %>% min()
      attr(min_dist, "units") <- NULL

      for (size in buffer_sizes) {
        metrics <- calculate_buffer_metrics(habitat, point, size)

        # Point metadata (field names vary between exports).
        surveyid <- if ("surveyid" %in% names(point)) as.character(point$surveyid) else
                    if ("site_id"  %in% names(point)) as.character(point$site_id) else
                    paste0("point_", i)
        survey_year <- if ("survey_yea"  %in% names(point)) as.numeric(point$survey_yea) else
                       if ("survey_year" %in% names(point)) as.numeric(point$survey_year) else year
        veg_type    <- if ("Veg_Type"   %in% names(point)) as.character(point$Veg_Type) else NA
        origin_year <- if ("Origin_Yea" %in% names(point)) as.numeric(point$Origin_Yea) else NA

        all_results <- rbind(all_results, data.frame(
          year = year, surveyid = surveyid, survey_year = survey_year,
          veg_type = veg_type, origin_year = origin_year,
          buffer_size = size, resolution = RESOLUTION, dist_to_habitat = min_dist,
          as.data.frame(t(unlist(metrics)))
        ))
      }

      # Interim save every 50 points (long runs).
      if (i %% 50 == 0)
        write.csv(all_results,
                  file.path(res_folder, paste0("hyp", hyp_num,
                            "_landscape_metrics_", RESOLUTION, "m_res_interim.csv")),
                  row.names = FALSE)
      if (i %% 10 == 0) gc(verbose = FALSE)
    }

    # Per-year save.
    write.csv(all_results[all_results$year == year, ],
              file.path(res_folder, paste0("hyp", hyp_num, "_year", year,
                        "_landscape_metrics_", RESOLUTION, "m_res.csv")),
              row.names = FALSE)
  }

  # Combined save (this is the file the occupancy analysis reads).
  output_file <- file.path(res_folder,
                  paste0("hyp", hyp_num, "_landscape_metrics_", RESOLUTION, "m_res.csv"))
  write.csv(all_results, output_file, row.names = FALSE)
  cat(sprintf("\nHypothesis %d done -> %s\n", hyp_num, output_file))

  all_results
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================
# Set these to your local paths. input_dir holds the Stage 1 outputs
# (dissolved_habitat_*.shp and patch_metrics_*.csv). 

input_dir  <- "gis_outputs"                      
output_dir <- "outputs/landscape_metrics"      

cat("Resolution:", RESOLUTION, "m | Buffers:",
    paste(BUFFER_SIZES, collapse = ", "), "m\n")

# First test just one year for one definition before launching the full run
# process_landscape_metrics(1, input_dir, output_dir, test_year = "2020")

# Full run: all three fragmentation definitions (hyp1 = A, hyp2 = B, hyp3 = C).
for (hyp in 1:3) {
  cat(sprintf("\n===== Hypothesis %d =====\n", hyp))
  process_landscape_metrics(hyp_num = hyp, input_dir = input_dir, output_dir = output_dir)
  gc(verbose = FALSE)
}

# Now copy the three hyp*_landscape_metrics_2m_res.csv files into "2_occupancy_analysis/data/" to run script 01
