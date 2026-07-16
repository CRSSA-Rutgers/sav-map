# ============================================================================= #
# Script Name:  02_build_sat-stack.R
# Description:  Processes PSS imagery from raw tiles into mosaics to create layers
#               to standardize bands, remove boat pixels, & address sun glint.
# Author:       Jess M. Stitt
# Date:         2026-03-31
# ============================================================================= #

# ============================================================================= #
# ---- TASK DESCRIPTION ----
# Input(s):
#  > [TIF] Satellite Image Tiles
#     (PlanetScope Satellite imagery collected for the AOI on a given date)
#  > [SHP] AOI Clip, Full Area
#     (Shapefile of just water pixels for BB-LEH to focus analyses)
#  ??> [SHP] AOI Clip, SAV Habitat Area
#     (Shapefile of just water pixels for BB-LEH to focus analyses)
# Output(s): 
#  > [, ] 
#     ()
# ============================================================================= #
# ---- LIBRARIES SPECIFIC TO THIS TASK ----
library(jsonlite)      #for parsing Planet metadata
# ============================================================================= #
# ---- USER-DEFINED VARIABLES FOR THIS SCRIPT ----

# > Identify outliers: find "Abnormally Bright" pixels using Green & NIR bands
#   **RELEVANT FOR LINES 174-176 BELOW IN SCRIPT
grn_thresh <- 0.6  # Define brightness threshold for green band (visible light)
nir_thresh <- 0.6  # High NIR over water usually means a floating object

# > Prepare for using moving window ('focal') method for boat detection
#   **RELEVANT FOR LINES 177-179 BELOW IN SCRIPT
boat_window <- 3 #moving window size (must be odd number, based on PSS resolution)
boat_fn <- "modal" #function to apply to window: can also try "max" or "mean"

# > Define maximum polygon size to mask (larger areas more likely to be sandbars)
#   **RELEVANT FOR LINES 188-191 BELOW IN SCRIPT
max_boat_area <- 5000

# > Define cutoffs to split boat polygon sizes into small, medium, and large
#   **RELEVANT FOR LINES 192-196 BELOW IN SCRIPT
md_boat_size <- 1500
sm_boat_size <- 500
# ============================================================================= #
# ---- PROCESSING STEPS ----
##------------------------------------------------------------------------------#
## 1. IMAGE GROUPING ----
## >> Consolidate all PSS image tiles taken on same day (& same time frame)
##------------------------------------------------------------------------------#
# > Pull datetime from filename to group by time (images taken <5 min apart)
scene_df <- data.frame(file_path = pss_tiles) %>%
  mutate(time_string = 
           regmatches(basename(file_path), #extract datetime from filename
                      regexpr("\\d{8}_\\d{6}", basename(file_path))),
         acq_time = as.POSIXct(time_string, 
                               format="%Y%m%d_%H%M%S", 
                               tz="UTC")) %>%
  arrange(acq_time) %>%
  mutate(time_diff = as.numeric(difftime(acq_time, 
                                         lag(acq_time, default = acq_time[1]), 
                                         units = "mins")),
         group_id  = cumsum(ifelse(time_diff > 5, 1, 0)))
scene_df <- scene_df[-c(1:3),] #removing duplicate date with bad imagery x3 tiles
# > Calculate time gaps & groups
scene_groups <- scene_df %>%
  group_by(group_id) %>%
  mutate(group_start = min(acq_time), 
         scene_count = n()
  ) %>%
  ungroup()
# > View summary by group
summary_table <- scene_groups %>%
  group_by(group_id, group_start) %>%
  summarize(scenes_in_cluster = n(), .groups = "drop")
message(paste("Satellite Imagery grouped by aquisition date \n > Images, total:",
              nrow(scene_groups),
              "\n > Grouped images (n.datetimes):", 
              nrow(summary_table))); print(summary_table)
# > List the PSS scenes by datetime grouping, to be mosaicked
group_list <- split(scene_df, scene_df$group_id)
# > Subset list to dates of interest (e.g., for image quality or by season)
focal_groups <- group_list[c(8:11,15:18)]
##----------------------------------------------------------------------------#
## 2. IMAGERY PROCESSING PIPELINE ----
## >> 
##----------------------------------------------------------------------------#
message(paste0(">> Beginning Image Processing pipeline: ",
               format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "..."))
inittime <- Sys.time()
# > Iterate through each datetime focal group
for (group in focal_groups) {
  # > Parameterize iteration labels
  curr_id  <- unique(group$group_id) #unique group ID (based on PSS datetime)
  file_timestamp  <- group$time_string[1] #time stamp to use as ID (1st in grp)
  file_ts <- substr(file_timestamp, 1, nchar(file_timestamp) - 2) #only HHmm
  stackdate <- as_date(substr(file_timestamp, 1, nchar(file_timestamp) - 7)) #YMD
  mid_time <- mean(group$acq_time) #average of all PSScene time stamps
  start_d <- format(mid_time - days(1), "%Y-%m-%d") #format dates for retrieval
  end_d   <- format(mid_time + days(1), "%Y-%m-%d") #add 1-day buffer per side
  starttime <- Sys.time()
  message(paste0(">> STARTING: IMAGERY PROCESSING, PSS ", file_ts, " | ", 
                 format(Sys.time(), "%T %Z")))
  ##----------------------------------------------------------------------------#
  ## A. MOSAIC PSS IMAGE TILES ----
  ## >> Join together tiles from same date and time; deal with overlapping areas
  ##----------------------------------------------------------------------------#
  message(paste("   >> Mosaicking PlanetScope satellite scenes...\n      ",
                format(Sys.time(), "%T %Z")))
  mosaic_filepath <- paste0(inrst, "/", file_ts, "_pss8b_00_RMosaicAvg.tif")
  if (file.exists(mosaic_filepath)) {
    # > If it exists, load previously processed mosaic imagery
    print(paste0(">> Existing mosaic found; loading raster stack for ", 
                 stackdate))
    img_focal <- rast(mosaic_filepath)
  } else {
    # > If not, process raw image tiles: load, mosaic, & re-project 8-band imagery
    print(paste0(">> Mosaic not found; initializing raster stitching for ",
                 stackdate))
    img_mosaic <- sprc(lapply(group$file_path, rast))  %>%
      mosaic(fun = "mean") #mosaic PSS tiles; take average of overlapping areas
    img_proj <- project(img_mosaic, proj_crs) #re-project to proper crs
    img_focal <- crop(img_proj, water_poly, mask = TRUE) #mask non-water pixels
    # > Create new mosaicked TIF for raw 8-band imagery (clipped to water extent)
    writeRaster(img_focal,
                file.path(inrst,
                          paste0(file_ts, "_pss8b_00_RMosaicAvg.tif")),
                gdal = c("COMPRESS=LZW", #zips the data to save space
                         "TFW=YES",      #prevents crash at 4GB
                         "BIGTIFF=YES"), #generates sidecar text file
                overwrite = TRUE)
  }
  # > Quantify mosaic area
  valid_mosaic_pixels <- global(img_focal[[1]], "notNA")$notNA
  mosaic_area_m2  <- valid_mosaic_pixels * pixel_area_m2       #square meters
  mosaic_area_km2 <- mosaic_area_m2 / 1000000    #sq.kilometer = 1,000,000 m2
  mosaic_area_ha  <- mosaic_area_m2 / 10000              #hectare = 10,000 m2
  # > Calculate coverage percentage
  pct_coverage <- (mosaic_area_m2 / total_bay_area_m2) * 100
  # > Print results for QAQC
  print("--- MOSAIC COVERAGE EXTENT ---")
  print(paste0("PSS Mosaic Area: ",
               round(mosaic_area_km2, 2), "km² (",
               format(round(mosaic_area_ha, 1), big.mark=","), "ha)"))
  print(paste0("Water Extent covered by Mosaic Footprint: ", 
               round(pct_coverage,0), "%"))
  ##----------------------------------------------------------------------------#
  ## B. DYNAMIC NORMALIZATION ----
  ## >> Scale the values of all bands from 0 to 1 for better comparisons
  ##----------------------------------------------------------------------------#
  message("   >> Dynamically scaling reflectance values from 0 to 1...\n      ",
          format(Sys.time(), "%T %Z"))
  normzd_filepath <- paste0(procrst, "/", file_ts, 
                            "_pss8b_mcf-norm01_normalized.tif")
  if (file.exists(normzd_filepath)) {
    # > If it exists, load previously processed scaled imagery
    print(paste0(">> Existing normalized mosaic found; loading raster stack for ", 
                 stackdate))
    img_norm <- rast(normzd_filepath)
  } else {
  # > Find max. natural brightness: use 99th percentile of Green Band as reference
  ref_val <- global(img_focal[[4]], fun=quantile, probs=0.99, na.rm=TRUE)[1,1]
  # > Normalize stack by ref_val to make 'bright sand' consistently ~0.8-0.9
  img_norm <- img_focal / ref_val
  img_norm <- clamp(img_norm, lower=0, upper=1) #clamp in case outliers >99thPct
  # > Create new mosaicked TIF for 8-band imagery with normalized reflectance
  writeRaster(img_norm, 
              file.path(procrst,
                        paste0(file_ts, "_pss8b_mcf-norm01_normalized.tif")),
              gdal = c("COMPRESS=LZW", #zips the data to save space
                       "TFW=YES",      #prevents crash at 4GB
                       "BIGTIFF=YES"), #generates sidecar text file
              overwrite = TRUE)
  plot(img_norm[[6]], main = paste0("Normalized PSS Mosaic: Green Band (B6) | ",
                                    stackdate))
  }
  ##----------------------------------------------------------------------------#
  ## C. BOAT & WAKE REMOVAL (MASK REFLECTANCE OUTLIERS) ----
  ## >> Watercraft create interference on surface that is highly reflective
  ##----------------------------------------------------------------------------#
  message("   >> Detecting & masking boats and wakes\n      ",
          format(Sys.time(), "%T %Z"))
  boatrm_filepath <- paste0(procrst, "/", file_ts, 
                            "_pss8b_mcf-norm02_boats-rm.tif")
  boatpoly_filepath <- paste0(procmsc, "/boat-outliers_w", boat_window, "_", 
                              file_ts, ".shp")
  if (file.exists(boatrm_filepath)) {
    # > If it exists, load previously processed boat-masked mosaic
    print(paste0(">> Existing boat-masked mosaic found; loading rasters for ", 
                 stackdate))
    img_clean <- rast(boatrm_filepath)
    boat_polys <- sf::st_read(boatpoly_filepath)
  } else {
  # > Create brightness mask for mosaic imagery, using user-defined thresholds
  #   **SEE LINES 027-030 TO SET THESE BRIGHTNESS THRESHOLDS (GREEN & NIR)
  bright_mask <- (img_norm[[4]] > grn_thresh) & (img_norm[[8]] > nir_thresh)
  # > Noise reduction (Erosion/Dilation)
  #   **SEE LINES 032-035 TO SET THE BOAT PARAMETERS (SIZE & FUNCTION)
  clean_mask <- focal(bright_mask, w=boat_window, fun=boat_fn, na.rm=TRUE)
  # > Export boat polygons (only vectorize if outliers were actually found)
  if(global(clean_mask, "sum", na.rm = TRUE)[1,1] > 0) {
    # > Convert bright pixels from raster to polygon
    bright_polys <- as.polygons(clean_mask, dissolve = TRUE) %>% st_as_sf() %>% 
      st_cast("POLYGON")
    st_crs(bright_polys) <- proj_crs
    # > Calculate area of each polygon & add to spatial attribute table
    bright_polys$Area_m2 <- as.numeric(st_area(bright_polys))
    # > Filter out polygons above maximum area (='sandbar' threshold) 
    #   **SEE LINES 037-039 TO SET THE MAXIMUM BOAT POLYGON SIZE
    boat_polys <- bright_polys %>%
      filter(Area_m2 < max_boat_area)
    # > Assign size classes
    #   **SEE LINES 041-044 TO SET THE MEDIUM & SMALL SIZE CLASSES
    boat_polys$Size_Class <- "Large" # Default everything to Large
    boat_polys$Size_Class[boat_polys$Area_m2 < md_boat_size] <- "Medium"
    boat_polys$Size_Class[boat_polys$Area_m2 <= sm_boat_size] <- "Small"
    # > Force as categorical factor
    boat_polys$Size_Class <- factor(boat_polys$Size_Class, 
                                    levels = c("Small", "Medium", "Large"))
    # > See how many boats fall into each category
    table(boat_polys$Size_Class)
    # > Export the polygons as shapefiles
    poly_out_path <- paste0(procmsc, "/boat-outliers_w", boat_window, "_", 
                            file_ts, ".shp")
    sf::st_write(boat_polys, poly_out_path, delete_dsn = TRUE,
                 append = FALSE, quiet = TRUE)
    # > Mask Mosaic CF with boat polygons
    img_clean <- mask(img_norm, boat_polys, inverse = TRUE)
    # > Calculate data loss/retention
    boat_pixels <- freq(clean_mask, value = 1)$count
    boat_area_m2 <- boat_pixels * pixel_area_m2
    boat_area_km2 <- boat_area_m2 / 1000000
    boat_area_ha  <- boat_area_m2 / 10000
    # > Calculate the Boat-Free retention percentage
    pct_lost_boat <- round((boat_pixels / valid_mosaic_pixels) * 100, 2)
    print("--- BOAT/WAKE REMOVAL ---")
    print(table(boat_polys$Size_Class))
    print(paste0("Removed ", nrow(boat_polys), 
                 " bright objects (Boats/Wakes) from imagery"))
    print(paste0("Total Area Masked: ", round(boat_area_km2, 2), "km² (",
                 format(round(boat_area_ha, 1), big.mark=","), "ha)"))
    print(paste0("Total extent lost to boats: ", 
                 round(pct_lost_boat, 2), "%"))
  } else {
    img_clean <- img_norm
    print("--- BOAT/WAKE REMOVAL ---")
    print("!! No Boats/Wakes detected for", paste0(file_ts))
  }
  # > Create new mosaicked TIF for 8-band imagery with boats/wakes masked
  writeRaster(img_clean, 
              file.path(procrst,
                        paste0(file_ts, "_pss8b_mcf-norm02_boats-rm.tif")),
              gdal = c("COMPRESS=LZW", #zips the data to save space
                       "TFW=YES",      #prevents crash at 4GB
                       "BIGTIFF=YES"), #generates sidecar text file
              overwrite = TRUE)
  }
  # gc()
  ##----------------------------------------------------------------------------#
  ## D. SUN GLINT CORRECTION ----
  ## >> Remove sun glint: regress Visible Bands against NIR over Deep Water
  ##----------------------------------------------------------------------------#
  message("   >> Performing sun glint correction...\n      ",
          format(Sys.time(), "%T %Z"))
  deglint_filepath <- paste0(procrst, "/", file_ts, 
                            "_pss8b_mcf-norm03_deglint.tif")
  if (file.exists(deglint_filepath)) {
    # > If it exists, load previously processed deglinted mosaic
    print(paste0(">> Existing deglinted mosaic found; loading raster stack for ", 
                 stackdate))
    img_deglint <- rast(deglint_filepath)
  } else {
  # > Find "Deep Water" pixels (= darkest; lowest 20% of NIR) to calibrate slope
  nir_band <- img_clean[[8]] #Near-Infrared Band (NIR)
  deep_stat <- global(nir_band, fun=quantile, probs=0.2, na.rm=TRUE)
  deep_threshold <- deep_stat[1,1]
  # > Calculate minimum NIR value to use as the "Glint-Free" reference
  min_nir  <- global(nir_band, "min", na.rm=TRUE)[1,1]
  # > Collect sample points from deep water & run regression to calculate slope
  sample_size <- 10000
  s_df <- spatSample(img_clean, size=sample_size, 
                     na.rm=TRUE, method="random") %>% as.data.frame()
  # > Filter only pixels where NIR < Threshold (Deep channels) to calibrate slope
  nir_col_name <- names(nir_band) # e.g., "B8"
  deep_water_samples <- s_df[s_df[[nir_col_name]] < deep_threshold, ]
  # > SAFETY CHECK: Ensure there are enough points
  if(nrow(deep_water_samples) < 100) {
    warning("Not enough pixels found for glint calibration; skipping correction")
    img_clean <- img_clean
  } else {
    print("--- SUN GLINT CORRECTION ---")
    print(paste0("Minimum reference value (NIR, Band 8): ", round(min_nir,4)))
    # > Correct visible bands: iterate through bands 1-6
    img_deglint <- img_clean #start with a clone (for structure)
    for(b_idx in 1:6) {
      vis_col <- names(img_deglint)[b_idx]
      nir_col <- nir_col_name
      # > Calculate slope (regression): Vis ~ NIR
      fit <- lm(as.formula(paste(vis_col, "~", nir_col)), 
                data = deep_water_samples)
      slope <- coef(fit)[2]
      if(slope < 0) slope <- 0
      # > Apply correction by band: New.val = Prior.val - (Slope*(NIR.val-minNIR))
      glint_correction <- slope * (nir_band - min_nir)
      corrected_band <- img_deglint[[b_idx]] - glint_correction
      # > Clamp to 0 to ensure no negative reflectance
      img_deglint[[b_idx]] <- clamp(corrected_band, lower=0.0001)
      print(paste0("Correction Coefficient for B", b_idx ,": ", 
                   round(slope, 4)))
    }
  }
  # > Create TIF for glint-corrected 8-band imagery
  writeRaster(img_deglint,  
              file.path(procrst,
                        paste0(file_ts, "_pss8b_mcf-norm03_deglint.tif")),
              gdal = c("COMPRESS=LZW", #zips the data to save space
                       "TFW=YES",      #prevents crash at 4GB
                       "BIGTIFF=YES"), #generates sidecar text file
              overwrite = TRUE)
  }
  ##----------------------------------------------------------------------------#
  ## 4. DIAGNOSTICS GRAPHICS: Generate Plots to Summarize Process Visually ----
  ##----------------------------------------------------------------------------#
  # CALCULATE AVERAGE BRIGHTNESS: Calculate mean reflectance for every band
  # > Normalized Reflectance (Initial PSS Mosaic)
  band_avgs_rnorm <- global(img_norm, fun = "mean", na.rm = TRUE)
  avg_reflect_rnorm <- as.vector(unlist(band_avgs_rnorm))
  names(avg_reflect_rnorm) <- names(img_norm)
  # > Boat & Wake Removal (Small Bright Outliers masked)
  band_avgs_noboat <- global(img_clean, fun = "mean", na.rm = TRUE)
  avg_reflect_noboat <- as.vector(unlist(band_avgs_noboat))
  names(avg_reflect_noboat) <- names(img_clean)
  # > Sun Glint Correction (Remaining Bright Pixels corrected)
  band_avgs_deglint <- global(img_deglint, fun = "mean", na.rm = TRUE)
  avg_reflect_deglint <- as.vector(unlist(band_avgs_deglint))
  names(avg_reflect_deglint) <- names(img_deglint)
  # > Summarize boat mask stats
  boat_stats <- as.data.frame(boat_polys)
  # > Select a band to show across plots (b1-8; only 1-6 were de-glinted)
  iband <- 4 
  print(paste0("--- AVERAGE REFLECTANCE VALUES (BAND ", iband, ") ---"))
  print(paste0("Normalized Mosaic: ", round(avg_reflect_rnorm[iband], 3)))
  print(paste0("Boat-Free Mosaic:  ", round(avg_reflect_noboat[iband], 3)))
  print(paste0("Deglinted Mosaic:  ", round(avg_reflect_deglint[iband], 3)))
  # DIAGNOSTIC PLOTS: Compare Normalized Original vs Corrected Imagery (band 4)
  p_normalized <- ggplot() +
    geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
    geom_spatraster(data = img_norm[[iband]],  
                    maxcell = 5e5, alpha = 0.8) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    # scale_fill_hypso_c(palette = "spain", direction = 1,
    scale_fill_whitebox_c(palette = "deep", direction = -1,
                          # scale_fill_viridis_c(option = "turbo", direction = 1,
                          guide = "none",
                          na.value = "transparent",
                          limits = c(0,1)) +
    labs(title = "Normalized Mosaic", 
         subtitle = paste0("Avg.Reflectance: ", 
                           round(avg_reflect_rnorm[iband], 3)),
         caption = "Normalized via 99th percentile of Green Band (b4)") +
    theme_void()
  p_boatmask <- ggplot() +
    geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
    geom_spatraster(data = img_clean[[iband]],
                    maxcell = 5e5, alpha = 0.8) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    # scale_fill_hypso_c(palette = "spain", direction = 1,
    scale_fill_whitebox_c(palette = "deep", direction = -1,
                          # scale_fill_viridis_c(option = "turbo", direction = 1,
                          guide = "none",
                          na.value = "transparent",
                          limits = c(0,1)) +
    geom_sf(data = boat_polys, fill = "orange2", color = "orange2", 
            linewidth = 1) +
    labs(title = "Boat & Wake Removal", 
         subtitle = paste0("Avg.Reflectance: ", 
                           round(avg_reflect_noboat[iband], 3)),
         caption = paste0("Data Loss: ", round(pct_lost_boat, 2),
                          "% area masked by boat removal")) +
    theme_void()
  p_deglint <- ggplot() +
    geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
    geom_spatraster(data = img_deglint[[iband]], 
                    maxcell = 5e5, alpha = 0.8) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    # scale_fill_hypso_c(palette = "c3t1", direction = 1,
    scale_fill_whitebox_c(palette = "deep", direction = -1,
                          # scale_fill_viridis_c(option = "turbo", direction = 1,
                          # guide = "none",
                          na.value = "transparent",
                          name = "Normalized \nReflectance",
                          limits = c(0,1)) +
    labs(title = "Sun Glint Correction", 
         subtitle = paste0("Avg.Reflectance: ", 
                           round(avg_reflect_deglint[iband], 3)), 
         caption = paste0("Calibrated using NIR band; slope = ", 
                          round(slope, 2))) +
    theme_void()
  # > Collate: Join plots and include relevant statistics
  diagnostics_imgcorr <- p_normalized + p_boatmask + p_deglint +
    plot_layout(guides = 'collect') +
    plot_annotation(title = paste0(
      "SAV Diagnostics 01 | Imagery Processing Pipeline"), 
      subtitle = paste0(
        "Barnegat Bay Seagrass Mapping, PlanetScope Satellite Imagery for ", 
        stackdate,
        "\n  >> SHOWING: Band ", iband)) +
    theme_void() 
  
  # SUPPLEMENTAL PLOTS
  # > Focus on Boat / Wake Removal by Size Distribution
  p_boats <- ggplot() +
    geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
    geom_spatraster(data = img_clean[[4]], maxcell = 5e5, alpha = 0.6) + 
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    # scale_fill_gradient(low = "gray10", high = "gray90", 
    #                   na.value = "transparent", guide = "none") +
    scale_fill_whitebox_c(palette = "deep", direction = -1,
                          # scale_fill_viridis_c(option = "turbo", direction = 1,
                          guide = "none",
                          na.value = "transparent",
                          limits = c(0,1)) +
    new_scale_fill() +
    geom_spatvector(data = boat_polys, linewidth = 1.2,
                    aes(fill = Size_Class, color = Size_Class)) +
    scale_fill_manual(values = c("Small" = "gold",
                                 "Medium" = "darkorange",
                                 "Large" = "firebrick"), guide = "none") +
    scale_color_manual(values = c("Small" = "gold",
                                  "Medium" = "darkorange",
                                  "Large" = "firebrick"), guide = "none") +
    labs(title = "Boat & Wake Detections by Size", 
         subtitle = "Categorized footprint of transient surface disruption") +
    theme_void()
  # > Plot distribution of polygon size as histogram
  boat_hist <- ggplot(boat_stats, aes(x = Area_m2, fill = Size_Class)) +
    geom_histogram(binwidth = 100, color = "black", alpha = 0.8) +
    geom_vline(xintercept = 500, linetype = "dashed", 
               color = "gray30", linewidth = 1) +
    geom_vline(xintercept = 1500, linetype = "dashed", 
               color = "gray30", linewidth = 1) +
    annotate("text", x = 0, y = Inf, label = " SMALL \n(<500m2)", 
             vjust = 1, hjust = -0.1, size = 3.5, color = "gray30") +
    annotate("text", x = 510, y = Inf, label = " MEDIUM \n(500 - 1500m2)", 
             vjust = 1, hjust = -0.1, size = 3.5, color = "gray30") +
    annotate("text", x = 1510, y = Inf, label = " LARGE \n(>1500m2)", 
             vjust = 1, hjust = -0.1, size = 3.5, color = "gray30") +
    scale_fill_manual(values = c("Small" = "gold", 
                                 "Medium" = "darkorange", 
                                 "Large" = "firebrick")) +
    labs(title = "Distribution of Transient Object Footprints",
         subtitle = "Boat & Wake Frequency by Size",
         x = bquote("Polygon Area ("*m^2*")"), #proper math notation for m-squared
         y = "Frequency (Number of Objects)",
         fill = "Size Class") +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 11),
          axis.title = element_text(face = "bold", size = 12)) +
    theme_minimal()
  # > Collate: Join plots and include relevant statistics
  diagnostics_boatmask <- boat_hist + p_boats + 
    plot_layout(guides = 'collect') +
    plot_annotation(title = 
    "SUPP Diagnostics 01  |  Boat & Wake Detection, Size Distribution", 
                    subtitle = paste0(
    "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
                      stackdate),
                    caption = paste0("BASE: PlanetScope Band ", iband, 
    " | *Only polygons <5,000m² included & masked (n=", nrow(boat_polys), ")")) +
    theme_void() 
  # print(diagnostics_boatmask)
  # > EXPORT GRAPHICS: Summary diagnostics
  ggsave(file.path(procimg, 
                   paste0("01_SAV-diagnostics_img-correction-analysis_",   
                          file_ts, ".png")), diagnostics_imgcorr,
         width = 11, height = 9, dpi = 300)
  ggsave(file.path(procimg, 
                   paste0("S01_SUPP-diagnostics_boat-detection-distr_",   
                          file_ts, ".png")), diagnostics_boatmask,
         width = 11, height = 9, dpi = 300)
  # CLEAN UP: Remove temporary objects & trigger garbage collection
  rm(img_norm, img_clean, img_deglint, boat_polys, boat_stats,
     band_avgs_rnorm, band_avgs_noboat, band_avgs_deglint) 
  gc()
  endtime <- Sys.time()
  # duration <- endtime - starttime
  loop_time <- round(as.numeric(difftime(endtime, starttime, units = "mins")), 2)
  message(paste0(">> FINISHED: IMAGERY PROCESSING, PSS ", file_ts, " | ", 
                 format(Sys.time(), "%T %Z"), "\n   Duration:", loop_time,
                 "minutes", "\n----------------------------------------"))
  beep(2) #----
}
elapsed <- round(as.numeric(difftime(endtime, inittime, units = "mins")), 2)
message(paste0(">> Imagery Processing completed for all dates on \n   ", 
               format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "  |  Elapsed Time: ", 
               elapsed, " minutes"))
Sys.sleep(2)
beep(10)

# ---- OUTPUTS FROM SCRIPT ----
# 

