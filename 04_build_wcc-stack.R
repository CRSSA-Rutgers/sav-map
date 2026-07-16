# ==============================================================================
# Script Name:  04_build_wcc-stack.R
# Description:  Creates water column correction layers based on PSS image stack to 
#               account for light attenuation, turbidity, and algae.
# Author:       Jess M. Stitt
# Date:         2026-03-31
# ==============================================================================

# ============================================================================= #
# ---- TASK DESCRIPTION ----
# Water Column Correction (Water Quality Metrics)
# Input(s):
#  > [TIF] Normalized, Boat-Masked PSS Image Mosaics (**with or w/o De-Glinting) 
#  > [SHP] Sand Polygons, drawn by hand in areas with known sand (no SAV) 
# Output(s): 
#  > [TIF] Raster Stack of Water Column Corrected Imagery: Sagawa x8, Lyzenga X1 
#  ** ADDITIONAL CODE FOR TURBIDITY & ALGAE AT BOTTOM, NOT CURRENTLY RUN
# ============================================================================= #
# ---- PROCESSING STEPS ----
##------------------------------------------------------------------------------#
pss_stack_paths <- list.files(procrst,
                              pattern = "boats-rm\\.tif$",
                              full.names = TRUE)
# > **FOR SUN DEGLINTED IMAGERY INSTEAD
# pss_stack_paths <- list.files(procrst,
#                               pattern = "deglint\\.tif$",
#                               full.names = TRUE)
pss_df <- data.frame(file_path = pss_stack_paths) %>%
  mutate(time_string = 
           regmatches(basename(file_path), #extract datetime from filename
                      regexpr("\\d{8}_\\d{4}", basename(file_path))))
pss_datetimes <- pss_df$time_string %>%
  unique(na.omit(pss_datetimes))

message(paste0(">> Beginning Water Column Correction pipeline: ",
               format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "..."))
inittime <- Sys.time()
for (value in pss_datetimes) {
  ##----------------------------------------------------------------------------##
  ## PARAMETERIZING ITERATION LABELS ----
  ##----------------------------------------------------------------------------##
  file_ts <- value #only HHmm 
  file_timestamp <- value %>%
    as.POSIXct(format="%Y%m%d_%H%M", tz="UTC")
  stackdate <- as_date(substr(value, 1, nchar(value) - 5)) #YMD
  starttime <- Sys.time()
  message(paste0(">> STARTING: WATER COLUMN PROCESSING, PSS ", 
                 file_ts, " | ", 
                 format(Sys.time(), "%T %Z")))
  # Point to existing pre-processed imagery file locations
  mosaic_filepath    <- paste0(procrst, "/", file_ts, 
                               "_pss8b_mcf-norm01_normalized.tif")
  cleaned_filepath   <- paste0(procrst, "/", file_ts,
                               "_pss8b_mcf-norm02_boats-rm.tif")
  # **deglinted_filepath <- paste0(procrst, "/", file_ts,
  #                              "_pss8b_mcf-norm03_deglint.tif")
  ##----------------------------------------------------------------------------##
  ## SAGAWA BENTHIC CORRECTION ----
  ##----------------------------------------------------------------------------##
  #  Rw    : Surface reflectance (glint-corrected band)
  #  Rinf  : Deep water reflectance (the 1st percentile proxy)
  #  k     : Attenuation coefficient (from your sand polygons)
  #  depth : Tide-adjusted LiDAR or Adjusted SDB
  message("   >> Performing Sagawa 8-band correction...\n      ",
          format(Sys.time(), "%T %Z"))
  # LOAD DATASETS
  img_norm <- rast(mosaic_filepath)
  img_clean <- rast(cleaned_filepath)
  # img_deglint <- rast(deglinted_filepath) #**DECIDE IF TO USE DEGLINT (See L073)
  #  > Rw = normalized + cleaned [+ sun-glint-corrected] PSS mosaic 8-band imagery
  #    **DECIDE WHICH IMAGERY STACK TO USE: CLEANED VS DEGLINTED FOR ALL WCC
  img_focal <- img_clean
  #  > k = shapefile of hand-drawn polygons around known sand areas (no SAV)
  sand_shps <- st_read(sand_poly_path) %>% 
    st_transform(crs(proj_crs))
  names(sand_shps) <- c("Shape_length","Sand Polygons", "geometry")
  sand_poly <- sand_shps[2]
  #  > depth = bathymetry raster stack (built in prior script)
  barnegat_depth_stack <- rast(file.path(procdep, paste0(
    "TCB_allDepths_", file_ts, ".tif")))
  tad_topo <-  barnegat_depth_stack[[3]] #tide-adjusted depth
  depth_rast <- tad_topo
  depth_aligned <- resample(depth_rast, img_focal, method = "bilinear")
  ext(depth_aligned) <- ext(img_clean)
  print(paste("Extents match?",
              ext(depth_aligned) == ext(img_focal[[1]])))
  plot(depth_aligned)
  plot(img_focal[[4]])
  # EXTRACT R_INFINITY (Deep Water Proxy)
  #  > Calculate the 1st percentile (darkest pixels) for each band
  r_inf_vec <- as.numeric(unlist(global(img_focal, fun = quantile, 
                                        probs = 0.01, na.rm = TRUE)))
  # CALIBRATE k (DECAY) USING SAND POLYGONS
  # Extract data from ground-referenced sand areas: reflectance & depth
  sand_data <- terra::extract(c(img_focal, depth_aligned), 
                              sand_poly, df=TRUE, na.rm=TRUE)
  colnames(sand_data) <- c("ID", paste0("B", 1:8), "Depth")
  # Filter for depths where the bottom signal is reliable (0.5m to 4m)
  sand_vals <- sand_data %>% filter(Depth > 0.5 & Depth < 4)
  sagawa_bands <- list()
  k_summary <- data.frame(Band = paste0("B", 1:8), 
                          k = numeric(8), R2 = numeric(8))
  # SAGAWA CORRECTION LOOP
  for(i in 1:8) {
    band_name <- paste0("B", i)
    r_inf_val <- r_inf_vec[i]
    # CRITICAL FILTER: Remove invalid values before log()
    # Removes pixels where reflectance is lower than the deep water baseline
    #  > Only want pixels where (Reflectance - R_inf) is a positive number
    #  > Also must ensure Depth is valid (>0) to avoid Inf.
    cal_data <- sand_vals %>% 
      filter(.data[[band_name]] > r_inf_val) %>% 
      filter(Depth > 0)
    # CONDITIONAL REGRESSION
    #  > Need a minimum amount of valid sand pixels (e.g., 50) 
    #  > And the band shouldn't be NIR (B7, B8) as it rarely works for Sagawa
    if(nrow(cal_data) > 50 & i < 8) {
      tryCatch({
        y_val <- log(cal_data[[band_name]] - r_inf_val)
        x_val <- cal_data$Depth
        # Run the linear model
        fit <- lm(y_val ~ x_val)
        k_val <- -coef(fit)[2] / 2
        r2_val <- summary(fit)$r.squared
        # Check for 'Physics Violations' (k should be positive)
        if(k_val <= 0) {
          k_val <- 0.0001 # Near-zero decay fallback
          r2_val <- 0
        }
      }, error = function(e) {
        k_val <- 0.5; r2_val <- 0 # Error fallback
      })
    } else {
      # DEFAULT FALLBACK for B8 or failed regressions (if band is too absorbed)
      #  > Set a high k-value for NIR because it attenuates almost instantly
      k_val <- ifelse(i >= 8, 2.5, 0.5)
      r2_val <- NA
    }
    # CALCULATE STABILIZED SAGAWA CORRECTION
    # Calculate Two-Way Transmittance (the exponential decay)
    #  > This represents the light surviving the trip to the bottom and back
    transmittance <- exp(-2 * k_val * depth_aligned)
    # Stability Buffer: If transmittance is < 5%, we stop boosting the signal
    #  > Prevents dividing by near-zero in deep water
    transmittance <- clamp(transmittance, lower = 0.05, upper = 1.0)
    # Apply Sagawa Formula: Rb = [Rw - Rinf * (1 - transmittance)] / transmittance
    #  > Numerator: Remove the light reflected by the water column itself 
    #  > Denominator: Compensate for the light lost during travel (aka depth)
    Rb <- (img_focal[[i]] - r_inf_val * (1 - transmittance)) / transmittance
    # Final Physical Clamp
    #  > Reflectance is a ratio (0-1)
    #  > Coastal reflectance rarely exceeds 0.6; clamp at 0.8 for safety
    sagawa_bands[[i]] <- clamp(Rb, lower = 0, upper = 0.8)
    # LOG TO MATRIX
    k_summary$k[i] <- k_val
    k_summary$R2[i] <- summary(fit)$r.squared
    message(paste0("   >> Sagawa Correction for ", 
                   band_name, ": Light decay constant (k) = ", 
                   round(k_val, 2), " | R2 = ", 
                   round(summary(fit)$r.squared, 2),
                   "\n      ", format(Sys.time(), "%T %Z")))
  }
  # WRAP OUTPUT INTO ONE MULTIBAND
  sagawa_corrected <- rast(sagawa_bands)
  # RENAME FOR CLARITY: <original band name>_Sc = Sagawa corrected
  names(sagawa_corrected) <- paste0(names(img_focal), "_Sc")
  print(list(Raster = sagawa_corrected, Metadata = k_summary, Rinf = r_inf_vec))
  plot(sagawa_bands[[4]], main = paste0("Sagawa-corrected B4, ", stackdate))
  # CALCULATE STATISTICS: Summarize parameters from each step of processing
  #  > AVERAGE BRIGHTNESS: Calculate the mean for every band in the stack
  #  > 1. Normalized Reflectance (Initial PSS Mosaic)
  band_avgs_rnorm <- global(img_norm, fun = "mean", na.rm = TRUE)
  avg_reflect_rnorm <- as.vector(unlist(band_avgs_rnorm))
  names(avg_reflect_rnorm) <- names(img_norm)
  #  > 2. Boat & Wake Removal (Small Bright Outliers masked)
  band_avgs_noboat <- global(img_clean, fun = "mean", na.rm = TRUE)
  avg_reflect_noboat <- as.vector(unlist(band_avgs_noboat))
  names(avg_reflect_noboat) <- names(img_clean)
  #  > 3. **Sun Glint Correction (Remaining Bright Pixels corrected)
  # band_avgs_deglint <- global(img_deglint, fun = "mean", na.rm = TRUE)
  # avg_reflect_deglint <- as.vector(unlist(band_avgs_deglint))
  # names(avg_reflect_deglint) <- names(img_deglint)
  #  > 4. Sagawa-Corrected Reflectance
  band_avgs_saga <- global(sagawa_corrected, fun = "mean", na.rm = TRUE)
  avg_reflect_saga <- as.vector(unlist(band_avgs_saga))
  names(avg_reflect_saga) <- names(sagawa_corrected)
  # DIAGNOSTIC PLOTS
  # Plots the relationship between depth and your final Sagawa Green
  sc_diag_df <- spatSample(c(sagawa_corrected[[4]], depth_aligned), 
                           2000, na.rm=TRUE)
  colnames(sc_diag_df) <- c("Rb_Green", "Z")
  #  > If the line trends UP or DOWN at great depths, the constants are wrong
  p_sag_stab <- ggplot(sc_diag_df, aes(x = Z, y = Rb_Green)) +
    geom_point(alpha = 0.2, color = "darkgreen") +
    geom_smooth(method = "gam", color = "red") +
    xlim(0,4) +
    labs(title = "Stability Check: Depth vs. Reflectance Correction (Green Band)",
         subtitle = paste0("R_bottom = R_obs * exp(k * depth) ",
                           "(k = ", round(k_summary$k[4],3),
                           ", R2 = ", round(k_summary$R2[4], 2), ")"),
         x = "Depth (m)", y = "Corrected Benthic Reflectance") +
    theme_minimal()
  # Compare Normalized Original vs Corrected Imagery (band 4)
  iband <- 4
  p_cleancomp <- ggplot() +
    geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
    geom_spatraster(data = img_clean[[iband]], 
                    maxcell = 5e5) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    scale_fill_viridis_c(option = "mako", 
                         # name = "Normalized \nReflectance",
                         guide = "none",
                         na.value = "transparent", direction = 1,
                         limits = c(0,1)) +
    labs(title = "Normalized Reflectance", 
         subtitle = paste0("Avg.Reflectance: ", 
                           round(avg_reflect_rnorm[iband], 3))) +
    theme_void()
  # **Compare Deglinted vs Corrected Imagery (band 4)  
  # p_deglcomp <- ggplot() +
  #   geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
  #   geom_spatraster(data = img_deglint[[iband]], 
  #                   maxcell = 5e5) +
  #   geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
  #           fill = "transparent") +
  #   scale_fill_viridis_c(option = "mako", 
  #                        name = "Normalized \nReflectance",
  #                        # guide = "none",
  #                        na.value = "transparent", direction = 1,
  #                        limits = c(0,1)) +
  #   labs(title = "Sun Glint Correction", 
  #        subtitle = paste0("Avg.Reflectance: ", 
  #                          round(avg_reflect_deglint[iband], 3)))+
  #   # , caption = paste0("Calibrated using NIR band; slope = ", 
  #   #                  round(slope, 2))) + theme_void()
  p_sagawa <- ggplot() +
    geom_sf(data = water_poly, fill = "grey60", color = "grey60") +
    geom_spatraster(data = sagawa_corrected[[iband]]) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    scale_fill_viridis_c(option = "mako", 
                         # name = "Normalized \nReflectance",
                         guide = "none",
                         na.value = "transparent", direction = 1,
                         limits = c(0,1)) +
    labs(title = "Depth-Calibrated Bottom Reflectance",
         subtitle = paste0("Avg.Reflectance: ", 
                           round(avg_reflect_saga[iband], 3)), 
         caption = paste0("Calibrated using Band ", iband, "; slope (k) = ", 
                          round(k_summary$k[iband], 2))) +
    # subtitle = paste0("R2 = ", round(summary(fit)$r.squared, 3))
    # ) +
    theme_void()
  #  >  Join plots and include relevant statistics
  # diagnostics_sagawa <- p_cleancomp + p_deglcomp + p_sagawa +
  # diagnostics_sagawa <- p_sag_stab + p_sagawa +
  diagnostics_sagawa <- p_cleancomp + p_sagawa + p_sag_stab + #replaced p_deglcomp
    # plot_layout(guides = 'collect') + 
    plot_annotation(title = 
    "WCC Diagnostics 05 | Sagawa Correction of Benthic Reflectance",
    subtitle = paste0(
    "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
                      as_date(stackdate), 
                      "\n  >> SHOWING: Band ", iband)) +
    # theme_void()
    theme_minimal()
  print(diagnostics_sagawa)
  # EXPORT GRAPHIC: Summary diagnostics
  ggsave(file.path(procimg,
                   paste0("05_WCC-diagnostics_sagawa-correction_B", iband, "_",
                          file_ts, ".png")), diagnostics_sagawa,
         width = 11, height = 9, dpi = 300)
  
  rm(sc_diag_df, sagawa_bands, transmittance, tad_topo, sand_data)
  ##----------------------------------------------------------------------------##
  ## WCC PROCESSING: LYZENGA DEPTH-INVARIANT INDEX (DII) ----
  ##----------------------------------------------------------------------------##
  message("   >> Creating Lyzenga Depth-Invariant Index (DII)...\n      ",
          format(Sys.time(), "%T %Z"))
  # CALCULATE DEPTH-INVARIANT INDEX: Creates band independent of depth (in theory)
  #  > Sourcing: Lyzenga (1981) - "Remote sensing of bottom reflectance..."
  #  > Lyzenga uses the Natural Log of a band ratio, typically: 
  #    > Blue (B2) and Green (B4) for clear water 
  #    > Green (B4) and Red (B6) for turbid water
  #  --> *USE TURBIDITY INDEX FROM ABOVE TO DETERMINE PATHWAY
  # LOAD IMAGERY & SAND: Load normalized PSS mosaic and a hand-drawn sand polygon
  # EXTRACT LOG-TRANSFORMED VALUES: Pull pixels from inside the sand polygon
  sand_dii <- terra::extract(img_focal, sand_poly)
  message(paste0("   >> Calibrating Lyzenga DII for CLEAR water... \n      ", 
                 format(Sys.time(), "%T %Z")))
  # CREATE DATAFRAME: Calculate standard Blue-Green DII ratio as separate layer
  #  > Prep data for regression (Log-Linear) & add constant to avoid log(0)
  lyz_df <- data.frame(
    ln_Blue  = log(sand_dii$blue + 0.001), # Use your Band 2 Name
    ln_Green = log(sand_dii$green + 0.001)  # Use your Band 4 Name
  ) %>% 
    filter(is.finite(ln_Blue) & is.finite(ln_Green))
  # CALCULATE THE ATTENUATION RATIO (SLOPE): Build linear regression model
  #  > Slope of the regression line for ln(Blue) vs ln(Green)
  #  > Represents a ratio of attenuation coefficients (ki/kj)
  #  > Produces a "dimensionless" value to map Bottom Type (aka Benthic Cover)
  # RUN LINEAR MODEL: Regress Y (blue) on X (green) to find the slope
  fit_lyz <- lm(ln_Blue ~ ln_Green, data = lyz_df)
  slope_bg <- coef(fit_lyz)[2] # This is ki/kj
  message(paste("   >> Attenuation slope for low turbidity (blue/green):", 
                round(slope_bg, 3)))
  # APPLY DII INDEX: Creates new "band" for the bottom but with depth removed
  #  > Formula: DII = ln(Blue) - slope * ln(Green)
  lyzenga_dii <- log(img_focal$blue + 0.001) - 
    (slope_bg * log(img_focal$green + 0.001))
  names(lyzenga_dii) <- "Lyzenga_DII_SandCalibrated"
  # NORMALIZE INDEX: Get stats while ignoring NA (masks)
  #  > Use the 2nd and 98th percentile to get a "robust" range
  stats_dii <- global(lyzenga_dii, fun=quantile, 
                      probs=c(0.02, 0.98), na.rm=TRUE)
  dii_min <- stats_dii[1,1]
  dii_max <- stats_dii[1,2]
  #  > Apply Min-Max Stretch
  dii_norm <- (lyzenga_dii - dii_min) / (dii_max - dii_min)
  #  > Clamp values to 0-1 to handle the outliers we ignored
  dii_norm <- clamp(dii_norm, lower=0, upper=1)
  # CREATE DIAGNOSTIC PLOTS: Map of the index across the bay & show regression
  p_dii <- ggplot() +
    geom_sf(data = water_poly, fill = "grey80", color = "grey60") +
    geom_spatraster(data = dii_norm, 
                    maxcell = 5e5) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    scale_fill_viridis_c(option = "cividis", name = "Normalized\nIndex",
                         # trans = "reverse",
                         na.value = "transparent", direction = 1) +
    labs(title = paste(
      "Proxy of Benthic Brightness"),
      subtitle =  "Lyzenga DII = ln(Blue) - slope * ln(Green)") +
    theme_void() 
  # Linear regression results
  #  > Should have linear "cigar" shape
  #  > if it looks like a blob, sand polygon may be problematic;
  #    could include other bottom types or depths are too uniform
  p_lyz <- ggplot(lyz_df, aes(x = ln_Green, y = ln_Blue)) +
    geom_point(alpha = 0.1, color = "darkblue") +
    geom_smooth(method = "lm", color = "red") +
    labs(title = "Lyzenga Bi-Plot Calibration",
         subtitle = paste("Slope (ki/kj) =", round(slope_bg, 3)),
         x = "ln(Green Band)", y = "ln(Blue Band)") +
    theme_minimal()
  # Join plots and include relevant statistics
  diagnostics_lyzenga <- p_lyz + p_dii + 
    # plot_layout(guides = 'collect') +
    plot_annotation(title = paste0(
      "WCC Diagnostics 06 | Lyzenga Depth Invariant Index (DII)"),
      subtitle = paste0(
        "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
        as_date(stackdate))) +
    theme_void()
  print(diagnostics_lyzenga)
  # EXPORT GRAPHIC: Summary diagnostics
  ggsave(file.path(procimg,
                   paste0("06_WCC-diagnostics_lyzenga-calibration_",
                          file_ts, ".png")), diagnostics_lyzenga,
         width = 11, height = 9, dpi = 300)
  ##----------------------------------------------------------------------------##
  ## PREPARING FINAL EXPORTS ----
  ##----------------------------------------------------------------------------##
  message("   >> Building corrected water column stack...\n      ",
          format(Sys.time(), "%T %Z"))
  # CLIP ALL LAYERS TO SAV DEPTH (0-3M)
  # *APPLY TO STACK: can mask murky areas before the Random Forest runs (optional)
  # CREATE WCC IMAGERY STACK: Combine WCC layers to create depth-corrected stack
  wcc_stack <- c(sagawa_corrected, dii_norm) #9 layers
  names(wcc_stack)[nlyr(wcc_stack)] <- c("Lyzenga_DII")
  wcc_filename <- paste0("wcc-stack_allDepths_", file_ts, ".tif") 
  writeRaster(wcc_stack, 
              file.path(procwcc, wcc_filename),
              gdal = c("COMPRESS=LZW", "TFW=YES"),
              overwrite=TRUE)
  # CLEAN UP: Remove temporary objects & trigger garbage collection
  #  > Leave only layers needed for next steps
  rm(lyzenga_dii, stats_dii, lyz_df, 
     sagawa_corrected, dii_norm, wcc_stack)
  gc()
  endtime <- Sys.time()
  loop_time <- round(as.numeric(difftime(endtime, starttime, units = "mins")), 2)
  message(paste0(">> FINISHED: WATER COLUMN CORRECTION PROCESSING, PSS ", 
                 file_ts, " | ", 
                 format(Sys.time(), "%T %Z"), "\n   Duration: ", loop_time,
                 " minutes", "\n  ----------------------------------------"))
  ##----------------------------------------------------------------------------##
  beep(2) #----
  ##----------------------------------------------------------------------------##
}
elapsed <- round(as.numeric(difftime(endtime, inittime, units = "mins")), 2)
message(paste0(">> Depth & Water Quality Processing completed for all dates on
    ", format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "  |  Elapsed Time: ", 
               elapsed, " minutes"))
Sys.sleep(2)
beep(4)

# ============================================================================= #
# ** ADDITIONAL STEPS, NOT RUN (STILL IN DEVELOPMENT) **
##----------------------------------------------------------------------------##
## FLOATING ALGAE DETECTION <NOT RUN> ----
##----------------------------------------------------------------------------##
# # LOGIC: Uses a "Surface Vegetation" difference: NIR - RedEdge
# message("   >> Detecting floating macroalgae...\n      ",
#         format(Sys.time(), "%T %Z"))
# # **Using the 'cleaned' mosaic imagery (not deglinted, only boats removed)
# img_clean <- rast(cleaned_filepath)
# # CALCULATE INDICES
# #  > Define bands (from PSS 8-Band Mosaic)
# blue      <- img_clean[[2]] #Band 2 (Blue)
# green     <- img_clean[[4]] #Band 4 (Green)
# red       <- img_clean[[6]] #Band 6 (Red)
# red_edge  <- img_clean[[7]] #Band 7 (Red Edge)
# nir       <- img_clean[[8]] #Band 8 (NIR)
# # A. Floating Algae Index (FAI): The "Baseline" shape
# #  Measures height of NIR peak above the Red-SWIR baseline
# #  > Planet has no SWIR, so use a Red-NIR baseline approximation
# #  > Q: Is there a spectral peak in NIR? (Physical presence of surface biomass)
# fai <- nir - (red + (nir - red) * 0.5) 
# # B. Surface Algal Bloom Index (SABI): The "Contrast" check
# #  > Excellent for separating surface vegetation from water
# #  > Q: Is object brighter in NIR than in B/G? (Optical contrast against water)
# sabi <- (nir - red) / (blue + green + 0.0001)
# # C. Normalized Difference Red Edge (NDRE): The "Chlorophyll" check
# #  (NIR - RedEdge) / (NIR + RedEdge) 
# #  > Floating algae will have distinct relationship compared to SAV
# #  > Q: Does vegetation edge shape indicate healthy surface chlorophyll?
# ndre <- (nir - red_edge) / (nir + red_edge + 0.0001)
# # NORMALIZE SCALES (0 to 1)
# #  > Use the 99th percentile to clamp outliers, ensuring a clean 0-1 scale
# normalize <- function(x) {
#   min_v <- global(x, "min", na.rm=TRUE)[1,1]
#   max_v <- global(x, fun=quantile, probs=0.99, na.rm=TRUE)[1,1]
#   return(clamp((x - min_v) / (max_v - min_v), 0, 1))
# }
# fai_norm  <- normalize(fai) 
# names(fai_norm) <- "FAI_Norm"
# sabi_norm <- normalize(sabi)
# names(sabi_norm) <- "SABI_Norm"
# ndre_norm <- normalize(ndre)
# names(ndre_norm) <- "NDRE_Norm"
# # 4. WEIGHTED VOTING SYSTEM
# # We give NDRE a strong vote because it leverages the unique PlanetScope band
# # FAI gets the highest vote as it's the standard for "Floating" detection
# algae_confA <- (fai_norm * 0.5) + (sabi_norm * 0.3) + (ndre_norm * 0.2)
# names(algae_confA) <- "Algae_Conf"
# # EXTRACT UNCERTAINTY ZONE
# # algae_uncertain <- algae_confA
# # algae_uncertain[algae_uncertain > 0.90 | algae_uncertain < 0.10] <- NA
# # 1. CALCULATE ALGAE PROBABILITY DISTRIBUTION & SPATIAL ERROR
# # Pull the specific Algae Probability layer from your RF Prediction
# prob_algae <- algae_confA
# # 2. DEFINE THRESHOLDS (Based on your Algae Confidence Logic)
# high_cutoff <- 0.75
# low_cutoff  <- 0.25
# pct90th_alg <- 0.9
# # Create Binary Mask: 1 = Clear/Reliable, 0 = High Macroalgae
# algae_reliability_mask <- prob_algae < high_cutoff
# # CALCULATE STATISTICS: Calculate the mean for every band in the stack
# #  > Weighted Algae Probability
# band_avg_algae <- global(algae_confA, fun = "mean", na.rm = TRUE)
# avg_conf_algae <- as.vector(unlist(band_avg_algae))
# # DIAGNOSTIC PLOTS: Visualize Indices and Algal Confidence
# # DENSITY "RESIDUAL" PLOT
# #  > This shows how "clean" the separation is: 
# #  > Should see 2 peaks at the ends & a 'valley' in the middle
# sample_vals <- spatSample(prob_algae, 50000, na.rm=TRUE)
# colnames(sample_vals) <- "prob"
# p_algprob <- ggplot(sample_vals, aes(x = prob)) +
#   geom_density(fill = "grey", alpha = 0.3) +
#   geom_vline(xintercept = c(low_cutoff, high_cutoff),
#              linetype = "dashed", color = "darkred", linewidth = 1) +
#   annotate("text", x = 0.1, y = 0.3, label = "Clear\nWater/Sand",
#            color = "blue") +
#   annotate("text", x = 0.9, y = 0.3, label = "Likely\nMacroalgae",
#            color = "red") +
#   annotate("text", x = 0.5, y = 0.3, label = "Uncertainty Zone\n(Mixed Pixels)",
#            color = "black") +
#   xlim(0,1) + 
#   # ylim(0,7) +
#   labs(title = "ALGAE PROBABILITY DISTRIBUTION",
#        subtitle = paste0("Avg. Weighted Score: ", 
#                          round(avg_conf_algae, 3)), 
#        x = "Probability of Macroalgae", y = "Sample Density") + theme_bw()
# # SPATIAL ERROR MASK: maps where model is uncertain (mid-range) OR certain
# # uncertainty_map <- ifel(prob_algae > low_cutoff & prob_algae < high_cutoff,
# #                         1, NA)
# uncertainty_map <- ifel(prob_algae > high_cutoff, 1, NA)
# plot(uncertainty_map , col = "darkred")
# # FAI
# p_fai <- ggplot() +
#   geom_sf(data = water_poly, color = "grey60", fill = "grey60") +
#   geom_spatraster(data = fai_norm) +
#   geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
#           fill = "transparent") +
#   scale_fill_whitebox_c(palette = "high_relief", direction = -1,
#                         # scale_fill_viridis_c(option = "rocket", direction = -1,
#                         # scale_fill_whitebox_c(palette = "viridi", direction = 1,
#                         guide = "none",
#                         limits = c(0,1),
#                         na.value = "transparent") +
#   labs(title = "FLOATING ALGAE INDEX",
#        subtitle = paste0(
#          "FAI = NIR/Red ratio")) +
#   # caption = "FLOATING ALGAE INDEX") +
#   theme_void()
# # SABI
# p_sabi <- ggplot() +
#   geom_sf(data = water_poly, color = "grey60", fill = "grey60") +
#   geom_spatraster(data = sabi_norm) +
#   geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
#           fill = "transparent") +
#   scale_fill_whitebox_c(palette = "high_relief", direction = -1,
#                         # scale_fill_viridis_c(option = "rocket", direction = -1,
#                         # scale_fill_whitebox_c(palette = "viridi", direction = 1,
#                         guide = "none",
#                         limits = c(0,1),
#                         na.value = "transparent") +
#   labs(title = "SURFACE ALGAL BLOOM INDEX",
#        subtitle = paste0(
#          "SABI = NIR-Red / Blue+Green ratio")) +
#   # caption = "SURFACE ALGAL BLOOM INDEX") +
#   theme_void()
# # NDRE
# p_ndre <- ggplot() +
#   geom_sf(data = water_poly, color = "grey60", fill = "grey60") +
#   geom_spatraster(data = ndre_norm) +
#   geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
#           fill = "transparent") +
#   scale_fill_whitebox_c(palette = "high_relief", direction = -1,
#                         # scale_fill_viridis_c(option = "rocket", direction = -1,
#                         # scale_fill_whitebox_c(palette = "viridi", direction = 1,
#                         # guide = "none", 
#                         name = "Normalized\n Values", limits = c(0,1),
#                         na.value = "transparent") +
#   labs(title = "NORM. DIFF. RED EDGE",
#        subtitle = paste0(
#          "NDRE = NIR/RedEdge ratio")) +
#   # caption = "NORMALIZED DIFFERENCE RED EDGE") +
#   theme_void()
# # p_uncertainty <- ggplot() +
# #   # Layer 1: Boat-masked classification (to provide background context)
# #   geom_spatraster(data = img_clean[[1]], alpha = 0.3, maxcell = 5e5) +
# #   scale_fill_gradient(low = "gray80", high = "gray20", 
# #                       guide = "none", na.value = "white") +
# #   # Layer 2: Algae Uncertainty Layer
# #   # We use a bright, high-contrast palette to make the uncertainty pop
# #   geom_spatraster(data = algae_uncertain, maxcell = 5e5) +
# #   scale_fill_viridis_c(option = "inferno", 
# #                        na.value = "transparent", 
# #                        name = "Uncertainty\n(30-60%)") +
# #   theme_minimal() +
# #   theme(legend.position = "right",
# #         panel.grid.major = element_line(color = "gray90", 
# #                                         linetype = "dashed")) +
# #   labs(title = "Macroalgae Classification Uncertainty",
# #        subtitle = "Highlighting zones of spectral confusion in prediction",
# #        x = "Longitude", y = "Latitude")
# # WEIGHTED CONFIDENCE SCORE
# p_algconf <- ggplot() +
#   geom_sf(data = water_poly, color = "grey60", fill = "grey60") +
#   geom_spatraster(data = algae_confA) +
#   geom_spatraster_contour(data = prob_algae, breaks = pct90th_alg, 
#                           alpha = 0.6,
#                           color = "darkred", linewidth = 1) +
#   geom_sf(data = water_poly, color = "grey40", linewidth = 0.4,
#           fill = "transparent") +
#   # , aes(fill = FAI_Norm)) +
#   scale_fill_hypso_c(palette = "arctic", direction = 1,
#                      # scale_fill_whitebox_c(palette = "deep", direction = 1,
#                      # scale_fill_viridis_c(option = "rocket", direction = -1,
#                      name = "Algae\nProb", limits = c(0,1),
#                      na.value = "transparent") +
#   labs(title = "CLASSIFICATION CONFIDENCE",
#        subtitle = paste0(
#          "Weighted Scores across FAI, SABI, & NDRE")) +
#   # caption =  +
#   theme_void()
# #  >  Join plots and include relevant statistics
# # DIAGNOSTIC A
# diagnostics_algInd <- p_fai + p_sabi + p_ndre +
#   plot_annotation(title = paste0(
#     "SAV Diagnostics 03a  |  Macroalgae Indices"), 
#     subtitle = paste0(
#       "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
#       stackdate)) +
#   theme_void()
# print(diagnostics_algInd)
# # DIAGNOSTIC B
# diagnostics_algae <- p_algprob + p_algconf +
#   # plot_layout(guides = 'collect') +
#   plot_annotation(title = paste0(
#     "WCC Diagnostics 03b  |  Macroalgae Classification & Confidence"), 
#     subtitle = paste0(
#       "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
#       stackdate),
#     caption = "RED = Confidence above 90th Percentile Threshold") +
#   theme_void() 
# print(diagnostics_algae)
# # EXPORT GRAPHICS: Summary diagnostics
# ggsave(file.path(procimg, 
#                  paste0("03a_SAV-diagnostics_algae-indices_",   
#                         file_ts, ".png")), diagnostics_algInd,
#        width = 11, height = 9, dpi = 300)
# ggsave(file.path(procimg, 
#                  paste0("03b_SAV-diagnostics_algae-conf-class_",   
#                         file_ts, ".png")), diagnostics_algae,
#        width = 11, height = 9, dpi = 300)
# # CREATE ALGAE STACK: Combine indices & confidence scores
# algae_stack <- c(fai_norm, sabi_norm, ndre_norm, algae_confA) # layers
# names(algae_stack) <- c("Algae_FAI", 
#                         "Algae_SABI", 
#                         "Algae_NDRE", 
#                         "Algae_Confidence")
# plot(algae_stack)
# # EXPORT OUTPUTS (ALL DEPTHS)
# alg_filename <- paste0("algae-stack_", file_ts, ".tif") 
# writeRaster(algae_stack, file.path(procwcc, alg_filename), 
#             overwrite = TRUE)
# # CLEAN UP: Remove temporary objects & trigger garbage collection
# #  > Leave only layers needed for Diagnostics Plots
# rm(fai, sabi, ndre, fai_norm, sabi_norm, ndre_norm, algae_confA, prob_algae,
#    algae_reliability_mask, uncertainty_map, sample_vals, algae_stack) 
# gc()
##----------------------------------------------------------------------------##
## TURBIDITY INDEX <NOT RUN>  ----
##----------------------------------------------------------------------------##
# message("   >> Calculating turbidity index...\n      ",
#         format(Sys.time(), "%T %Z"))
# # Point to existing pre-processed imagery file locations
# # **Using the deglinted mosaic imagery to calculate predicted depth
# # img_deglint <- rast(deglinted_filepath)
# img_focal
# # CALCULATE TURBIDITY INDEX: New layer based on raster calculation
# #  > Normalized Difference Turbidity Index (NDTI) = Red/Green Ratio
# #  > Higher values = more suspended sediment
# #  > Formula: (Red - Green) / (Red + Green)
# turbidity_idx <- (img_clean[[6]] - img_clean[[4]]) /
#   (img_clean[[6]] + img_clean[[4]])
# names(turbidity_idx) <- "Turbidity_Index"
# # GENERATE STATISTICS
# #  > Use 2nd & 98th percentile as anchors (avoids scaling on extreme outliers)
# stats_turb <- global(turbidity_idx, fun = quantile, 
#                      probs = c(0.02, 0.98), na.rm = TRUE)
# t_min <- stats_turb[1,1]
# t_max <- stats_turb[1,2]
# t_avg <- ((t_min + t_max) / 2)
# # STANDARDIZE THE INDEX
# # Apply Min-Max Scaling (Formula: (x - min) / (max - min)) & clamp to 0-1 range
# turbidity_index_std <- clamp((turbidity_idx - t_min) / (t_max - t_min), 0, 1)
# names(turbidity_index_std) <- "Standardized_Turbidity"
# plot(turbidity_index_std)
# # CALCULATE MEAN TURBIDITY: Determine single value for entire study area
# avg_turbidity <- global(turbidity_index_std, "mean", na.rm = TRUE)$mean
# # DEFINE TURBIDITY RELIABILITY THRESHOLD: Calculate 90th percentile of turbidity
# #  > Anything in the top 10% of "murkiness" is considered unreliable
# # turb_threshold <- global(turbidity_index_std, fun=quantile,
# #                          probs=0.90, na.rm=TRUE)[1,1]
# turb_threshold <- 0.9
# # DEFINE CONFIDENCE FLAG
# if(avg_turbidity > turb_threshold) {
#   message(paste0(
#     "   !! WARNING !! Turbidity above acceptable limits for Benthic Mapping (",
#     round(avg_turbidity, 3), 
#     "):\n   HIGH TURBIDITY protocol applied to Water Column Correction for ",
#     file_ts))
#   is_reliable <- FALSE
# } else {
#   message(paste0(
#     "  >> Turbidity within acceptable limits for Benthic Mapping (", 
#     round(avg_turbidity, 3), 
#     "):\n     CLEAR protocol applied to Water Column Correction \n     ", 
#     format(Sys.time(), "%T %Z")))
#   is_reliable <- TRUE
# }
# # Create Binary Mask: 1 = Clear/Reliable, 0 = Too Murky
# t_reliability_mask <- turbidity_index_std < turb_threshold
# plot(t_reliability_mask)
# # CALCULATE "PERCENT DATA LOSS"
# t_stats <- freq(t_reliability_mask)
# # 0 is the "Murky" class, 1 is the "Clear" class
# clear_pixels <- t_stats$count[t_stats$value == 1]
# murky_pixels <- t_stats$count[t_stats$value == 0]
# pct_lost_turb <- (murky_pixels / (clear_pixels + murky_pixels)) * 100
# message(paste0("   >> Quality Check: ", round(pct_lost_turb, 2), 
#                "% of habitat marked by high turbidity for ", 
#                file_ts))
# # VISUALIZE RESULTS: Create maps of turbidity and reliability
# #  > Map reliability as a binary: show exactly where seagrass data is "Trusted"
# p_t_reliability <- ggplot() +
#   geom_sf(data = water_poly, fill = "grey80", color = "grey60") +
#   geom_spatraster(data = t_reliability_mask) +
#   # geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
#   #         fill = "transparent") +
#   scale_fill_gradient(low = "red4", high = "goldenrod",
#                       na.value = "transparent", guide = "none") +
#   labs(title = "RELIABILITY THRESHOLD",
#        subtitle = paste0("RED = High Turbidity (above 90th Percentile)",
#                          "\n             >> ", round(pct_lost_turb, 2),
#                          "% of total area obscured")) +
#   theme_void()
# #  > Map NDTI across extent
# p_turbidity <- ggplot() +
#   geom_sf(data = water_poly, fill = "grey80", color = "grey60") +
#   geom_spatraster(data = turbidity_index_std, 
#                   aes(fill = Standardized_Turbidity)) +
#   geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
#           fill = "transparent") +
#   # scale_fill_whitebox_c(palette = "muted", direction = 1, 
#   scale_fill_viridis_c(option = "mako", direction = -1,
#                        name = "Normalized\nTurbidity",
#                        na.value = "transparent") +
#   labs(title = "NDTI",
#        subtitle = paste0("Red/Green ratio = suspended sediment proxy", 
#                          "\n   Avg. Turbidity: ", 
#                          round(avg_turbidity, 2))) +
#   # caption =  +
#   theme_void()
# #  >  Join plots and include relevant statistics
# diagnostics_turbidity <- p_turbidity + p_t_reliability +
#   # plot_layout(guides = 'collect') +
#   plot_annotation(title = 
#                     "SAV Diagnostics 04  |  Water Column Turbidity & Reliability",
#                   subtitle = paste0(
#                     "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
#                     stackdate)) +
#   # caption = paste0("Data Quality: ", round((100-pct_lost_turb), 2), 
#   # "% of standardized pixels fall below the threshold for this date")) +
#   theme_void()
# print(diagnostics_turbidity)
# # EXPORT GRAPHIC: Summary diagnostics
# ggsave(file.path(procimg,
#                  paste0("04_SAV-diagnostics_turbidity-reliability-check_",
#                         file_ts, ".png")), diagnostics_turbidity,
#        width = 11, height = 9, dpi = 300)
# # # In normalized 0-1 space, a Turbidity Index > 0.25 is often the tipping point
# # #   where bottom reflectance is lost to backscatter.
# # # Alternatively, use the 90th percentile for a data-driven approach:
# # turb_threshold <- global(turbidity_idx, fun=quantile, 
# #                          probs=0.90, na.rm=TRUE)[1,1]
# # Create Binary Mask: 1 = Clear/Reliable, 0 = Too Murky
# p_turbidityQC <- ggplot() +
#   geom_sf(data = water_poly, fill = "grey80", color = "grey60") +
#   geom_spatraster(data = turbidity_index_std) +
#   # Highlight the "Unreliable" zones with a red contour or hash
#   geom_spatraster_contour(data = turbidity_index_std, breaks = turb_threshold, 
#                           color = "red2", linewidth = 0.2) +
#   scale_fill_viridis_c(option = "mako", name = "Normalized\nTurbidity",
#                        na.value = "transparent", direction = -1) +
#   labs(title = paste0("SUPP 03  |  Normalized Turbidity & Reliability for ", 
#                       stackdate),
#        subtitle = paste0("Average turbidity: ", round(avg_turbidity, 2), 
#                          # "  |  90th Percentile Threshold: ",
#                          # round(turb_threshold, 2), 
#                          "\n",
#                          "RED = areas above threshold (bottom reflectance may be lost to backscatter)"),
#        caption = paste0("*Data Quality: ", round(pct_lost_turb, 2), 
#                         "% of pixels obscured by sediment on this date")) +
#   theme_minimal(); p_turbidityQC
# ggsave(file.path(procimg,
#                  paste0("S03_SAV-diagnostics_turbidity-reliability-checkQC_",
#                         file_ts, ".png")), p_turbidityQC,
#        width = 8, height = 14, dpi = 300)
# # CREATE ALGAE STACK: Combine indices & confidence scores
# turb_stack <- c(turbidity_index_std, t_reliability_mask) # layers
# names(turb_stack) <- c("Turbidity", 
#                        "Turb_Reliability")
# plot(turb_stack)
# # EXPORT OUTPUTS (ALL DEPTHS)
# turb_filename <- paste0("turbidity-stack_", file_ts, ".tif") 
# writeRaster(turb_stack, file.path(procwcc, turb_filename), 
#             overwrite = TRUE)
# 
# 
# # *APPLY TO STACK: can mask murky areas before the Random Forest runs (optional)
# # img_cleaned <- mask(r_stack, t_reliability_mask) #r_stack=load of FinalStack
# 
# # CLEAN UP: Remove temporary objects & trigger garbage collection
# #  > Leave only layers needed for next steps
# rm(turbidity_idx, turbidity_index_std, t_reliability_mask, turb_stack)
# gc()
# 
# # DETERMINE TURBIDITY PROCESSING 
# if(avg_turbidity < turb_threshold) {
#   message(paste0(
#     "   >> Calibrating Lyzenga DII for CLEAR water... \n      ", 
#         format(Sys.time(), "%T %Z")))
#   is_reliable <- TRUE
#   # CREATE DATAFRAME: Calculate standard Blue-Green DII ratio as sep layer
#   #  > Prep data for regression (Log-Linear) & add constant to avoid log(0)
#   lyz_df <- data.frame(
#     ln_Blue  = log(sand_dii$blue + 0.001), # Use your Band 2 Name
#     ln_Green = log(sand_dii$green + 0.001)  # Use your Band 4 Name
#   ) %>% 
#     filter(is.finite(ln_Blue) & is.finite(ln_Green))
#   # CALCULATE THE ATTENUATION RATIO (SLOPE): Build linear regression model
#   #  > Slope of the regression line for ln(Blue) vs ln(Green)
#   #  > Represents a ratio of attenuation coefficients (ki/kj)
#   #  > Produces a "dimensionless" value to map Bottom Type (aka Benthic Cover)
#   # RUN LINEAR MODEL: Regress Y (blue) on X (green) to find the slope
#   fit_lyz <- lm(ln_Blue ~ ln_Green, data = lyz_df)
#   slope_bg <- coef(fit_lyz)[2] # This is ki/kj
#   message(paste("   >> Attenuation slope for low turbidity (blue/green):", 
#                 round(slope_bg, 3)))
#   # APPLY DII INDEX: Creates new "band" for the bottom but with depth removed
#   #  > Formula: DII = ln(Blue) - slope * ln(Green)
#   lyzenga_dii <- log(img_focal$blue + 0.001) - 
#     (slope_bg * log(img_focal$green + 0.001))
#   names(lyzenga_dii) <- "Lyzenga_DII_SandCalibrated"
#   # NORMALIZE INDEX: Get stats while ignoring NA (masks)
#   #  > Use the 2nd and 98th percentile to get a "robust" range
#   stats_dii <- global(lyzenga_dii, fun=quantile, 
#                       probs=c(0.02, 0.98), na.rm=TRUE)
#   dii_min <- stats_dii[1,1]
#   dii_max <- stats_dii[1,2]
#   #  > Apply Min-Max Stretch
#   dii_norm <- (lyzenga_dii - dii_min) / (dii_max - dii_min)
#   #  > Clamp values to 0-1 to handle the outliers we ignored
#   dii_norm <- clamp(dii_norm, lower=0, upper=1)
#   # CREATE DIAGNOSTIC PLOTS: Map of the index across the bay & show regression
#   p_dii <- ggplot() +
#     geom_sf(data = full_shp, fill = "grey80", color = "grey60") +
#     geom_spatraster(data = dii_norm, 
#                     maxcell = 5e5) +
#     geom_sf(data = bbay_shp, color = "grey60", linewidth = 0.4,
#           fill = "transparent") +
#     scale_fill_viridis_c(option = "cividis", name = "Normalized\nIndex",
#                          # trans = "reverse",
#                          na.value = "transparent", direction = 1) +
#     labs(title = paste(
#       "Proxy of Benthic Brightness"),
#       subtitle =  "Lyzenga DII = ln(Blue) - slope * ln(Green)") +
#     theme_void() 
#   # Linear regression results
#   #  > Should have linear "cigar" shape
#   #  > if it looks like a blob, sand polygon may be problematic...
#   #    could include other bottom types or depths are too uniform
#   p_lyz <- ggplot(lyz_df, aes(x = ln_Green, y = ln_Blue)) +
#     geom_point(alpha = 0.1, color = "darkblue") +
#     geom_smooth(method = "lm", color = "red") +
#     labs(title = "Lyzenga Bi-Plot Calibration",
#          subtitle = paste("Slope (ki/kj) =", round(slope_bg, 3)),
#          x = "ln(Green Band)", y = "ln(Blue Band)") +
#     theme_minimal()
# } else {
#   message(paste0(
#     "   !! WARNING !! Turbidity above acceptable limits for Benthic Mapping (",
#     round(avg_turbidity, 3), 
#     "):\n   HIGH TURBIDITY protocol applied to Lyzenga DII Calibration for ",
#     file_ts))
#   is_reliable <- FALSE
#   # *IF HIGH TURBIDITY, DII RUNS WITH RED/GREEN INSTEAD
#   # CREATE DATAFRAME: Calculate standard Blue-Green DII ratio as separate layer
#   #  > Prep data for regression (Log-Linear) & add constant to avoid log(0)
#   lyzTb_df <- data.frame(
#     ln_Green  = log(sand_dii$green_Sc + 0.001), #name of Band 4
#     ln_Red = log(sand_dii$red_Sc + 0.001)  #name of Band 6
#   ) %>%
#     filter(is.finite(ln_Green) & is.finite(ln_Red))
#   # CALCULATE THE ATTENUATION RATIO (SLOPE): Build linear regression model
#   #  > Slope of the regression line for ln(Green) vs ln(Red)
#   #  > Represents a ratio of attenuation coefficients (ki/kj)
#   #  > Produces a "dimensionless" value to map Bottom Type (aka Benthic Cover)
#   # RUN LINEAR MODEL: Regress Y (green) on X (red) to find the slope
#   fit_lyzTb <- lm(ln_Green ~ ln_Red, data = lyzTb_df)
#   slope_gr <- coef(fit_lyzTb)[2] # This is ki/kj
#   message(paste("   >> Lyzenga Slope for !High Turbidity! (Red/Green):",
#                 round(slope_gr, 3), "for", stackdate))
#   # APPLY DII INDEX: Creates new "band" for the bottom but with depth removed
#   #  > Formula: DII = ln(Green) - slope * ln(Red)
#   lyzenga_diiTb <- log(img_focal$green + 0.001) - 
#     (slope_gr * log(img_focal$red + 0.001))
#   names(lyzenga_diiTb) <- "Lyzenga_DII_HighTurbidity"
#   # NORMALIZE INDEX: Get stats while ignoring NA (masks)
#   #  > Use the 2nd and 98th percentile to get a "robust" range
#   stats_diiTB <- global(lyzenga_diiTB, fun=quantile, 
#                       probs=c(0.02, 0.98), na.rm=TRUE)
#   dii_min <- stats_diiTB[1,1]
#   dii_max <- stats_diiTB[1,2]
#   #  > Apply Min-Max Stretch
#   dii_normTB <- (lyzenga_diiTB - dii_min) / (dii_max - dii_min)
#   #  > Clamp values to 0-1 to handle the outliers we ignored
#   dii_normTB <- clamp(dii_norm, lower=0, upper=1)
#   # CREATE DIAGNOSTIC PLOTS: Map of the index across the bay & show regression
#   p_dii <- ggplot() +
#     geom_sf(data = full_shp, fill = "grey80", color = "grey60") +
#     geom_spatraster(data = dii_normTB, 
#                     maxcell = 5e5) +
#     geom_sf(data = bbay_shp, color = "grey60", linewidth = 0.4,
#           fill = "transparent") +
#     scale_fill_viridis_c(option = "cividis", name = "Normalized\nIndex",
#                          # trans = "reverse",
#                          na.value = "transparent", direction = 1) +
#     labs(title = paste(
#       "HIGH TURBIDITY Proxy of Benthic Brightness"),
#       subtitle =  "Lyzenga DII = ln(Green) - Slope * ln(Red)") +
#     theme_void() 
#   p_lyz <- ggplot(lyzTb_df, aes(x = ln_Red, y = ln_Green)) +
#     geom_point(alpha = 0.1, color = "darkblue") +
#     geom_smooth(method = "lm", color = "red") +
#     labs(title = "Lyzenga Bi-Plot Calibration, Turbid Conditions",
#          subtitle = paste("Slope (ki/kj) =", round(slope_gr, 3)),
#          x = "ln(Red Band)", y = "ln(Green Band)") +
#     theme_minimal()
# }