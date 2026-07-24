# ============================================================================= #
# Script Name:  03_build_dep-stack.R
# Description:  Creates water depth data layers based on tide data & PSS imagery. 
#               Compares tide-adjusted topobathy against image-based predictions.
# Author:       Jess M. Stitt
# Date:         2026-03-31
# ============================================================================= #

# ============================================================================= #
# ---- TASK DESCRIPTION ----
# Input(s):
#  > [] 
#     ()
#  > [SHP] Water Depth Raster
#     (Topobathymetric lidar model of BB-LEH, 2012 | USGS EAARL-B)
# Output(s): 
#  > [, ] 
#     ()
# ============================================================================= #
# ---- LIBRARIES SPECIFIC TO THIS TASK ----
library(curl)          #for reading URLs (web requests)
library(dataRetrieval) #for accessing USGS datasets (tidal elevation)
# ============================================================================= #
# ---- USER-DEFINED VARIABLES FOR SCRIPT ----
# <NONE>
# ============================================================================= #
# ---- PROCESSING STEPS ----
##------------------------------------------------------------------------------#
## 1. Pull in USGS Site IDs for Barnegat Bay to reference for tide height ----
##------------------------------------------------------------------------------#
# > List all relevant USGS Monitoring Station IDs
site_ids <- c("USGS-01408168", #Mantoloking
              "USGS-01408748", #BB Rt37 (HAS CONT BUT NOT DAILY DATA)
              "USGS-01409124", #Waretown
              "USGS-01409125", #Barnegat Light
              "USGS-01409146", #Ship Bottom
              "USGS-01409335") #Little Egg Inlet
# > Fetch station locations & project to CRS
site_metadata <- dataRetrieval::read_waterdata_monitoring_location(
  monitoring_location_id = site_ids)
stations_sf <- st_as_sf(site_metadata) %>% st_transform(proj_crs) %>%
  mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() #convert back to a standard data frame
write.csv(
  site_metadata, 
  file = paste0(proccsv, "/USGS_tidal-stations_BBLEH.csv"),
  row.names = FALSE
)
##------------------------------------------------------------------------------#
## 2. Alignment Check (plot water poly, topobathy, and stations together) ----
##------------------------------------------------------------------------------#
# > Check that water depth is in the correct metric CRS (UTM 18N)
if(crs(lidar_res, proj=TRUE) !=
   "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs") {
  lidar_res <- project(lidar_res, proj_crs)
}
message("Base datasets standardized to EPSG:32618 & clipped to water boundary.")
# > Visualize the clipped depth map alongside the tidal monitoring stations
align_plot <- ggplot() +
  geom_spatraster(data = topob_dem) +
  scale_fill_hypso_b(palette = "colombia_bathy",
                     n.breaks = 12,
                     name = "Water Depth (m)") +
  geom_sf(data = sav_poly_sf, fill = NA, color = "darkorange", linewidth = 0.6) +
  geom_spatvector(data = water_poly, fill = "transparent", color = "black", 
                  linewidth = 0.8) +
  # > Add the Station Points (green dots with black outline)
  geom_point(data = stations_sf, aes(x = X, y = Y),
             shape = 21, fill = "green", color = "black", size = 3,
             stroke = 1) +
  # > Optional: Label the points with their observed tide value
  geom_label(data = stations_sf,
             aes(x = X, y = Y, label = monitoring_location_id),
             nudge_y = -1200, nudge_x = -1200,
             color = "black", fill = "grey80", alpha = 0.7,
             size = 2.5) +
  labs(
    title = "Lidar-Derived Bathymetry, USGS EAARL-B",
    subtitle = "Barnegat Bay-Little Egg Harbor: EPSG 32618",
    caption = "Orange: SAV Suitability Zone (0-2m) | Green: USGS Tide Stations") +
  theme(legend.position = "right", plot.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "gray90", linetype = "dashed")) +
  theme_minimal(); align_plot
ggsave(paste0(procimg, "/00_USGS_depth-tide_align-check.png"), align_plot,
       width = 8, height = 10)
message("Alignment complete. Check Processing Snapshots for file:
  00_USGS_depth-tide_align-check.png")
##------------------------------------------------------------------------------#
## 3. IMAGE GROUPING & LOOP SETUP ----
##------------------------------------------------------------------------------#
# Pull datetime from filename to group by time (images taken <5 min apart)
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
scene_df <- scene_df[-c(1:3),]
# Calculate time gaps & groups
scene_groups <- scene_df %>%
  group_by(group_id) %>%
  mutate(group_start = min(acq_time), 
         scene_count = n()
  ) %>%
  ungroup()
# View summary by group
summary_table <- scene_groups %>%
  group_by(group_id, group_start) %>%
  summarize(scenes_in_cluster = n(), .groups = "drop")
message(paste("Satellite Imagery grouped by aquisition date: \n > Images, total:",
              nrow(scene_groups),
              "\n > N.Groups:", nrow(summary_table))); print(summary_table)
# > List the PSS scenes by datetime grouping, to be mosaicked
group_list <- split(scene_df, scene_df$group_id)
# > Subset list to dates of interest (e.g., for image quality or by season)
focal_groups <- group_list[c(8:11,15:18)]
##------------------------------------------------------------------------------#
## 5. WATER DEPTH PROCESSING PIPELINE ----
##------------------------------------------------------------------------------#
message(paste0(">> Beginning Water Depth Processing pipeline: ",
               format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "..."))
inittime <- Sys.time()
for (group in focal_groups) {
  ##----------------------------------------------------------------------------##
  ## PARAMETERIZING ITERATION LABELS ----
  ##----------------------------------------------------------------------------##
  curr_id  <- unique(group$group_id) #unique group ID (based on PSS datetime)
  file_timestamp  <- group$time_string[1] #time stamp to use as ID (1st in grp)
  file_ts <- substr(file_timestamp, 1, nchar(file_timestamp) - 2) #only HHmm
  stackdate <- as_date(substr(file_timestamp, 1, nchar(file_timestamp) - 7)) #YMD
  mid_time <- mean(group$acq_time) #average of all PSScene time stamps
  start_d <- format(mid_time - days(1), "%Y-%m-%d") #format dates for retrieval
  end_d   <- format(mid_time + days(1), "%Y-%m-%d") #add 1-day buffer per side
  starttime <- Sys.time()
  message(paste0(">> STARTING: WATER DEPTH PROCESSING, PSS ", 
                 file_ts, " | ", 
                 format(Sys.time(), "%T %Z")))
  ##----------------------------------------------------------------------------##
  ## LOAD PROCESSED PSS IMAGERY  ----
  ## >> Point to the existing pre-processed imagery file locations
  ##----------------------------------------------------------------------------##
  mosaic_filepath    <- paste0(procrst, "/", file_ts, 
                               "_pss8b_mcf-norm01_normalized.tif")
  cleaned_filepath   <- paste0(procrst, "/", file_ts,
                               "_pss8b_mcf-norm02_boats-rm.tif")
  deglinted_filepath <- paste0(procrst, "/", file_ts,
                               "_pss8b_mcf-norm03_deglint.tif")
  # > Using the cleaned mosaic imagery to calculate predicted depth
  img_clean   <- rast(cleaned_filepath)
  set.crs(img_clean, crs(proj_crs))
  img_focal <- img_clean
  # img_focal   <- crop(img_clean, water_poly, mask=TRUE)
  # > Alternately can use deglinted mosaic imagery to calculate predicted depth
  # img_deglint <- rast(deglinted_filepath)
  # set.crs(img_deglint, crs(proj_crs))
  # img_focal <- img_deglint
  # img_focal <- crop(img_deglint, water_poly, mask=TRUE)
  plot(img_focal[[6]])
  ##----------------------------------------------------------------------------##
  ## TIDE FETCHING & SURFACE INTERPOLATION ----
  ##----------------------------------------------------------------------------##
  # LOGIC: Pull tide elevation data to correct water column depth
  message(paste("   >> Interpolating tidal elevation...\n      ",
                format(Sys.time(), "%T %Z")))
    # Fetch tide elevation from USGS Tidal Data (NAVD88), recorded every 6-15 min
  raw_tide <- dataRetrieval::read_waterdata_continuous( #pull tide info by time
    monitoring_location_id = site_ids, 
    parameter_code = "72279", #tidal elevation (ft)
    time = c(start_d, end_d))
  # Average records within window of time for a single value per station location
  tide_summary <- raw_tide %>%
    filter(time >= (mid_time - 600) & time <= (mid_time + 600)) %>%
    group_by(monitoring_location_id) %>%
    summarize(z_meters = mean(value, na.rm=TRUE) * 0.3048) #convert ft to meters
  # Process USGS station data & join spatial coordinates to tide info
  stations_df <- stations_sf %>% 
    inner_join(tide_summary, by = "monitoring_location_id") %>%
    # mutate(X = st_coordinates(.)[,1], Y = st_coordinates(.)[,2]) %>%
    st_drop_geometry() #convert back to a standard data frame
  # Prepare interpolation variables for given datetime across the bay
  tide_grid <- rast(ext(img_focal), res = 100, crs = proj_crs) 
  grid_coords <- crds(tide_grid, df = TRUE) #create df of grid cell coordinates
  colnames(grid_coords) <- c("X", "Y") #rename to match stations_df exactly
  fit_idw <- gstat(formula = z_meters ~ 1,  #define Inverse Distance Weight model
                   locations = ~X+Y, 
                   data = stations_df, 
                   set = list(idp = 2.0))
  tide_pred <- predict(fit_idw, grid_coords) #predict tide ht across grid points
  tide_surf <- tide_grid #put predicted values back into a raster
  values(tide_surf) <- tide_pred$var1.pred
  # Resample to match mosaic res & then mask to water area
  tide_surf_res <- resample(tide_surf, img_focal, method = "bilinear") %>%
    mask(img_focal[[1]])
  set.crs(tide_surf_res, crs(proj_crs)) #EXPLICITLY set CRS to match
  plot(tide_surf_res, 
       main=paste0("Instantaneous Tidal Surface (m): ", 
                   file_ts, " UTC"))
  ##----------------------------------------------------------------------------##
  ## WATER DEPTH CALIBRATION (TEMPORALLY-CALIBRATED BATHYMETRY = TCB) ----
  ## >> Calculate depth relative to the water surface at time of capture
  ##----------------------------------------------------------------------------##
  message("   >> Calibrating water depth by time & spectral index...\n      ",
          format(Sys.time(), "%T %Z"))
  # Calculate Observed Depth (from EAARL-B topobathy lidar, 2012)
  #  > Instantaneous Depth (m) = Tide Elevation (m) - Benthic Elevation (m)
  surf_tide_aligned <- project(tide_surf_res, 
                               img_focal[[1]], method = "bilinear")
  plot(surf_tide_aligned)
  print(paste("Extents match?",
              ext(surf_tide_aligned) == ext(img_focal[[1]])))
  print(paste("Resolutions match?",
              res(surf_tide_aligned)[1] == res(img_focal[[1]])[1]))
  # Combine topobathy with tide data & calculate instantaneous water depth
  topobathy <- topob_dem
  set.crs(topobathy, crs(proj_crs)) #ensure topob & tide surface have same CRS
  obs_topo_aligned <- extend(crop(topobathy, surf_tide_aligned),
                             tide_surf_res)
  print(paste("Extents match?",
              ext(surf_tide_aligned) == ext(obs_topo_aligned)))
  print(paste("Resolutions match?",
              res(surf_tide_aligned)[1] == res(obs_topo_aligned)[1]))
  # Force the tide extent to perfectly mirror the topo extent
  ext(surf_tide_aligned) <- ext(obs_topo_aligned)
  observed_depth <- surf_tide_aligned - obs_topo_aligned  
  plot(observed_depth, 
       main = paste0("Tide-Adjusted Depth (m): ", file_ts, " UTC"))
  # Calculate Predicted Depth (from PSS mosaicked imagery, 2023)
  #  > Use a spectral ratio (Stumpf Log-Ratio Algorithm) to model depth:
  #    pseudo-Satellite Derived Bathymetry (pSDB) = ln(Blue) / ln(Green)
  n <- 1000 #add a constant multiplier to avoid log(zero) or negative issues
  psdb <- log(n * img_focal[[2]]) / log(n * img_focal[[4]])
  ext(psdb) <- ext(observed_depth)
  plot(psdb, main = paste0("pseudo-Satellite Derived Bathymetry (pSDB): ", 
                           file_ts, " UTC"))
  # Build Linear Regression to model water depth from observed AND psdb data
  #  > Create dataframe of random pts where there is BOTH lidar & satellite data
  depth_df <- spatSample(c(psdb, observed_depth), size = 5000,
                         na.rm = TRUE, method = "random", values = TRUE)
  colnames(depth_df) <- c("Predicted", "Observed")
  cor_val <- cor(depth_df$Predicted, depth_df$Observed) #check correlation
  # Fit the linear model: Depth = m1 * (Ratio) + m0
  tcb_model <- lm(Observed ~ Predicted, data = depth_df)
  # summary(tcb_model)
  # Generate final temporally calibrated depth layer
  #  > Use the regression coefficients to predict depth everywhere
  m1 <- coef(tcb_model)[2] #slope = attenuation accuracy (1:1?)
  m0 <- coef(tcb_model)[1] #intercept = datum offset (how far off from 0?)
  tcb_m <-(m1 * psdb) + m0 #apply model to full mosaic
  plot(tcb_m)
  # Calculate Residual Error for Depth: (pSDB PRED - Tide-Adjusted Lidar OBS)
  #  > Positive = PRED thinks it's deeper than Observed (darker = possible SAV)
  #  > Negative = PRED thinks it's shallower than Observed (lighter = turbidity)
  tcb_residuals <- tcb_m - observed_depth   #where do pixels agree vs disagree?
  names(tcb_residuals) <- "ResidualError_m"
  ##----------------------------------------------------------------------------##
  ## DIAGNOSTICS: Calculate Summary Statistics ----
  ##----------------------------------------------------------------------------##
  # Linear Model & Residual Error measures
  r_sq   <- summary(tcb_model)$r.squared #how well do vars explain variability
  rmse_m <- sqrt(mean(values(tcb_residuals)^2, #sensitive to outliers
                      na.rm = TRUE)) 
  mae_m  <- mean(abs(values(tcb_residuals)), #robust average measure of model
                 na.rm = TRUE)
  # Average Tide Elevation Adjustment across extent
  band_avg_tide <- global(surf_tide_aligned, fun = "mean", na.rm = TRUE)
  avg_tide_elev <- as.vector(unlist(band_avg_tide))
  band_sd_tide <- global(surf_tide_aligned, fun = "sd", na.rm = TRUE)
  sd_tide_elev <- as.vector(unlist(band_sd_tide))
  # Relative Residual Error: ((Predicted - Observed) / Observed) * 100
  #  > Make RE relative in order to be able to compare variation across time 
  #  > Mask shallow depths <10cm (0.1m) to prevent errors along shoreline
  observed_safe <- ifel(observed_depth < 0.1, NA, observed_depth)
  relative_error_pct <- ((tcb_m - observed_safe) / observed_safe) * 100
  band_avg_re <- global(relative_error_pct, fun = "mean", na.rm = TRUE)
  avg_rel_err <- as.vector(unlist(band_avg_re))
  # print(avg_rel_err)
  resid_avg <- global(tcb_residuals, fun = "mean", na.rm = TRUE)
  resid_avg <- as.vector(unlist(resid_avg))
  resid_sd <- global(tcb_residuals, fun = "sd", na.rm = TRUE)
  resid_sd <- as.vector(unlist(resid_sd))
  # Create a log entry for given group
  log_entry <- data.frame(
    group_id = curr_id,       
    file_ts = file_ts,        
    tide_avg = avg_tide_elev,
    tide_sd = sd_tide_elev,
    tcb_rsquared = r_sq,      #higher values indicate better model performance
    slope = m1,               
    intercept = m0,           
    RMSE = rmse_m,            #lower values indicate better model performance
    MAE = mae_m,              #lower values indicate better model performance
    corr = cor_val,
    RE_avg = resid_avg,
    RE_sd = resid_sd
  )
  # Write/Append log entry to CSV file
  log_file <- file.path(proccsv, "TCB_processing_log.csv")
  write.table(log_entry, 
              log_file,
              sep = ",", col.names = !file.exists(log_file), 
              row.names = FALSE, append = file.exists(log_file))
  message(paste("--- Validation Report ---", 
                "\nCorrelation between SDB & Tidal-Adjusted LiDAR Depths: ",
                round(cor_val, 3),
                "\nR-Squared: ", round(r_sq, 3),
                "\nRMSE:      ", round(rmse_m, 3), "meters",
                "\nMAE:       ", round(mae_m, 3), "meters", 
                "\n   >> Completed Water Depth Calibration at", 
                format(Sys.time(), "%T %Z")))
  # EXPORT OUTPUTS (ALL DEPTHS)
  # > Stack the depth layers into a single multiband raster
  completed_depth_stack <- c(
    tcb_m,             # Band 1: Temp-Calibrated Bathymetry (m) = final depth map
    surf_tide_aligned, # Band 2: Instantaneous Tide (NAVD88) = vertical anchor
    observed_depth,    # Band 3: Tide-Adjusted Depth (m) = tide elev - topob lidar
    tcb_residuals      # Band 4: TCB Residuals = error map between pSDB & TAD
  )
  names(completed_depth_stack) <- c("TCB_Derived_Depth_m",
                                   "Tide_NAVD88_m", 
                                   "Tide_Adj_Depth_m",
                                   "Residuals_m")
  # > Evaluate alignment with imagery-based raster stacks
  if (compareGeom(img_focal, completed_depth_stack, stopOnError = FALSE)) {
    message("Geometries match perfectly. Skipping resample.")
    final_depth_stack <- completed_depth_stack
  } else {
    message("Geometries do not match. Resampling depth to imagery stack...")
    # > If extents do not match, adjust alignment on depth stack
    final_depth_stack <- terra::resample(completed_depth_stack, img_focal, 
                                    method = "bilinear")
  }
  # Write Temporally Calibrated Bathymetry (TCB) final raster stack (8-band)
  out_filename <- paste0("TCB_allDepths_", file_ts, ".tif") 
  writeRaster(final_depth_stack, file.path(procdep, out_filename), 
              overwrite = TRUE)
  message(">> FINISHED Satellite Depth Pre-Processing for PSS Timestamp: ",
          stackdate, " \n   on ",
          format(Sys.time(), "%A, %b %d %Y %I:%M %p"))
  ##----------------------------------------------------------------------------##
  ## > VISUALS & DIAGNOSTICS
  ##----------------------------------------------------------------------------##
  message("   >> Loading calibrated depth layers...\n      ",
          format(Sys.time(), "%T %Z"))
  # LOAD DEPTH DATA: Read in water depth raster stack (built during pre-proc)
  # final_depth_stack <- rast(file.path(procwcc,
  #                                     paste0("TCB_allDepths", file_ts, ".tif")))
  tcb_final <- final_depth_stack[[1]]   #temporally-calibrated bathymetry map
  tide_map <-  final_depth_stack[[2]]   #interpolated tidal elev map
  tad_topo <-  final_depth_stack[[3]]   #tide-adjusted depth
  residuals <- final_depth_stack[[4]]   #residual error of water depth
  # plot(final_depth_stack)
  # title(main = paste0("BBLEH Depth Stack: ",
  #                     stackdate), line = 11)
  ##----------------------------------------------------------------------------##
  ## DEPTH ANOMALY PROCESSING: SPATIAL RESIDUAL CLUSTERING ----
  ##----------------------------------------------------------------------------##
  # METHOD: Local Moran's I (LISA) on TCB Residuals
  # LOGIC: Use depth errors for Turbidity & Seagrass Detection
  message("   >> Analyzing depth residuals...\n      ",
          format(Sys.time(), "%T %Z"))
  # AGGREGATE: Coarsen to ~30m pixels for broad pattern detection
  res_coarse <- aggregate(residuals, fact = 10, fun = "mean", na.rm = TRUE)
  # STANDARDIZE RESIDUALS: calculate the Z-scores for your coarse residuals
  res_coarse_z <- scale(res_coarse)
  plot(res_coarse_z)
  #  > Perfectly standardized Z-score raster will always have mean = 0 & sd =1
  print(global(res_coarse_z, fun = c("mean", "sd"), na.rm = TRUE))
  # CONVERT TO VECTOR: Raster into polygons to define "neighbors"
  res_poly <- as.polygons(res_coarse_z, values = TRUE, na.rm = TRUE) %>% 
    st_as_sf()
  colnames(res_poly)[1] <- "resid_val"
  plot(res_poly)
  # DEFINE NEIGHBORS: Queen's Contiguity = sharing a border or point
  nb <- poly2nb(res_poly) #construct neighbors list from polygons
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE) #spatial weights for nbs
  # CALCULATE LOCAL MORAN'S I: High error surrounded by high error = clustering
  loc_m <- localmoran(res_poly$resid_val, lw, zero.policy = TRUE)
  res_poly$Ii <- loc_m[,1] #Local Moran's Index
  res_poly$Pr <- loc_m[,5] #P-value (significance)
  # CLASSIFY CLUSTERS: Group based on Residual Error value + LocM Index
  #  > High-High = Turbidity Plume or Seagrass (consistently shallower than truth)
  #  > Low-Low   = Deep Bias (consistently deeper than truth)
  res_poly <- res_poly %>%
    mutate(
      type = case_when(
        Pr > 0.05 ~ "Insignificant (p > 0.05)",
        resid_val > 0 & Ii > 1 ~ "1. Turbidity (High-High)",
        resid_val > 0 & Ii < 1 ~ "2. Shallow Bias (High-Low)\n    *Potential SAV",
        resid_val < 0 & Ii < 1 ~ "4. Deep Bias (Low-Low)",
        resid_val < 0 & Ii > 1 ~ "5. Low-High (Uncertain)",
        TRUE ~ "3. Non-Anomalous (Expected)"
      )
    )
  plot(res_poly)
  # DIAGNOSTIC PLOTS: Evaluate depth processing pipeline
  p_depth_der <-ggplot() +
    geom_sf(data = water_poly, fill = "grey90", color = "grey60") +
    geom_spatraster(data = tcb_final, maxcell = 5e5) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    scale_fill_whitebox_c(palette = "deep", direction = 1,
                          na.value = "transparent",
                          name = "Derived \nDepth \n(m)") +
    labs(title = "Temporally-Calibrated Bathymetry (m)",
         subtitle = "Modeled from observed & predicted depths",
         caption = sprintf("R-squared: %.2f  |  RMSE: %.2f m", r_sq, rmse_m)
    ) +
    theme_void()
  p_tide <- ggplot() + 
    geom_sf(data = water_poly, fill = "grey90", color = "grey60") +
    geom_spatraster(data = tide_map, maxcell = 5e5) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    scale_fill_whitebox_c(palette = "muted",
                          name = "Tidal \nElev \n(m)") +
    labs(title = "Instantaneous Tide Surface (NAVD88)",
         subtitle = 
           "Interpolated across USGS stations (+/-10min)",
         caption = paste0("Avg. Tidal Adjustment: ", 
                          round(avg_tide_elev, 2), "m (SD = ", 
                          round(sd_tide_elev, 2), ")"
         )
    ) +
    theme_void()
  p_reserr <- ggplot() +
    geom_sf(data = water_poly, fill = "grey90", color = "grey60") +
    geom_spatraster(data = residuals, maxcell = 5e5) +
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    scale_fill_gradient2(
      low = "#0571b0",      #blue: underestimating
      mid = "lightyellow",  #light: accurate
      high = "#ca0020",     #red: overestimating
      midpoint = 0,
      na.value = "transparent",
      name = "Residual \nError \n(m)"
    ) +
    labs(title = "Spatial Bias in Derived Depth",
         subtitle = "Red = Predicted too deep | Blue = Predicted too shallow",
    ) +
    theme_void()
  p_depth_error <- ggplot() +
    geom_sf(data = water_poly, fill = "grey90", color = "grey60") +
    geom_spatraster(data = relative_error_pct, maxcell = 5e5) +
    scale_fill_gradient2(low = "dodgerblue4", 
                         mid = "lightyellow", 
                         high = "firebrick", 
                         midpoint = 0, 
                         limits = c(-100, 100), 
                         oob = squish, 
                         name = "Relative\nError \n(%)",
                         na.value = "transparent") + 
    geom_sf(data = water_poly, color = "grey60", linewidth = 0.4,
            fill = "transparent") +
    labs(title = "Spatial Bias in Derived Depth Error",
         subtitle =
           "Predicted (pSDB) vs. Observed (Lidar)",
         caption = paste0("Avg. Residual Error: ", 
                          round(resid_avg, 2), "m (SD = ", 
                          round(resid_sd, 2), ")"
         )
    ) +
    theme_void()
  #  >  Join plots and include relevant statistics
  diagnostics_depth <- p_depth_der + p_tide + p_depth_error +
    # plot_layout(guides = 'collect') +
    plot_annotation(title = paste0(
      "SAV Diagnostics 02  |  Depth Anomaly Calculations"),
      subtitle = paste0(
        "Barnegat Bay SAV Mapping, PlanetScope Satellite Imagery for ", 
        as_date(stackdate))
# caption = "Positive Error (red) = predicted deeper, potential SAV with dark signal |  Negative Error (blue) = predicted shallower, turbidity masking deeper water"
    ) +
    theme_void(); diagnostics_depth
  # EXPORT GRAPHIC: Summary diagnostics
  ggsave(file.path(procimg, 
                   paste0("02b_SAV-diagnostics_depth-calibration_allDepths_",   
                          file_ts, ".png")), diagnostics_depth,
         width = 11, height = 9, dpi = 300)
  # MODEL VALIDATION
  #  > Create the text label
  plot_label <- sprintf("R-squared: %.2f\nRMSE: %.2f m", r_sq, rmse_m)
  # BUILD HEXBIN PLOT
  p_model_validation <- ggplot(
    depth_df, aes(x = Observed, y = Predicted)) +
    # A. The Density Hexagons (Solves the overplotting problem)
    geom_hex(bins = 100) +
    scale_fill_viridis_c(option = "plasma", name = "Pixel Count", 
                         trans = "log10", labels = comma) +
    # B. The 1:1 Reference Line (The "Perfect Model" baseline)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray30", linewidth = 1) +
    # C. Your Actual Regression Line (How your model actually performed)
    geom_smooth(method = "lm", color = "firebrick", 
                linetype = "solid", se = FALSE, linewidth = 1.2) +
    # D. Add the Metrics to the top left corner
    annotate("text", x = max(depth_df$Observed, na.rm=TRUE), 
             y = min(depth_df$Predicted, na.rm=TRUE), 
             label = plot_label, hjust = 1, vjust = 0.5, 
             fontface = "bold", size = 5) +
    # E. Formatting for Publication
    theme_minimal() +
    labs(title = "SUPP Diagnostics 02  |  Depth Model Validation", 
         subtitle = 
           "Predicted (pSDB) vs. Observed (Tide-Adjusted Topobathy Lidar)",
         x = "Observed TCB Depth (m)",
         y = "Predicted Satellite Depth (m)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          axis.title = element_text(face = "bold", size = 12),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  p_model_validation
  # EXPORT GRAPHIC: Depth clusters map as PNG
  ggsave(file.path(procimg, 
                  paste0("S02_SUPP-diagnostics_depth-model-validation_allDepths_",
                          file_ts, ".png")), p_model_validation,
         width = 9, height = 9, dpi = 300)
  # MAP THE ERROR CLUSTERS
  #  > Error Clusters Map
  p_residuals <- ggplot() +
    geom_sf(data = water_poly, color = "grey50", fill="grey90") +
    geom_sf(data = res_poly, aes(fill = type), color = NA) +
    scale_fill_manual(values = c(
      "1. Turbidity (High-High)" = "orchid2",
      "2. Shallow Bias (High-Low)\n    *Potential SAV" = "steelblue3",
      "3. Non-Anomalous (Expected)" = "honeydew2",
      "4. Deep Bias (Low-Low)" = "sienna",
      "5. Low-High (Uncertain)" = "gold",
      "Insignificant (p > 0.05)" = "navajowhite3"
    )) +
    geom_sf(data = water_poly, color = "grey30", linewidth = 0.4,
            fill = "transparent") +
    labs(title =  
           "SUPP Diagnostics 03  |  Depth Anomaly-Based Clustering", 
         subtitle = paste0("Turbidity Error Types for ", stackdate),
         fill = "Error\nType") +
    theme_minimal(); print(p_residuals)
  # EXPORT GRAPHIC: Depth clusters map as PNG
  ggsave(file.path(procimg, 
                   paste0("S03_SUPP-diagnostics_depth-error-clusters_allDepths_",
                          file_ts, ".png")), p_residuals,
         width = 9, height = 16, dpi = 300)
}
  ##----------------------------------------------------------------------------##
  ## BUILD DEPTH ANOMALY POLYGONS ----
  ##----------------------------------------------------------------------------##
  # message("   >> Exporting depth anomaly polygons...\n      ",
  #          format(Sys.time(), "%T %Z"))
  # # LOGIC: Positive Depth Residual (Darker than expected) AND Clustered (Ii > 0)
  # #  > These are the "SAV Bed Detection" arguments
  # turbidity_clusters <- res_poly %>%
  #   filter(resid_val > 0.5) %>%  #threshold: Residual > 0.5m (significant error)
  #   filter(Ii > 0) %>%              #clustered (high surrounded by high)
  #   filter(Pr < 0.1)           #statistically significant
  # seagrass_clusters <- res_poly %>%
  #   filter(resid_val < 0) %>%  #threshold: Residual > 0.5m (significant error)
  #   filter(Ii > 0) %>%              #clustered (high surrounded by high)
  #   filter(Pr < 0.1)           #statistically significant
  # # LOGIC: Negative Depth Residual (Brighter than expected) AND Clustered (Ii>0)
  # #  > This can locate deep water bias and possible turbidity
  # spatial_outliers <- res_poly %>%
  #   filter(resid_val < -0.5) %>% #threshold: Residual <-0.5m (significant error)
  #   filter(Ii < 0) %>%              #clustered (high surrounded by high)
  #   filter(Pr < 0.1)               #clustered (High surrounded by High)
  # # VISUALIZE RESULTS: Depth Residuals Map
  # ggplot() +
  #   geom_sf(data = water_poly, color = "transparent", fill = "cornflowerblue") +
  #   # geom_sf(data = turbidity_clusters, aes(fill = resid_val),
  #   #         fill = "navy", color = "transparent", linewidth = 0.4,
  #   # )  +
  #   geom_sf(data = seagrass_clusters, aes(fill = resid_val),
  #           fill = "olivedrab2", color = "olivedrab2"
  #   )  +
  #   labs(title = "TCB-derived Errors: Potential SAV & Deep Water",
  #        subtitle = paste("Depths Shallower & Deeper than Expected for",
  #                         file_ts)) +
  #   theme_minimal()
  # # EXPORT OBJECTS: Benthic cluster polygons as SHPs: Seagrass & Deep Water
  # sf::st_write(seagrass_clusters, 
  #          file.path(paste0(
  #            procben, "/anomaly-clustering_shallow_SAV_", 
  #            file_ts, ".shp")),
  #          append = FALSE)
  # # sf::st_write(spatial_outliers, 
  # #          file.path(paste0(
  # #            procben, "/anomaly-clustering_xdeep_Turbid_", 
  # #            file_ts, ".shp")),
  # #          append = FALSE)
  ##----------------------------------------------------------------------------##
  # # CLEAN UP: Remove temporary objects & trigger garbage collection
  # #  > Leave only layers needed for Algae & Turbidity modeling
  # rm(tide_grid, tide_pred, tide_surf, topobathy, observed_depth, psdb,
  #    tcb_model, tcb_m, res_poly, res_coarse, df_points, tide_surf_res,
  #    tcb_residuals, seagrass_clusters, img_focal, grid_coords) 
  # gc()
##----------------------------------------------------------------------------
  endtime <- Sys.time()
  # duration <- endtime - starttime
  loop_time <- round(as.numeric(difftime(endtime, starttime, units = "mins")), 2)
  message(paste0(">> FINISHED: DEPTH & WATER QUALITY PROCESSING, PSS ", 
                 file_ts, " | ", 
                 format(Sys.time(), "%T %Z"), "\n   Duration: ", loop_time,
                 " minutes", "\n  ----------------------------------------"))
  ##----------------------------------------------------------------------------##
  beep(2) #----
  ##----------------------------------------------------------------------------##
elapsed <- round(as.numeric(difftime(endtime, inittime, units = "mins")), 2)
message(paste0(">> Depth & Water Quality Processing completed for all dates on
    ", format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "  |  Elapsed Time: ", 
               elapsed, " minutes"))
Sys.sleep(2)
beep(4)
# ============================================================================= #
