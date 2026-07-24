# ============================================================================= #
# Script Name:  05_clean_pred-vars.R
# Description:  Pulls together raster stacks to build daily modeling input. Clips 
#               the extent to area of interest (based on water depth from topob)
# Author:       Jess M. Stitt
# Date:         2026-03-31
# ============================================================================= #

# ============================================================================= #
# ---- TASK DESCRIPTION ----
# Input(s):
#  > [CSV] Cleaned Field Reference Dataset 
#     (A subset of plot data, to be used in training models)
#  > [SHP] SAV Polygon
#     (Vector)    
#  > [TIF] Processed Predictor Variable Raster Stacks
#     (Separate Imagery, Depth, and Water Quality raster stacks to join)
# Output(s): 
#  > [, ] 
#     ()
# ============================================================================= #
# ---- LIBRARIES SPECIFIC TO THIS TASK ----
#  <NONE>
# ============================================================================= #
# ---- PROCESSING STEPS ----
##------------------------------------------------------------------------------#
## 1. 
##------------------------------------------------------------------------------#
pss_stack_paths <- list.files(procwcc,
                              pattern = "^wcc-stack_(.*)\\.tif$",
                              full.names = TRUE)
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
  ## A. PARAMETERIZE ITERATION LABELS ----
  ##----------------------------------------------------------------------------##
  file_ts <- value #only HHmm 
  file_timestamp <- value %>%
    as.POSIXct(format="%Y%m%d_%H%M", tz="UTC")
  stackdate <- as_date(substr(value, 1, nchar(value) - 5)) #YMD
  starttime <- Sys.time()
  message(paste0(">> STARTING: WATER COLUMN PROCESSING, PSS ", 
                 file_ts, " | ", 
                 format(Sys.time(), "%T %Z")))
  ##----------------------------------------------------------------------------##
  ## B. DEFINE PREDICTOR VARIABLES ----
  ##----------------------------------------------------------------------------##
  # > Point to existing pre-processed imagery file locations
  full_stack <- paste0(procstk, "/", file_ts, "_pss8b_rf-full-stack_2m.tif")
  if (file.exists(full_stack)) {
    # LOAD PREVIOUSLY PROCESSED RASTER STACK
    print(paste0(">> Existing file found; loading raster stack for ",
                 stackdate))
    predictors2m <- rast(full_stack)
    predictors <- predictors2m[[c(3:5,11:13,17:29)]]
  } else {
    # PROCESS RAW IMAGE TILES: load, mosaic, reproject 8-band imagery
    print(paste0(">> File not found; creating raster stack for ",
                 stackdate))
    # Read in the raster stacks for all bands plus depth stack
    rfs_clean <- rast(file.path(procrst, paste0(
      file_ts, "_pss8b_mcf-norm02_boats-rm.tif")))
    rfs_glint <- rast(file.path(procrst, paste0(
      file_ts, "_pss8b_mcf-norm03_deglint.tif")))
    names(rfs_glint) <- paste0(names(rfs_clean), "_dg")
    rfs_depth <- rast(file.path(procdep, paste0("TCB_allDepths_",
                                                file_ts, ".tif")))
    rfs_wqual <- rast(file.path(procwcc, paste0("wcc-stack_allDepths_",
                                                file_ts, ".tif")))
    if (compareGeom(rfs_wqual, rfs_depth, stopOnError = FALSE)) {
      message("Geometries match perfectly. Skipping resample.")
      # > Combine into comprehensive raster stack (ALL depths)
      predictor_stack <- c(
        rfs_clean,             #select all 8 bands: boats removed
        rfs_glint,             #select all 8 bands: glint-corrected
        rfs_depth,             #select all 4 bands: depth-based
        rfs_wqual)             #select all 9 bands: water column correction
    } else {
      message("Geometries do not match. Resampling depth to imagery stack...")
      # > If extents do not match, adjust alignment on depth stack
      #   *NOTE: Using method = "bilinear" for Continuous data (near for Cat)
      rfs_depth_aligned <- resample(rfs_depth, rfs_wqual, method = "bilinear")
      writeRaster(rfs_depth_aligned,
                  file.path(procdep,
                            "TCB_allDepths_aligned",
                            file_ts, ".tif"),
                  gdal = c("COMPRESS=LZW", #zips the data to save space
                           "TFW=YES",      #prevents crash at 4GB
                           "BIGTIFF=YES"), #generates sidecar text file
                  overwrite = TRUE)
      # > Combine into comprehensive raster stack (ALL depths)
      predictor_stack <- c(
        rfs_clean,             #select all 8 bands: boats removed
        rfs_glint,             #select all 8 bands: glint-corrected
        rfs_wqual,             #select all 9 bands: water column correction
        rfs_depth_aligned)     #select all 4 bands: depth-based (aligned)
    }
    # > Create new raster stack with all layers
    writeRaster(predictor_stack,
                file.path(procstk,
                          paste0(file_ts, "_pss8b_rf-full-stack_allDepths.tif")),
                gdal = c("COMPRESS=LZW", #zips the data to save space
                         "TFW=YES",      #prevents crash at 4GB
                         "BIGTIFF=YES"), #generates sidecar text file
                overwrite = TRUE)
    # > Clip the layers to only suiatble SAV habitat (0-2m depths)
    predictors2m <- mask(predictor_stack, sav_poly)
    # > Create new raster stack with all layers
    writeRaster(predictors2m,
                file.path(procstk,
                          paste0(file_ts, "_pss8b_rf-full-stack_2m.tif")),
                gdal = c("COMPRESS=LZW", #zips the data to save space
                         "TFW=YES",      #prevents crash at 4GB
                         "BIGTIFF=YES"), #generates sidecar text file
                overwrite = TRUE)
    # > Remove any layers not of interest for modeling
    predictors <- predictors2m[[c(3:5,11:13,17:29)]] #omit norm+dg spectral bands 
    # names(predictors)
  }
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
## 2.  APPEND TRAINING PREDICTORS ----
## >>  Description:  Extracts spatial predictor variables at survey locations and 
##     binds them directly to the training dataframe.
##----------------------------------------------------------------------------##
# 1. LOAD SPATIAL DATA
##----------------------------------------------------------------------------##
print("Loading points and raster data...")
#  > Pull in training plot coordinates as vector
train_vect <- vect(train_plots)
#  > Load both seasonal composites
sp_comp <- rast(paste0(procagg, "/pss8b_SP_composite-stack_2m-clean.tif"))
fa_comp <- rast(paste0(procagg, "/pss8b_FA_composite-stack_2m-clean.tif"))
#  > Load the static Bathymetry raster
bathy_raster <- topob_dem
bathy.col <- colorRampPalette(c("darkblue", "aquamarine2"))
plot(bathy_raster, col=c(bathy.col(20)))
plot(train_vect, add=TRUE, col="orangered", pch=16, cex=0.4)
##----------------------------------------------------------------------------##
# 2. EXTRACT DYNAMIC SEASONAL VARIABLES
##----------------------------------------------------------------------------##
print("Extracting PlanetScope SPRING signatures...")
sp_extracted <- extract(sp_comp, train_vect, bind = FALSE)
colnames(sp_extracted)[-1] <- paste0("SP_", colnames(sp_extracted)[-1])

print("Extracting PlanetScope FALL signatures...")
fa_extracted <- extract(fa_comp, train_vect, bind = FALSE)
colnames(fa_extracted)[-1] <- paste0("FA_", colnames(fa_extracted)[-1])
##----------------------------------------------------------------------------##
# 3. EXTRACT STATIC ENVIRONMENTAL DATA
##----------------------------------------------------------------------------##
print("Extracting Bathymetry...")
bathy_extracted <- extract(bathy_raster, train_vect, bind = FALSE) #depths, EAARL
# Rename the 2nd column (the actual data column) so it looks clean in the spreadsheet
colnames(bathy_extracted)[2] <- "Bathy_eaarlb"
##----------------------------------------------------------------------------##
# 4. MERGE EVERYTHING TOGETHER
##----------------------------------------------------------------------------##
print("Merging datasets...")
train_df <- as.data.frame(train_vect)
train_df$ID <- 1:nrow(train_df)
# Use dplyr to staple all three extractions to original plot data
final_training_data <- train_df %>%
  left_join(bathy_extracted, by = "ID") %>% #can add addl static vars same way
  left_join(sp_extracted, by = "ID") %>%
  left_join(fa_extracted, by = "ID") %>%
  select(-ID)
# Check columns to verify
print(colnames(final_training_data))
# head(final_training_data)
##----------------------------------------------------------------------------##
# 5. SAVE FOR MODELING
##----------------------------------------------------------------------------##
write.csv(final_training_data, 
          paste0(proccsv, "/bb23_sav_train70_multiseason_AGG.csv"), 
          row.names = FALSE)
print("Extraction complete. Training data is ready for Random Forest modeling.")
##----------------------------------------------------------------------------##
beep(10) #----
##----------------------------------------------------------------------------##
# ============================================================================= #