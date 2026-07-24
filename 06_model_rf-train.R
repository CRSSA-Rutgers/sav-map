# ==============================================================================
# Script Name:  06_model_rf-train.R
# Description:  Trains seasonal Random Forest models and generates daily 
#               seagrass classification maps for Barnegat Bay.
# Author:       Jess M. Stitt
# Date:         2026-03-31
# ==============================================================================

# ---- TASK DESCRIPTION ----
# Input(s):
#  > [] 
#     ()
##------------------------------------------------------------------------------#
# ---- LIBRARIES SPECIFIC TO THIS TASK ----
library(terra)
library(randomForest)
# 0.50 is standard. 0.65 means 65% of the RF trees must agree it is seagrass.
sav_threshold <- 0.7
# Output(s): 
#  > [, ] 
#     ()

# Output(s): 
#  > [, ] 
#     ()
# ============================================================================= #
# ---- USER-DEFINED VARIABLES FOR SCRIPT ----
# <NONE>
# ============================================================================= #
# ---- PROCESSING STEPS ----
##------------------------------------------------------------------------------#
## 1. AUTOMATED DATE SORTER ----
##------------------------------------------------------------------------------#
print("Scanning for processed PlanetScope daily imagery stacks...")
# Load path & pattern to match exactly where pre-processed raster stacks live
# > Scan folder for all relevant predicotr variable stack files
all_files <- list.files(procstk,
                        pattern = "rf-full-stack_2m.tif$",
                        full.names = TRUE)
# > Extract full datetimes from filenames
all_datetimes <- regmatches(basename(all_files),#extract date & time from file
                            regexpr("\\d{8}_\\d{4}", basename(all_files)))
# > Define seasonal split (e.g., August 1st)
season_split <- "0801"
# > Initialize empty named list
image_dates <- list(SP = c(), FA = c())
# > Loop through the found dates and sort
for (current_date in all_datetimes) {
  current_month <- substr(current_date, 5, 6)
  # --- AUGUST FILTER ---
  if (current_month == "08") next # Skip August
  # Extract the MMDD part for the seasonal split
  mmddhhmm <- substr(current_date, 5, 10)
  # Sort based on split logic
  if (mmddhhmm < season_split) {
    image_dates$SP <- c(image_dates$SP, current_date) #append to the Spring list
  } else {
    image_dates$FA <- c(image_dates$FA, current_date) #append to the Fall list
  }
}
message(paste0(
  ">> Successfully loaded and sorted the following datetimes:\n   SPRING: ", 
               image_dates [1], "\n   FALL: ", image_dates[2]))
##------------------------------------------------------------------------------#
## 2. LOAD STATIC ASSETS ----
##------------------------------------------------------------------------------#
# Load the training data subset built in script 01
full_train_data <- read.csv(paste0(
  proccsv, "/bb23_sav_train70_multiseason_AGG.csv"))
# str(full_train_data)
# Load the Bathymetry raster (we will need to staple this to the daily images!)
bathy_raster <- topob_dem
# ==============================================================================
# 3. OUTER LOOP: TRAIN SEASONAL MODELS
# ==============================================================================
inittime <- Sys.time()
message(paste0(">> Beginning Random Forest Modeling: ",
               format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "..."))
for (season in names(image_dates)) {
  starttime <- Sys.time()
  season_dates <- image_dates[[season]]
  if (length(season_dates) == 0) {
    print(paste("! No dates to process for", season, " - Skipping."))
    next
  }
  print(paste("========== STARTING SEASON:", season, "=========="))
  # --- A. Dynamic Data Filtering ---
  print("Isolating and formatting seasonal training data...")
  target_col <- paste0("Occ", season)
  season_train_data <- full_train_data %>%
    select(
      all_of(target_col),  # Dynamically grab the correct Occ column
      Bathy_eaarlb,        # Static environmental variable
      starts_with(paste0(season, "_")) # Dynamically grab SP_ or FA_ columns
    ) %>%
  rename(Presence = all_of(target_col)) # Rename OccSP/OccFA to 'Presence'
  # --- B. The Data Leakage Fix ---
  # Strip "SP_" or "FA_" prefix so the model just learns generic band names
  colnames(season_train_data) <- str_replace(
    colnames(season_train_data), 
    pattern = paste0("^", season, "_"), 
    replacement = ""
  )
  # --- C. Train the Random Forest Model ---
  print("Training Random Forest model (ntree = 500)...")
  rf_model <- randomForest(
    as.factor(Presence) ~ ., 
    data = season_train_data, 
    ntree = 500, 
    importance = TRUE,
    na.action = na.omit # Safety catch for any empty pixels
  )
  print("Model Training Complete! Commencing Daily Predictions...")
  # ============================================================================
  # 4. INNER LOOP: DAILY PREDICTIONS
  # ============================================================================
  for (current_date in season_dates) {
    stackdate <- as_date(substr(current_date, 1,8))
    print(paste("   >> Building SAV Predictions for", stackdate))
    # --- A. Load the Daily Image ---
    daily_stack <- rast(paste0(procstk, "/", current_date, 
                               "_pss8b_rf-full-stack_2m.tif"))
    # --- B. Add Bathymetry to the Stack ---
    # 1. Snap the static Bathymetry layer to the exact grid/extent of mosaic
    # Method = "bilinear" because depth is continuous (use "near" for categorical)
    bathy_matched <- resample(x = bathy_raster, y = daily_stack, 
                              method = "bilinear")
    names(bathy_matched) <- "Bathy_eaarlb"
    # CRITICAL: Because the RF model was trained with Bathymetry, it expects to 
    # see a layer with same name ("Bathy_eaarlb") when it makes prediction.
    predictor_stack <- c(daily_stack, bathy_matched)
    # --- C. Predict the Map (Probability Mode) ---
    occ_map <- terra::predict(
      object = predictor_stack, 
      model = rf_model, 
      type = "prob",
      na.rm = TRUE
    )
    # Save the Model Object
    saveRDS(occ_map, paste0(procres, "/RFmodel-pOCC_", current_date,".rds"))
    # To map the likelihood of seagrass, you specifically have to extract 
    # the "Present" column (usually column 2) to feed into terra.
    predicted_probs <- occ_map[[2]]
    plot(predicted_probs, main = paste0("SAV Probability for ", stackdate))
    Sys.sleep(2)
    # --- C.2. Apply Custom Threshold & Isolate SAV ---
    # The ifel() function acts exactly like an ifelse() statement for rasters:
    # If the probability is >= your threshold, write a 1. Otherwise, write NA.
    sav_only_map <- ifel(predicted_probs >= sav_threshold, 1, NA)
    plot(predicted_probs, main = paste0(
      "Predcited SAV for ", stackdate, 
      "\n >>Confidence Threshold: ", sav_threshold))
    # sav_only_map <- subst(predicted_map, from = 0, to = NA)
    # --- D. Save the Final Product ---
    out_name <- paste0(procres, "/rf-predict_conf", sav_threshold*100, "_SAV_",
                       current_date, ".tif")
    writeRaster(sav_only_map, out_name, overwrite = TRUE)
    message("   >> Completed Probability Map for ", stackdate, "...\n      ",
            format(Sys.time(), "%T %Z"))
    gc() #force memory cleanup after every heavy spatial prediction
  }
  message("   >> RF Classification Complete\n      ",
          format(Sys.time(), "%T %Z"))
  endtime <- Sys.time()
  loop_time <- round(as.numeric(difftime(endtime, starttime, units = "mins")), 2)
  message(paste0(">> FINISHED: RF MODELING ALL DATES FOR SEASON: ", season,  
                 format(Sys.time(), "%T %Z"), "\n   Duration: ", loop_time,
                 " minutes", "\n  ----------------------------------------"))
  beep(3)
  print(paste("========== FINISHED ALL DATES FOR", season, "=========="))
}
elapsed <- round(as.numeric(difftime(endtime, inittime, units = "mins")), 2)
message(paste0(">> Random Forest modeling completed for all dates on 
    ", format(Sys.time(), "%A, %b %d %Y %I:%M %p"), "  |  Elapsed Time: ", 
               elapsed, " minutes"))
print("PIPELINE COMPLETE.")
beep(4)
# ============================================================================= #



# ==============================================================================
# ---- INPUTS FOR SCRIPT ----
# 
rfstack_filepaths <- list.files(paste0(paste0(procstk, "/tempOLD")),
                                pattern = "rf-full-stack.tif$",
                                full.names = TRUE)
rf_filenames <- list(basename(rfstack_filepaths))
rf_datetimes <- substr(basename(rfstack_filepaths), 0, #isolate PSS timestamp
                       nchar(basename(rfstack_filepaths)) - 24)
rf_dates <- substr(basename(rfstack_filepaths), 0, #isolate PSS timestamp
                   nchar(basename(rfstack_filepaths)) - 29)
survey_dates <- unique(na.omit(rf_datetimes))
season_split <- as.Date(paste0("2023-08-01")) #delineation for SPRING vs FALL
dates_grid <- expand_grid(
  Date = survey_dates
)
target_classes <- c( #CLIP BY SEASON: SPRING (May-Jul) VS FALL (Aug-Oct)
  "Occ"    #Presence/Absence
  # "Class"  #SAV Density Classes
  # "Rupp",  #Ruppia-Only
  # "Zost",  #Zostera-Only
  # "DomCov"
)
# # Group predictor columns into a list of scenarios & subset by variables
# scenarios <- list(
#   "Optics_only" = predvars[c(1:8)],
#   # "Optics_and_Depth" = predvars[c(9:16,28,30:31)],
#   # "Depth_and_WQ" = predvars[c(17:18,27:28,30:31)],
#   "Sagawa_Depth_WQ" = predvars[c(17:28,30:31)],
#   # "Sagawa_only" = predvars[c(19:26)],
#   "Full" = predvars
# )
mod_grid <- expand_grid(
  Target = target_classes, #4 total per season
  # Scenario_Name = names(scenarios),
  Date = survey_dates #11 total per year (SP=4, FA=7)
)
# ====
# Initialize empty lists to store results of each experiment grid iteration
results_list <- list()
importance_list <- list()
cm_list <- list()
rf_map_list <- list()

for (i in 1:nrow(dates_grid)) {
  starttime <- Sys.time()
  target_file <- rfstack_filepaths[i] 
  file_ts <- mod_grid$Date[i]
  message(paste0(">> STARTING: RF PROCESSING, PSS ", file_ts,
                 " | ", 
                 i, " of ", nrow(dates_grid),
                 " | ",
                 format(Sys.time(), "%T %Z")
                 ))
  stackdate <- as_date(substr(file_ts, 1,8)) #reduce to only date (YMD)
  # --- SEASONAL SPLIT LOGIC ---
  if (month(stackdate) == "8") {
    print(paste("Skipping August dates:", stackdate))
    next # Instantly jumps back to the top of the loop for the next file
    print(paste("Processing:", stackdate))
  }
  if (stackdate < season_split) {
    message("    Running Spring/Early Summer Model (SAV Presence)")
    # Define season tag to use later for saving files
    season_tag <- "SP"
    true_label <- train_plots$OccSP
    valid_label <- val_plots$OccSP
    # sav_train_density <- train_data_clean %>% 
    #   filter(OccSP == "SAV")
    # print(length(true_label))
    # print(length(valid_label))
  } else {
    message("    Running Late Summer/Fall Model (SAV Presence)")
    # Define season tag to use later for saving files
    season_tag <- "FA"
    true_label <- train_plots$OccFA
    valid_label <- val_plots$OccFA
    # sav_train_density <- train_data_clean %>% 
    #   filter(OccFA == "SAV")
    # print(paste("n.Training Plots:", length(true_label)))
    # print(summary(true_label))
    # print(paste("n.Validation Plots:", length(valid_label)))
  }
  # print(paste0("Processing dataset ", i, " out of ", nrow(dates_grid)))
}
  # Extract the target and scenario for this specific iteration
  current_target <- paste0(mod_grid$Target[i], season_tag)
  # Initialize empty lists to store results of each experiment grid iteration
  results_list <- list()
  importance_list <- list()
  cm_list <- list()
  rf_map_list <- list()
  print(paste("Processing:", stackdate))
  ##----------------------------------------------------------------------------##
  ## DEFINE PREDICTOR VARIABLES ----
  ##----------------------------------------------------------------------------##
  full_stack2m <- paste0(procstk, "/",
                         file_ts, "_pss8b_rf-full-stack_2m.tif")
  if (file.exists(full_stack2m)) {
    # LOAD PREVIOUSLY PROCESSED RASTER STACK
    print(paste("Processing:", stackdate))
    predictors_2m <- rast(full_stack2m)
  } else {
    # PROCESS RAW IMAGE TILES: load, mosaic, reproject 8-band imagery
    warning(paste0(">> File not found; please create raster stack for ",
                   stackdate))
    # Read in the raster stacks for all bands plus depth stack
    full_stack <- paste0(procstk, "/tempOLD",
                         "/", file_ts, "_pss8b_rf-full-stack.tif")
    predictors_ALL <- rast(full_stack)
    predictors_2m <- mask(predictors_ALL, sav_poly)
    # > Create new raster stack with all layers
    writeRaster(predictors_2m,
                file.path(procstk,
                          paste0(file_ts, "_pss8b_rf-full-stack_2m.tif")),
                gdal = c("COMPRESS=LZW", #zips the data to save space
                         "TFW=YES",      #prevents crash at 4GB
                         "BIGTIFF=YES"), #generates sidecar text file
                overwrite = TRUE)
  }
  predictors <- predictors_2m[[c(3:5,11:13,17:25,27,30:31)]] #n=18 layers
  names(predictors)
  plot(predictors[[c(12,14,15,8)]])
  predvars <- names(predictors)
  ##----------------------------------------------------------------------------##
  ## PREPARE TRAINING DATA ----
  ## > Extract data from each raster layer by plot location
  ##----------------------------------------------------------------------------##
  #  > First for the training plots (70%)
  train_extract <- terra::extract(x = predictors, 
                                  y = train_plots, #location data
                                  ID = TRUE, df = TRUE) #remove pixels not in 
  train_data <- merge(train_plots, train_extract, by = "ID") %>%
    mutate(across(c(15:20), as.factor)) %>%
    select(-c(3:5,8:9,13:14))  #remove unused columns 
  train_data$geometry <- NULL
  names(train_data)
  #  > Then for the validation plots (30%)
  test_extract <- terra::extract(x = predictors, 
                                 y = test_plots,
                                 ID = TRUE, df = TRUE)
  test_dataprep <- test_plots %>%
    mutate(ID = 1:length(plotid))
  test_data <- merge(test_dataprep, test_extract, by = "ID") %>%
    mutate(across(c(15:20), as.factor)) %>%
    select(-c(3:5,8:9,13:14))  #remove unused columns 
  test_data$geometry <- NULL
  names(test_data)
  #  > Set the required class threshold based on target name
  required_classes <- 2
  conf_pos <- "SAV"
  target_type <- "Occurrence"
  # Count the unique classes actually present in the data for this date
  unique_classes <- unique(na.omit(train_data[[current_target]]))
  if (length(unique_classes) < required_classes) {
    message(paste("Skipping", current_target, "on", stackdate, 
                  "- Expected", required_classes, 
                  "classes but found only", length(unique_classes)))
    results_list[[i]] <- data.frame(
      Target_Classification = current_target,
      Date = stackdate,
      Season = season_tag,
      OOB_Error_Rate = NA,
      Overall_Accuracy = NA,
      Kappa = NA,
      Note = "Skipped: Only 1 class present" #add note column to explain the NAs
    )
    next #use to immediately jump to the next iteration of the loop 
  }
  ##----------------------------------------------------------------------------##
  ## TRAIN THE RANDOM FOREST ----
  ##----------------------------------------------------------------------------##
  message("   >> Training Random Forest Model...\n      ",
          format(Sys.time(), "%T %Z"))
  set.seed(29)
  current_features <- predvars
  train_data_clean <- train_data %>%
    drop_na(all_of(c(current_target, current_features))) #sav_train
  test_data_clean <- test_data %>%
    drop_na(all_of(c(current_target, current_features))) #sav_test
  # Construct model formula 
  occ_formula_str <- paste(current_target, "~",
                           paste(current_features, collapse = " + "))
  occ_model_formula <- as.formula(occ_formula_str)
  # Train the Probability Forest
  occ_model <- ranger(
    formula = model_formula, 
    data = train_data_clean, 
    probability = TRUE, 
    importance = "permutation"
  )
  print(occ_model)
  # To map the likelihood of seagrass, you specifically have to extract 
  # the "Present" column (usually column 2) to feed into terra.
  predicted_probs <- occ_model$predictions[, 2]
  ###
  # Save the Model Object
  saveRDS(occ_model, paste0(procres, "/", file_ts, "_RFmodel-pOCC_",
                            current_target, ".rds"))
  message("   >> Building Probability Map for ", stackdate, "...\n      ",
          format(Sys.time(), "%T %Z"))
  # Predict the probability map
  prob_map <- terra::predict(
    object = predictors, 
    model = occ_model, 
    fun = function(model, data, ...) {
      # Extract specifically the probability of "1" (Presence)
      # predict(model, data, ...)$predictions[, "SAV"] 
      p <- predict(model, data, ...)
      return(p$predictions)
    },
    na.rm = TRUE
  )
  ###
  # ------------------------------------------------------------------------------
  # STEP 2: THRESHOLD & CLIP THE RASTER STACK
  # ------------------------------------------------------------------------------
  # Convert the continuous probabilities into a hard binary mask (1 or NA)
  # Anything below your optimized threshold becomes NA (empty space)
  optimal_threshold <- 0.8
  presence_mask <- ifel(prob_map >= optimal_threshold, 1, NA)
  plot(prob_map)
  plot(presence_mask, add=TRUE)
  # Clip the original environmental predictor stack to remove low prob areas
  masked_predictor_stack <- mask(predictors, presence_mask)
  ###
  # ------------------------------------------------------------------------------
  # STEP 3: TRAIN & PREDICT DENSITY (CATEGORICAL)
  # ------------------------------------------------------------------------------
  # CRITICAL: Filter your tabular training data so the Density model 
  # ONLY learns from points where SAV actually exists
  sav_train_density <- train_data_clean %>% 
    filter(.data[[current_target]] == "SAV")
  # Construct model formula 
  den_formula_str <- paste0("Class", season_tag, "~",
                           paste(current_features, collapse = " + "))
  den_model_formula <- as.formula(den_formula_str)
  # Train the hard categorical forest
  density_model <- ranger(
    formula = den_model_formula, 
    data = sav_train_density, 
    probability = FALSE, 
    importance = "permutation"
  )
  # Predict the final density map using the clipped raster stack
  final_density_map <- terra::predict(
    object = masked_predictor_stack, 
    model = density_model, 
    fun = function(model, data, ...) {
      predict(model, data, ...)$predictions
    },
    na.rm = TRUE
  )
  # Re-attach your categorical labels
  levels(final_density_map) <- data.frame(
    id = c(1, 2, 3), 
    Category = c("Sparse", "Moderate", "Dense")
  )
  ###
  ###
  ###
  # # Extract the Out-Of-Bag (OOB) classification error
  # oob_error <- rf_model$prediction.error
  # # Generate the caret confusion matrix object
  # predicted_classes <- rf_model$predictions
  # actual_classes <- train_data_clean[[current_target]]
  # predicted_classes <- factor(predicted_classes,
  #                             levels = levels(actual_classes))
  # conf_matrix <- confusionMatrix(
  #   data = predicted_classes,
  #   positive = conf_pos,
  #   reference = actual_classes
  # )
  # # Extract core metrics
  # overall_acc <- conf_matrix$overall["Accuracy"]
  # print(overall_acc)
  # kappa_stat <- conf_matrix$overall["Kappa"]
  # # Save the raw confusion matrix table to our new list so you can plot it later
  # cm_list[[i]] <- list(
  #   Datetime = file_ts,
  #   Target = current_target,
  #   # Scenario = current_scenario_name,
  #   Matrix_Table = conf_matrix$table
  # )
  # # Store accuracy metrics and metadata in a temporary dataframe
  # results_list[[i]] <- data.frame(
  #   Datetime = file_ts,
  #   Target_Class = current_target,
  #   OOB_Error_Rate = round(oob_error, 4),
  #   Kappa = round(as.numeric(kappa_stat), 3),
  #   Accuracy = round((1 - oob_error)*100, 2)
  # )
  # # Store Variable Importance as dataframe
  # varimp_df <- enframe(rf_model$variable.importance,
  #                      name = "Feature", value = "Importance") %>%
  #   mutate(
  #     Datetime = file_ts,
  #     Target_Class = current_target,
  #   )
  # importance_list[[i]] <- varimp_df
  # # Generate Confusion Matrix & Build a Heat Map
  # conf_matrix <- rf_model$confusion
  # print(conf_matrix)
  # cm_df <- as.data.frame(conf_matrix)
  # cm_plot <- ggplot(cm_df, aes(x = predicted, y = true)) +
  #   geom_tile(aes(fill = Freq), color = "white", linewidth = 1) +
  #   geom_text(aes(label = Freq), vjust = 0.5,
  #             size = 5, fontface = "bold",
  #             color = ifelse(cm_df$Freq > mean(cm_df$Freq), "white", "black")) +
  #   scale_fill_gradient(low = "navajowhite",
  #                       high = "mediumpurple1",
  #                       name = "Pixel Count") +
  #   scale_x_discrete(position = "top") +
  #   scale_y_discrete(limits = rev) +
  #   labs(
  #     title = paste0("Classification Confusion Matrix: ",
  #                    stackdate),
  #     subtitle = paste0("Classification Type = ",
  #                       current_target, "  |  Accuracy: ",
  #                       round(overall_acc*100, 2), "%"),
  #     x = "RF Prediction",
  #     y = "Field Reference Data"
  #   ) +
  #   theme_minimal(base_size = 14); print(cm_plot)
}



# ---- PROCESSING STEPS ----
# 1. ----

# 2. ---- 

# 3.----
##----------------------------------------------------------------------------##
## LOAD PREDICTOR VARIABLES ----
##----------------------------------------------------------------------------##
full_stack <- paste0(procstk, "/", file_ts, "_pss8b_rf-full-stack_2m.tif")
if (file.exists(full_stack)) {
  # LOAD PREVIOUSLY PROCESSED RASTER STACK
  print(paste("Processing:", stackdate))
  predictors2m <- rast(full_stack)
  predictors <- predictors2m[[c(3:5,11:13,17:25,27,30:31)]]
} else {
  # PROCESS RAW IMAGE TILES: load, mosaic, reproject 8-band imagery
  warning(paste0(">> File not found; please create raster stack for ",
               stackdate))
  # Read in the raster stacks for all bands plus depth stack
  # rfs_pred  <- rast(target_file) #read in optical-based bands
  # rfs_depth <- rast(file.path(procwcc, paste0("TCB_", file_ts, ".tif")))
  # # rfs_depth <- rfs_depth[[3:4]]
  # # Adjust alignment on depth stack
  # depth_aligned <- resample(rfs_depth, rfs_pred, method = "bilinear")
  # print(paste("Extents match?", ext(depth_aligned) == ext(rfs_pred[[1]])))
  # # Combine into comprehensive raster stack (clipped to 0-3m depths)
  # predictors <- c(rfs_pred, depth_aligned)
  # names(predictors)[9:16] <- paste0(names(predictors)[1:8], "_dg")
  # predictors <- mask(predictors, sav_poly)
  # # plot(predictors[[c(4,12,22,29)]])
  # # > Create new raster stack with all layers
  # writeRaster(predictors,
  #             file.path(procstk,
  #                       paste0(file_ts, "_pss8b_rf-full-stack_2m.tif")),
  #             gdal = c("COMPRESS=LZW", #zips the data to save space
  #                      "TFW=YES",      #prevents crash at 4GB
  #                      "BIGTIFF=YES"), #generates sidecar text file
  #             overwrite = TRUE)
}
##----------------------------------------------------------------------------##
## RUN TURBIDITY MASK TO MASK PIXELS >90% TURBID ----
##----------------------------------------------------------------------------##
# full_stackTM <- paste0(procstk, "/", file_ts, "_pss8b_rf-full-stack_2mT90.tif")
full_stack <- paste0(paste0(procmap, "/07_mod-stack"),
                     "/", file_ts, "_pss8b_rf-full-stack.tif")
# if (file.exists(full_stackTM)) {
#   # LOAD PREVIOUSLY PROCESSED RASTER STACK
#   print(paste0(">> Existing file found; loading turbidity-masked rasters for ",
#                stackdate))
# predictors_tmALL <- rast(full_stackTM)
# } else {
#   # PROCESS RAW IMAGE TILES: load, mosaic, reproject 8-band imagery
#   print(paste0(">> File not found; creating turbidity-masked rasters for 
#                 ",
#                stackdate, "  |  ", format(Sys.time(), "%T %Z")))
#   # Read in the raster stacks for all bands plus depth stack
#   rfs_pred  <- rast(full_stack)
#   turb_threshold <- 0.9
#   turbidity_index_std <- rfs_pred[[18]]
#   turbidity_mask <- turbidity_index_std < turb_threshold
#   predictors_tmALL <- mask(rfs_pred, turbidity_mask, 
#                            maskvalues=FALSE,
#                            updatevalue=NA)
#   # > Create new raster stack with all layers
#   writeRaster(predictors_tmALL,
#               file.path(procstk,
#                         paste0(file_ts, "_pss8b_rf-full-stack_2mT90.tif")),
#               gdal = c("COMPRESS=LZW", #zips the data to save space
#                        "TFW=YES",      #prevents crash at 4GB
#                        "BIGTIFF=YES"), #generates sidecar text file
#               overwrite = TRUE)
# }
predictors_ALL <- rast(full_stack)
predictors <- mask(predictors_ALL, sav_poly)
predictors <- predictors[[c(3:5,11:13,17:25,27,30:31)]] #n=18 layers

# predictors <- predictors_tmALL[[c(3:5,11:13,17:25,27,30:31)]] #n=18 layers
names(predictors)
plot(predictors[[c(12,14,15,8)]])
# mtext(paste0("PSS Bands & Indices (0-2m), Turbidity-Masked  |  Barnegat Bay: ",
#              stackdate), 2)
predvars <- names(predictors)
# 4. ----

# 5. Convert to spatial dataframe ----

# 6. Filter plots to 0-2m depth ----

# 7. Split the dataset into Training & Validation subsets----

# 8. Visualize results ----

# ---- OUTPUTS FROM SCRIPT ----
# 

