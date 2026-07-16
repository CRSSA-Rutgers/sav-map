# ============================================================================= #
# Script Name:  01_clean_ref-data.R
# Description:  Cleans field reference datasets to use in training & validating 
#               Random Forest models and predictions for mapping seagrass (SAV).
# Author:       Jess M. Stitt
# Date:         2026-03-31
# ============================================================================= #

# ============================================================================= #
# ---- TASK DESCRIPTION ----
# Input(s):
#  > [CSV] Field Reference Dataset ["raw"]
#     (SAV plot data from spring & fall surveys in BB-LEH, 2023 | DEP & Stockton)
#  > [SHP] SAV Suitable Habitat Polygon
#     (All areas 0-2m deep across BB-LEH, suitable for seagrass | USGS + CRSSA)
# Output(s): 
#  > [CSV, SHP] Field Reference Dataset ["clean"], subset by Training & Validation
#     (Datasets organized with variables usable for RF modeling; split the 
#     plots into Training and Validation subsets to teach then test models)
# ============================================================================= #
# ---- LIBRARIES SPECIFIC TO THIS SCRIPT ----
library(readxl)     #for reading in excel files
library(ggalluvial) #for creating alluvial plots (compare spring vs. fall data)
# ============================================================================= #
# ---- USER-DEFINED VARIABLES FOR SCRIPT ----
#  <NONE>
# ============================================================================= #
# ---- PROCESSING WORKFLOW ----
##------------------------------------------------------------------------------#
## 1. LOAD & ORGANIZE ----
## >> Load field data as dataframes & organize columns for use in modeling
##------------------------------------------------------------------------------#
tidy_field_data <- field_data_subplots %>%
  # > Match values to names to pivot into rows (long-form table)
  pivot_longer(
    cols = matches("_([bmrz])([1-8])$"), 
    names_to = c("season", "cover_init", "subplot"),
    names_pattern = "(.*)_([bmrz])([1-8])",
    values_to = "pct_cov"
  ) %>%
  # > Map initials to ecological names for clean reporting
  mutate(
    cover_type = case_when(
      cover_init == "b" ~ "BARE",
      cover_init == "m" ~ "MACR",
      cover_init == "r" ~ "RUPP",
      cover_init == "z" ~ "ZOST"
    )
  ) %>%
  # > Pivot cover types back into columns for math
  pivot_wider(
    names_from = cover_type,
    values_from = pct_cov
  )
##------------------------------------------------------------------------------#
## 2. SUMMARIZE PLOT DATA ---- 
## >> Average subplots (n=8) per XY & Date for one non-zero %COVER value by type
##------------------------------------------------------------------------------#
# > Calculate summary data of all subplots, by plotID and season
sav_plot_summary <- tidy_field_data %>%
  group_by(plotid, season) %>%
  summarise(
    # Non-zero mean cover percentage (ignoring 0s and NAs)
    across(c(BARE, MACR, RUPP, ZOST), 
           ~ {
             non_zeros <- .x[.x > 0 & !is.na(.x)]
             if(length(non_zeros) == 0) 0 else mean(non_zeros)
           }, 
           .names = "avgNZ_{.col}"),
    # Heterogeneity (~patchiness) via standard deviation
    across(c(RUPP, ZOST), 
           ~ sd(.x, na.rm = TRUE), 
           .names = "sd_{.col}"),
    .groups = "drop"
  ) %>%
  # Calculate Total SAV and Dominant Class with the new columns
  mutate(
    sav_RUZO = avgNZ_ZOST + avgNZ_RUPP,
    dominant_class = case_when(
      avgNZ_ZOST > avgNZ_RUPP & 
        avgNZ_ZOST > avgNZ_MACR & 
        avgNZ_ZOST > avgNZ_BARE ~ "Zostera",
      avgNZ_RUPP > avgNZ_ZOST & 
        avgNZ_RUPP > avgNZ_MACR & 
        avgNZ_RUPP > avgNZ_BARE ~ "Ruppia",
      avgNZ_MACR > avgNZ_ZOST & 
        avgNZ_MACR > avgNZ_RUPP & 
        avgNZ_MACR > avgNZ_BARE ~ "Macroalgae",
      avgNZ_BARE > avgNZ_ZOST & 
        avgNZ_BARE > avgNZ_RUPP & 
        avgNZ_BARE > avgNZ_MACR ~ "Bare",
      sav_RUZO > avgNZ_BARE & 
        sav_RUZO >= avgNZ_MACR ~ "Mixed_SAV",
      TRUE ~ "Transition"
    )
  ) %>%
  drop_na(sd_RUPP)
# > Pivot table to wide format again (collapse rows)
wide_training_points <- sav_plot_summary %>%
  select(
    plotid,
    season,
    dominant_class,
    pRUZO = sav_RUZO,
    pMACR = avgNZ_MACR,
    pBARE = avgNZ_BARE,
    pRUPP = avgNZ_RUPP,
    pZOST = avgNZ_ZOST,
    sd_RUPP,
    sd_ZOST
  ) %>%
  pivot_wider(
    # ANCHORS (id_cols)
    id_cols = plotid, 
    names_from = season,
    values_from = c(
      dominant_class,
      # sav_RUZO,
      pRUZO,
      pBARE, pMACR, pRUPP, pZOST
    ),
    # Rename columns using values of interest
    names_glue = "{season}_{.value}"
  )
##------------------------------------------------------------------------------#
## 3. CATEGORIZE SAV ---- 
## >> Separate seagrass into density classes (sparse, mod, dense) by %COVER value
##------------------------------------------------------------------------------#
wide_training_all <- wide_training_points %>%
  mutate(
    sp23_density_class = case_when(
      is.na(sp23_pRUZO ) ~ NA_character_,
      sp23_pRUZO  < 10 ~ "Bare",
      sp23_pRUZO  >= 10  & sp23_pRUZO < 40 ~ "Sparse",
      sp23_pRUZO  >= 40 & sp23_pRUZO  < 80 ~ "Moderate",
      sp23_pRUZO  >= 80 ~ "Dense"
    ),
    fa23_density_class = case_when(
      is.na(fa23_pRUZO) ~ NA_character_,
      fa23_pRUZO  < 10 ~ "Bare",
      fa23_pRUZO  >= 10  & fa23_pRUZO < 40 ~ "Sparse",
      fa23_pRUZO  >= 40 & fa23_pRUZO  < 80 ~ "Moderate",
      fa23_pRUZO  >= 80 ~ "Dense"
    )
  ) %>%
  mutate(
    across(
      c(sp23_density_class, fa23_density_class),
      ~ factor(., levels = c("Bare", "Sparse", "Moderate", "Dense"))
    )
  )
# > Verify the new columns exist and the NAs were handled correctly
head(wide_training_all %>%
       select(plotid, sp23_density_class, fa23_density_class), n = 8)
##------------------------------------------------------------------------------#
## 4. COMPARE SEASONAL DIFFERENCES ---- 
## >> Contrast spring and fall survey data across plots
##------------------------------------------------------------------------------#
# > Isolate Paired Plots: Keep only rows that have data for BOTH seasons
paired_plots_den <- wide_training_all %>%
  drop_na(sp23_density_class, fa23_density_class)
paired_plots_dom <- wide_training_all %>%
  drop_na(sp23_dominant_class, fa23_dominant_class)
# > Build Transition Matrix
transition_dominance <- paired_plots_dom %>%
  count(sp23_dominant_class, fa23_dominant_class, name = "nplots_TYPE") %>%
  # > Calculate what percentage of re-sampled plots each transition represents
  mutate(
    pctTYPE = round((nplots_TYPE / sum(nplots_TYPE)) * 100, 1)) %>%
  arrange(desc(nplots_TYPE))
transition_density <- paired_plots_den %>%
  count(sp23_density_class, fa23_density_class, name = "nplots_DENSITY") %>%
  mutate(
    pctDENSITY = round((nplots_DENSITY / sum(nplots_DENSITY)) * 100, 1)) %>%
  arrange(desc(nplots_DENSITY))
##------------------------------------------------------------------------------#
## 5. CONVERT TO SPATIAL DATAFRAME ---- 
## >> Make dataset spatially explicit for mapping & model development
##------------------------------------------------------------------------------#
## > RENAME & REORDER: Move key fields to front and give new names
plots_sf <- field_data_subplots %>%
  select(
    plotid,
    # season,
    split_assign,
    longitude,
    latitude,
    yr = year, #newName = oldName
    sp23_sampTF = sp23_sampling,
    sp23_mm,
    sp23_dd,
    sp23_htm = sp23_depth_m,
    fa23_sampTF = fa23_sampling,
    fa23_mm,
    fa23_dd,
    fa23_htm = fa23_depth_m
  )
## > PREP GEOMETRY: organize coordinates for spatial manipulation
field_points_tbl <- plots_sf %>%
  mutate(
    X = coalesce(across(matches("^long|^lon|^x$", ignore.case = TRUE))),
    Y = coalesce(across(matches("^lat|^y$", ignore.case = TRUE)))
  ) %>%
  mutate(
    X = as.numeric(X[[1]]),
    Y = as.numeric(Y[[1]])
  ) %>%
  filter(!is.na(X) & !is.na(Y)) %>%
  as_tibble()
if(any(field_points_tbl$X > 0)) {
  warning("!! ALERT !! Positive Longitudes detected...
Barnegat Bay requires negative Longitude (e.g., -74.1); check your CSV")
}
## > CONVERT TO SF OBJECT
points_sf <- st_as_sf(field_points_tbl, coords = c("longitude", "latitude"), 
                      crs = 4326)
## > MERGE BENTHIC COVER DATA WITH SPATIAL GPS POINTS
spatial_sf <- points_sf %>%
  inner_join(wide_training_all, by = "plotid") %>%
  select(
    plotid,
    everything() #keep all other columns
  )
points_proj <- st_transform(spatial_sf, st_crs(sav_poly)) 
##------------------------------------------------------------------------------#
## 6. FILTER PLOTS TO SAV ZONE ---- 
## >> Select only plot locations that are found between 0 to 2m below sea level
##------------------------------------------------------------------------------#
sav_st <- sav_poly_sf
points_clipped <- st_filter(points_proj, sav_st) #0-2m depth for SAV Zone
points_omitted <- points_proj[lengths(st_intersects(points_proj, sav_st)) == 0, ]
n_original <- nrow(points_proj) #all points
n_kept     <- nrow(points_clipped) #clipped points (kept for modeling)
n_omitted  <- nrow(points_omitted) #removed points
message(paste("Successfully loaded all field reference plots for 2023",
              "\n >> Original Points:", n_original,
              "\n >> Points in Water:", n_kept,
              "\n >> Dropped (Deep/Bad GPS):", n_original - n_kept))
##------------------------------------------------------------------------------#
## 7. SPLIT PLOTS INTO TRAINING & VALIDATION SUBSETS ---- 
## >> Divide the dataset by 70% & 30% into Training & Validation subsets
##------------------------------------------------------------------------------#
poly_data <- as.data.frame(points_clipped)
poly_data$ID <- 1:nrow(poly_data)
classify_plotsPROC <- poly_data %>%
  group_by(plotid) %>%
  mutate(date_s23 = make_date(yr, sp23_mm, sp23_dd)) %>%
  mutate(date_f23 = make_date(yr, fa23_mm, fa23_dd)) %>%
  select(
    ID,
    plotid,
    splt_rf = split_assign,
    date_s23,
    deepm_s23 = sp23_htm,
    s23_TF = sp23_sampTF,
    sand_s23 = sp23_pBARE,
    savruzo_s23 = sp23_pRUZO,
    macroalga_s23 = sp23_pMACR,
    ruppia_s23 = sp23_pRUPP,
    zostera_s23 = sp23_pZOST, 
    dmcov_s23 = sp23_dominant_class,
    den_s23 = sp23_density_class,
    date_f23,
    deepm_f23 = fa23_htm,
    f23_TF = fa23_sampTF,
    sand_f23 = fa23_pBARE,
    savruzo_f23 = fa23_pRUZO,
    macroalga_f23 = fa23_pMACR,
    ruppia_f23 = fa23_pRUPP,
    zostera_f23 = fa23_pZOST,
    dmcov_f23 = fa23_dominant_class,
    den_f23 = fa23_density_class,
    X,
    Y,
    geometry
  ) %>%
  mutate(OccSP = case_when(
    savruzo_s23 >= 10 ~ "SAV",
    .default = "No_SAV" 
  )) %>% mutate(OccSP = as.factor(OccSP)) %>%
  mutate(OccFA = case_when(
    savruzo_f23 >= 10 ~ "SAV",
    .default = "No_SAV" 
  )) %>% mutate(OccFA = as.factor(OccFA)) %>%
  mutate(ClassSP = case_when(
    savruzo_s23 >= 80 ~ "SAV_dense",
    savruzo_s23 >= 40 & savruzo_s23 < 80 ~ "SAV_moderate",
    savruzo_s23 >= 10 & savruzo_s23 < 40 ~ "SAV_sparse",
    # sand_s23 >= 90 ~ "Bare",
    .default = "No_SAV" 
  )) %>% mutate(ClassSP = as.factor(ClassSP)) %>%
  mutate(ClassFA = case_when(
    savruzo_f23 >= 80 ~ "SAV_dense",
    savruzo_f23 >= 40 & savruzo_f23 < 80 ~ "SAV_moderate",
    savruzo_f23 >= 10 & savruzo_f23 < 40 ~ "SAV_sparse", 
    # sand_f23 > 90 ~ "Bare",
    .default ="No_SAV"
  )) %>% mutate(ClassFA = as.factor(ClassFA)) %>%
  mutate(DomCovSP = case_when(
    dmcov_s23 == "Ruppia" ~ "RUPP",
    dmcov_s23 == "Zostera" ~ "ZOST",
    dmcov_s23 == "Mixed_SAV" ~ "SAVmx",
    .default = "Not_SAV" 
  )) %>% mutate(DomCovSP = as.factor(DomCovSP)) %>%
  mutate(DomCovFA = case_when(
    dmcov_f23 == "Ruppia" ~ "RUPP",
    dmcov_f23 == "Zostera" ~ "ZOST",
    dmcov_f23 == "Mixed_SAV" ~ "SAVmx",
    .default = "Not_SAV" 
  )) %>% mutate(DomCovFA = as.factor(DomCovFA)) %>%
  mutate(across(c("OccFA",
                  "ClassFA",
                  # "RuppFA",
                  # "ZostFA",
                  "DomCovFA"),
                ~ case_when(
                  is.na(f23_TF) ~ .x,
                  f23_TF == FALSE ~ NA,
                  TRUE ~ .x
                )))
classify_plots <- classify_plotsPROC %>%
  select(-c(sand_s23,savruzo_s23,macroalga_s23,ruppia_s23,zostera_s23,
            sand_f23,savruzo_f23,macroalga_f23,ruppia_f23,zostera_f23))
# Subset the rows based exactly on your column text
train_plots <- classify_plots[classify_plots$splt_rf == "training", ] %>%
  select(-splt_rf)
val_plots <- classify_plots[classify_plots$splt_rf == "validation", ] %>%
  select(-splt_rf)
n_train <- nrow(train_plots) 
n_val   <- nrow(val_plots) 
##------------------------------------------------------------------------------#
## 8. VISUALIZE RESULTS & GENERATE DIAGNOSTICS ---- 
## >> Develop graphics for evaluation and visualization
##------------------------------------------------------------------------------#
### > PLOT 00a: map of plots against 0-2m depth cutoff ----
plots_clip2023 <- ggplot() +
  geom_spatvector(data = water_poly, 
                  aes(fill = "Deep Water (>2m)"), color = "cadetblue4") +
  geom_spatvector(data = sav_poly, 
                  aes(fill = "Shallow/SAV (0-2m)"), color = "cadetblue3") +
  scale_fill_manual(name = "Water Extent",
                    values = c(
                      "Deep Water (>2m)" = "cadetblue4", 
                      "Shallow/SAV (0-2m)" = "cadetblue3")) +
  geom_spatvector(data = water_poly, fill = "transparent", color = "black", 
                  linewidth = 0.6) +
  geom_sf(data = points_proj, color = "darkred", size = 1, alpha = 0.8) +
  geom_sf(data = points_clipped, color = "chartreuse", size = 1) +
  labs(title = "Barnegat Bay 2023 Survey Points, Clipped",
       subtitle = "Green = Kept (<2m Depth)  |  Red = Dropped (Outside SAV Zone)",
       caption = paste(n_original - n_kept, "points removed.")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    plot.title = element_text(face = "bold")); plots_clip2023
### > PLOT 00b: map of reference data plot locations ----
plots_maps2023 <- ggplot() +
  geom_spatvector(data = water_poly, 
                  aes(fill = "Deep Water (>2m)"), color = "cadetblue4") +
  geom_spatvector(data = sav_poly, 
                  aes(fill = "Shallow/SAV (0-2m)"), color = "cadetblue3") +
  scale_fill_manual(name = "Water Extent",
                    values = c(
                      "Deep Water (>2m)" = "cadetblue4", 
                      "Shallow/SAV (0-2m)" = "cadetblue3")) +
  geom_spatvector(data = water_poly, fill = "transparent", color = "black", 
                  linewidth = 0.6) +
  geom_sf(data = points_clipped, color = "darkorange2", size = 1) +
  labs(title = "Barnegat Bay Reference Datasets 2023, Alignment Check", 
       subtitle = paste0(
         "Orange: Reference Plots (n = ", 
         nrow(points_clipped), ")"),
       caption = paste(n_kept, "points retained.")) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    plot.title = element_text(face = "bold")); plots_maps2023
### > PLOT 00c: map of training vs validation plot locations ----
tvsplit_maps2023 <- 
  ggplot() +
  geom_spatvector(data = water_poly, 
                  aes(fill = "Deep Water (>2m)"), color = "cadetblue4") +
  geom_spatvector(data = sav_poly, 
                  aes(fill = "Shallow/SAV (0-2m)"), color = "cadetblue3") +
  scale_fill_manual(name = "Water Extent",
                    values = c(
                      "Deep Water (>2m)" = "cadetblue4", 
                      "Shallow/SAV (0-2m)" = "cadetblue3")) +
  geom_spatvector(data = water_poly, fill = "transparent", color = "black", 
                  linewidth = 0.6) +
  geom_sf(data = points_clipped, 
          aes(color = split_assign), size = 1) +
  scale_color_manual(name = "Plot Assignment",
                     values = c("training" = "gold1", 
                                "validation" = "maroon2")
                     ) +
  labs(title = "Barnegat Bay 2023, Assignment for Random Forest Modeling",
       subtitle = paste0(
         "All plots (n=", n_kept,
         ") assigned to model training [Yellow] or validation [Magenta]"),
       caption = paste0(" n.Training = ", n_train, 
                        "  |  n.Validation = ", n_val),
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    plot.title = element_text(face = "bold")); tvsplit_maps2023
### > PLOT 00d: seasonal differences by SAV density ----
p_transition_density <- ggplot(
  data = transition_density,
  aes(axis1 = sp23_density_class, 
      axis2 = fa23_density_class, 
      y = nplots_DENSITY)) +
  # Draw the flowing bands between the seasons
  geom_alluvium(aes(fill = sp23_density_class), alpha = 0.75) +
  # Draw the stacked bars for each season
  geom_stratum(width = 0.25, fill = "gray80", color = "black") +
  # Add the text labels to the bars
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)), fontface = "bold") +
  scale_fill_manual(
    values = c("Bare" = "#f7fcb9", "Sparse" = "#addd8e", 
               "Moderate" = "#31a354", "Dense" = "#005a32")
  ) +
  scale_x_discrete(limits = c("Spring Survey", "Fall Survey"), 
                   expand = c(0.15, 0.05)) +
  labs(
    title = "Seasonal SAV Density Transitions",
    subtitle = "Tracking 60 paired survey plots across Barnegat Bay",
    y = "Number of Plots"
  ) +
  theme_minimal() +
  theme(legend.position = "none"); p_transition_density
### > PLOT 00e: seasonal differences by dominant cover type ----
p_transition_dominance <- ggplot(
  data = transition_dominance,
  aes(axis1 = sp23_dominant_class, 
      axis2 = fa23_dominant_class, 
      y = nplots_TYPE)) +
  # Draw the flowing bands between the seasons
  geom_alluvium(aes(fill = sp23_dominant_class), alpha = 0.75) +
  # Draw the stacked bars for each season
  geom_stratum(width = 0.25, fill = "gray80", color = "black") +
  # Add the text labels to the bars
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum)), fontface = "bold") +
  scale_fill_manual(
    values = c("Bare" = "gold3", "Macroalgae" = "tomato3", 
               "Ruppia" = "orchid2", "Zostera" = "dodgerblue2")
  ) +
  scale_x_discrete(limits = c("Spring Survey", "Fall Survey"), 
                   expand = c(0.15, 0.05)) +
  labs(
    title = "Seasonal Dominant Cover Transitions",
    subtitle = "Tracking 60 paired survey plots across Barnegat Bay",
    y = "Number of Plots"
  ) +
  theme_minimal() +
  theme(legend.position = "none"); p_transition_dominance
### > DIAGNOSTICS REPORT: export maps and stats that summarize the data ----
# > 00a = Map of all plots & which were omitted
ggsave(file.path(procimg, paste0("00a_2023ref-plots-XY_all", ".png")),
       plots_clip2023, width = 9, height = 16, dpi = 300)
# > 00b = Map of plots 0-2m depth only (used for modeling)
ggsave(file.path(procimg, paste0("00b_2023ref-plots-XY_0to2mDepth", ".png")),
       plots_maps2023, width = 9, height = 16, dpi = 300)
# > 00c = Map of plots by model assignment (training vs validation)
ggsave(file.path(procimg, paste0("00c_2023ref-plots-XY_mod-split", ".png")),
       tvsplit_maps2023, width = 9, height = 16, dpi = 300)
# > 00d = Alluvial plot showing seasonal transitions in density 
ggsave(file.path(procimg, paste0("00d_2023ref-plots_seas-tdensity", ".png")),
       p_transition_density, width = 12, height = 9, dpi = 300)
# > 00e = Alluvial plot showing seasonal transitions in dominance
ggsave(file.path(procimg, paste0("00e_2023ref-plots_seas-tdominance", ".png")),
       p_transition_dominance, width = 12, height = 9, dpi = 300)
# ============================================================================= #
# ---- OUTPUTS FROM SCRIPT ----
##-----------------------------------------------------------------------------##
## Training dataset (70% of total plots)
# > CSV of training plots data (= attribute table with coords)
write.csv(train_plots, 
  file = file.path(proctrn, "sav23_training-plots.csv"),
  row.names = FALSE)
# > Shapefile of training points for the modeling script to use later
st_write(train_plots, file.path(proctrn, "sav23_training-plots.shp"),
         quiet = TRUE, append = FALSE)
## Validation dataset (30% of total plots)
# > CSV of validation plots data (= attribute table with coords)
write.csv(val_plots, 
  file = file.path(procval, "sav23_validation-plots.csv"),
  row.names = FALSE)
# > Shapefile of validation points for the modeling script to use later
st_write(val_plots, file.path(procval, "sav23_validation-plots.shp"),
         quiet = TRUE, append = FALSE)
## REF PLOT DATA -- CLIPPED 0-2M AOI ONLY
# > CSV of plots within modeling scope (= attribute table with coords)
write.csv(points_clipped, 
  file = file.path(procmod, "2023fieldPoints_2m-only.csv"),
  row.names = FALSE)
# > Shapefile of the points for the modeling script to use later
st_write(points_clipped,
         file.path(procmod, "2023fieldPoints_2m-only.shp"),
         quiet = TRUE, append=FALSE)
## ALL REF PLOT DATA -- UNCLIPPED
# > CSV of all plots (= attribute table with coords)
write.csv(points_proj, 
          file = file.path(procmsc, "2023fieldPoints_allSurveyed.csv"),
          row.names = FALSE)
# > Shapefile of all points for reference
st_write(points_proj,
         file.path(procmsc, "2023fieldPoints_allSurveyed.shp"),
         quiet = TRUE, append=FALSE)
## OMITTED REF PLOT DATA -- PLOTS REMOVED
# > CSV of the omitted plots (= attribute table with coords)
write.csv(points_omitted, 
          file = file.path(procmsc, "2023fieldPoints_omitted.csv"),
          row.names = FALSE)
# > Shapefile of omitted points for reference
st_write(points_omitted,
         file.path(procmsc, "2023fieldPoints_omitted.shp"),
         quiet = TRUE, append=FALSE)
# ============================================================================= #
