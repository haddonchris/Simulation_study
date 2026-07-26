rm(list = ls())
graphics.off()


# LOAD LIBRARIES: ####
# =============== #

{
  library(tidyverse)
  library(caret)
  library(pROC)
  library(ranger)
  library(parallel)
  library(doSNOW)
  library(terra)
  library(sf)
  library(sdmTMB)
  library(viridis)
}

# READ IN THE DATA: ####
# ================= #

sdf <- readRDS("./Data/Simple_50_data_frames_with_models.rds")
# meep <- readRDS("./Data/simple_data_frames.rds")

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# CREATE RASTER STACKS OF THE ENVIRONMENTAL VARIABLES & ADD STATS: ####
# ================================================================ #

{
  # NDVI:
  # ----
  cat("Getting NDVI data...", "\n")

  # Read in the files as a data frame:
  ndvi_files <- data.frame(
    files = list.files(
      # List all of the file names and add as a column in a data frame
      path = "./Satellite_images/NDVI",
      pattern = "\\.tif$",
      full.names = TRUE
    )
  )

  # Create a Mean NDVI value raster
  ndvi_mean <- mean(terra::rast(ndvi_files$files)) # Read in the files as a raster stack
  names(ndvi_mean) <- "ndvi"
  rm(ndvi_files)

  # EVI
  # ---
  cat("Getting EVI data...", "\n")
  evi_files <- data.frame(
    files = list.files(
      # List all of the file names and add as a column in a data frame
      path = "./Satellite_images/EVI",
      pattern = "\\.tif$",
      full.names = TRUE
    )
  )
  evi_mean <- mean(terra::rast(evi_files$files))
  names(evi_mean) <- "evi"
  rm(evi_files)

  # Leaf Area Index:
  # ----------------
  cat("Getting LAI data...", "\n")
  lai_files <- data.frame(
    files = list.files(
      # List all of the file names and add as a column in a data frame
      path = "./Satellite_images/LAI",
      pattern = "\\.tif$",
      full.names = TRUE
    )
  )
  lai_mean <- mean(terra::rast(lai_files$files))
  names(lai_mean) <- "lai"
  rm(lai_files)

  # Percent Canopy Cover:
  # ---------------------
  cat("Getting PTC data...", "\n")
  ptc_files <- data.frame(
    files = list.files(
      # List all of the file names and add as a column in a data frame
      path = "./Satellite_images/percent_canopy_cover",
      pattern = "\\.tif$",
      full.names = TRUE
    )
  )

  # Create a Raster Stack:
  ptc_mean <- mean(terra::rast(ptc_files$files)) # Read in the files as a raster stack
  names(ptc_mean) <- "ptc"
  rm(ptc_files)

  # Percent non-tree vegetation:
  # ----------------------------
  cat("Getting PNT data...", "\n")
  pnt_files <- data.frame(
    files = list.files(
      # List all of the file names and add as a column in a data frame
      path = "./Satellite_images/percent_non_tree",
      pattern = "\\.tif$",
      full.names = TRUE
    )
  )

  # Create a Raster Stack:
  pnt_mean <- mean(terra::rast(pnt_files$files)) # Read in the files as a raster stack
  names(pnt_mean) <- "pnt"
  rm(pnt_files)

  # Percent non-vegetation
  # ----------------------
  cat("Getting PNV data...", "\n")
  pnv_files <- data.frame(
    files = list.files(
      # List all of the file names and add as a column in a data frame
      path = "./Satellite_images/percent_non_veg",
      pattern = "\\.tif$",
      full.names = TRUE
    )
  )

  pnv_mean <- mean(terra::rast(pnv_files$files))
  names(pnv_mean) <- "pnv"
  rm(pnv_files)

  # Terrain variables
  # -----------------
  cat("Getting terrain data...", "\n")

  elev <- terra::rast("./Satellite_images/Elevation/elevation_export.tif")
  names(elev) <- "elev"
  slope <- terra::rast("./Satellite_images/Slope/slope_export.tif")
  names(slope) <- "slope"

  # Soil Variables:
  # ---------------
  cat("Getting soil data...", "\n")

  soil_files <- list.files(
    path = "./Satellite_images/Soil_variables",
    pattern = "\\.tif$",
    full.names = TRUE
  )

  # Make a raster Stack
  soil_rst <- terra::rast(soil_files)
  names(soil_rst) <- c("clay", "sand", "silt")
  rm(soil_files)

  # Distance to water:
  # ------------------
  cat("Getting water distance data...", "\n")
  rivers <- st_set_crs(
    st_read(
      "./Gorongosa_files/Rivers_pans/rivers_permanent_water.shp",
      quiet = T
    ),
    4326
  )
  rivers <- st_transform(rivers, 32736)

  lake <- st_transform(
    st_read("./Gorongosa_files/Rivers_pans/lake_urema.shp", quiet = T),
    32736
  )

  # Read in a random raster:
  rst <- terra::rast("./Satellite_images/Elevation/elevation_export.tif")

  # Rasterize the images:
  lake_rst <- terra::rasterize(lake, rst, background = NA, touches = T)

  rivers_rst <- terra::rasterize(rivers, rst, background = NA, touches = T)

  # Combine the rasters:
  water <- cover(lake_rst, rivers_rst)

  # Extract the distance to water:
  water_dist <- terra::distance(water)
  names(water_dist) <- "water_dist"
  # plot(water_dist, col = inferno(150))
  rm(rivers, lake, rst, lake_rst, rivers_rst, water)

  # Resample rasters to the same resolution:
  cat("Resampling rasters...", "\n")
  rst <- terra::rast("./Satellite_images/Elevation/elevation_export.tif")

  ndvi <- resample(x = ndvi_mean, y = rst)
  evi <- resample(x = evi_mean, y = rst)
  lai <- resample(x = lai_mean, y = rst)
  ptc <- resample(x = ptc_mean, y = rst)
  pnt <- resample(x = pnt_mean, y = rst)
  pnv <- resample(x = pnv_mean, y = rst)
  elev <- resample(x = elev, y = rst)
  slope <- resample(x = slope, y = rst)
  soil <- resample(x = soil_rst, y = rst)
  water_dist <- resample(x = water_dist, y = rst)

  rm(
    rst,
    ndvi_mean,
    evi_mean,
    lai_mean,
    ptc_mean,
    pnt_mean,
    pnv_mean,
    soil_rst
  )

  # Load park boundary:
  gnp <- st_read(
    "./Gorongosa_files/GNP_shp/WDPA_WDOECM_Sep2025_Public_801_shp-polygons.shp",
    quiet = T
  )

  # Edit/Clean the data
  gnp <- gnp[1, 31]
  gnp <- st_transform(gnp, "epsg:32736") # Set the CRS

  # Crop and mask the stack:
  cat("Masking and cropping rasters...", "\n")
  ndvi <- mask(crop(ndvi, gnp), gnp)
  evi <- mask(crop(evi, gnp), gnp)
  lai <- mask(crop(lai, gnp), gnp)
  ptc <- mask(crop(ptc, gnp), gnp)
  pnt <- mask(crop(pnt, gnp), gnp)
  pnv <- mask(crop(pnv, gnp), gnp)
  elev <- mask(crop(elev, gnp), gnp)
  slope <- mask(crop(slope, gnp), gnp)
  soil <- mask(crop(soil, gnp), gnp)
  water_dist <- mask(crop(water_dist, gnp), gnp)

  stack <- c(
    ndvi,
    evi,
    lai,
    ptc,
    pnt,
    pnv,
    elev,
    slope,
    soil,
    water_dist
  )

  prob_rst <- c(ndvi, lai, water_dist, ptc)

  cat("Calculating raster statistics...", "\n")
  ndvi_stats <- global(ndvi, fun = c("mean", "sd"), na.rm = TRUE)
  water_stats <- global(water_dist, fun = c("mean", "sd"), na.rm = TRUE)
  ptc_stats <- global(ptc, fun = c("mean", "sd"), na.rm = TRUE)
  lai_stats <- global(lai, fun = c("mean", "sd"), na.rm = TRUE)

  # Pre-extract stats columns and convert to vectors for faster access:
  ndvi_mean <- ndvi_stats[1, 1]
  ndvi_sd <- ndvi_stats[1, 2]

  water_mean <- water_stats[1, 1]
  water_sd <- water_stats[1, 2]

  ptc_mean <- ptc_stats[1, 1]
  ptc_sd <- ptc_stats[1, 2]

  lai_mean <- lai_stats[1, 1]
  lai_sd <- lai_stats[1, 2]

  # Remove unnecessary objects
  rm(
    elev,
    evi,
    pnt,
    pnv,
    slope,
    soil
  )
  cat("Done!!!")
}


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# CREATE SAMPLE POINTS TO EXTRACT PROBABILITIES FROM THE RASTERS: ####
# =============================================================== #

# Create a hexagonal grid of sampling points to extract the probabilities from the rasters:
hex <- st_make_grid(
  gnp,
  cellsize = 2500, # 2.5 km
  square = FALSE
)

pts <- st_centroid(hex) # Get the centroids of the hexagonal grid cells to use as sampling points
pts <- as.data.frame(st_coordinates(pts[gnp])) # Convert the matrix to a data frame
colnames(pts) <- c("x", "y") # Rename the columns to "

rm(hex)

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# EXTRACT PROBABILITIES FROM THE RASTERS FOR EACH MODEL: ####
# ====================================================== #

# Create a data frame to store the probabilities for each model:
probs <- data.frame(
  mod_id = factor(),
  true_p = numeric(),
  lr_p = numeric(),
  rf_p = numeric()
)

# Define a function to extract the probabilities from the random forest model:
rf_fun <- function(model, data, ...) {
  p <- predict(model, data = data, type = "response")$predictions

  if (is.matrix(p) || is.data.frame(p)) {
    p <- p[, "1"] # replace with the presence-class name if needed
  }
  as.numeric(p)
}

system.time({
  for (i in 1:length(sdf)) {
    # Print the progress of the model fitting and evaluation
    cat("Calculating probabilities for model", i, "of", length(sdf), "\n")

    # -------------------------------------------------------------------------- #

    # Make a vector containing the model ID for the current iteration:
    mod_id <- sdf[[i]]$mod_id

    # -------------------------------------------------------------------------- #

    # CALCULATE THE TRUE PROBABILITIES FOR THE CURRENT MODEL:
    # ------------------------------------------------------- #

    # Scale the coefficients by the standard deviation of each raster:
    b_ndvi <- sdf[[i]]$equation$beta_ndvi / ndvi_sd
    b_water <- sdf[[i]]$equation$beta_water / water_sd
    b_ptc <- sdf[[i]]$equation$beta_ptc / ptc_sd
    b_lai <- sdf[[i]]$equation$beta_lai / lai_sd

    logit_rst <- (sdf[[i]]$equation$alpha +
      (b_ndvi * stack$ndvi) +
      (b_water * stack$water_dist) +
      (b_ptc * stack$ptc)) +
      (b_lai * stack$lai)

    p_true_rst <- 1 / (1 + exp(-logit_rst))

    names(p_true_rst) <- "true_p"
    varnames(p_true_rst) <- "true_p"

    vect_p_true <- terra::vect(pts, geom = c("x", "y"), crs = "epsg:32736")

    var <- terra::extract(p_true_rst, vect_p_true, ID = FALSE)

    ture_p <- var$true_p

    # -------------------------------------------------------------------------- #

    # CALCULATE THE LOGISTIC REGRESSION PROBABILITIES FOR THE CURRENT MODEL:
    # ---------------------------------------------------------------------- #

    lr_rst <- terra::predict(
      stack,
      sdf[[i]]$lr_mod,
      type = "response",
      cores = 7
    )
    names(lr_rst) <- "lr_p"
    varnames(lr_rst) <- "lr_p"

    # plot(lr_rst, col = inferno(150))

    var <- terra::extract(lr_rst, vect_p_true, ID = FALSE)

    lr_p <- var$lr_p

    # -------------------------------------------------------------------------- #

    # CALCULATE THE LOGISTIC REGRESSION PROBABILITIES FOR THE CURRENT MODEL:
    # ---------------------------------------------------------------------- #

    rf_rst <- terra::predict(
      stack,
      sdf[[i]]$rf_mod,
      type = "response",
      fun = rf_fun,
      cores = 7,
      cpkgs = "ranger",
      na.rm = TRUE
    )

    na_mask <- is.na(stack[[1]])
    rf_rst[na_mask] <- NA
    # rf_rst <- rf_rst[[2]]

    names(rf_rst) <- "rf_p"
    varnames(rf_rst) <- "rf_p"

    var <- terra::extract(rf_rst, vect_p_true, ID = FALSE)

    rf_p <- var$rf_p

    # -------------------------------------------------------------------------- #

    # ADD THE PROBABILITIES TO THE DATA FRAME:
    # --------------------------------------- #

    temp <- data.frame(
      mod_id = mod_id,
      true_p = ture_p,
      lr_p = lr_p,
      rf_p = rf_p
    )

    probs <- rbind(probs, temp)
  }
})

write.csv(probs, "./Data/simple_50_probabilities.csv", row.names = FALSE)

rm(
  gnp,
  lai_stats,
  ndvi_stats,
  ptc_stats,
  water_stats,
  pts,
  temp,
  var,
  b_lai,
  b_ndvi,
  b_ptc,
  b_water,
  logit_rst,
  p_true_rst,
  lr_rst,
  rf_rst,
  vect_p_true,
  i,
  lai,
  ndvi,
  ptc,
  water_dist,
  lai_mean,
  lai_sd,
  ndvi_mean,
  ndvi_sd,
  ptc_mean,
  ptc_sd,
  water_mean,
  water_sd,
  ture_p,
  lr_p,
  rf_p,
  na_mask,
  prob_rst,
  mod_id
)

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# CALCULATE THE ERRORS FOR EACH MODEL: ####
# ==================================== #

probs$lr_error <- (probs$true_p - probs$lr_p)
probs$rf_error <- (probs$true_p - probs$rf_p)


test <- probs |>
  dplyr::mutate(
    lr_error = true_p - lr_p,
    rf_error = true_p - rf_p,
    mod_id = factor(mod_id)
  ) |>
  dplyr::group_by(mod_id) |>
  dplyr::summarise(
    lr_rmse = sqrt(mean(abs(lr_error)^2)),
    lr_mae = mean(abs(lr_error)),
    lr_bias = mean(lr_error),
    lr_cor = cor(lr_p, true_p),
    rf_rmse = sqrt(mean(abs(rf_error)^2)),
    rf_mae = mean(abs(rf_error)),
    rf_bias = mean(rf_error),
    rf_cor = cor(rf_p, true_p),
    rmse_error = lr_rmse - rf_rmse,
    mae_error = lr_mae - rf_mae
  )


# RMSE Stats:
mean(test$lr_rmse)
mean(test$rf_rmse)

sd(test$lr_rmse)
sd(test$rf_rmse)

# MAE Stats:
mean(test$lr_mae)
mean(test$rf_mae)
sd(test$lr_mae)
sd(test$rf_mae)

# Bias Stats:
mean(test$lr_bias)
mean(test$rf_bias)
sd(test$lr_bias)
sd(test$rf_bias)

# Correlation Stats:
mean(test$lr_cor)
mean(test$rf_cor)
sd(test$lr_cor)
sd(test$rf_cor)

range(probs$true_p)
range(probs$lr_error)
range(probs$rf_error)

test_longer <- pivot_longer(
  test,
  cols = c(lr_rmse, rf_rmse),
  names_to = "rmse",
  values_to = "value"
)

ggplot(test_longer, aes(x = value, fill = rmse)) +
  geom_density(alpha = 0.5)


ggplot(test_longer, aes(x = value, fill = rmse)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of RMSE values for Logistic Regression and Random Forest models",
    x = "RMSE",
    y = "Density",
    fill = "Model"
  ) +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Logistic Regression", "Random Forest")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )


ggplot(test_longer, aes(x = rmse, y = value, fill = rmse)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of RMSE values for Logistic Regression and Random Forest models",
    x = "Model",
    y = "RMSE",
    fill = "Model"
  ) +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Logistic Regression", "Random Forest")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggplot(test_longer, aes(x = value, fill = rmse)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 15)

test_longer_bias <- pivot_longer(
  test,
  cols = c(lr_bias, rf_bias),
  names_to = "bias",
  values_to = "value"
)

ggplot(test_longer_bias, aes(x = value, fill = bias)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of Bias values for Logistic Regression and Random Forest models",
    x = "Bias",
    y = "Density",
    fill = "Model"
  ) +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Logistic Regression", "Random Forest")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

test_longer_mae <- pivot_longer(
  test,
  cols = c(lr_mae, rf_mae),
  names_to = "mae",
  values_to = "value"
)

ggplot(test_longer_mae, aes(x = value, fill = mae)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of MAE values for Logistic Regression and Random Forest models",
    x = "MAE",
    y = "Density",
    fill = "Model"
  ) +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Logistic Regression", "Random Forest")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

test_longer_cor <- pivot_longer(
  test,
  cols = c(lr_cor, rf_cor),
  names_to = "cor",
  values_to = "value"
)

ggplot(test_longer_cor, aes(x = value, fill = cor)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of Correlation values for Logistic Regression and Random Forest models",
    x = "Correlation",
    y = "Density",
    fill = "Model"
  ) +
  scale_fill_manual(
    values = c("blue", "red"),
    labels = c("Logistic Regression", "Random Forest")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )
