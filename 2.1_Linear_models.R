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

sdf <- readRDS("./Data/Nonlinear_200_data_frames.rds")
# sdf <- readRDS("./Data/Nonlinear_1000_data_frames.rds")

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# EXTRACT STATELLITE VARIABLES AND EXTRACT DATA FOR EACH DATA FRAME: ####
# ================================================================== #

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

  # Remove unnecessary objects
  rm(
    elev,
    evi,
    pnt,
    pnv,
    slope,
    soil,
    gnp,
    lai,
    ndvi,
    ptc,
    water_dist
  )

  cat("Done!!!")
}


for (i in 1:length(sdf)) {
  # Extract satellite imagery data (vectorized):
  coords <- st_coordinates(sdf[[i]]$df)
  vals <- terra::extract(stack, terra::vect(coords), df = TRUE)
  sdf[[i]]$df <- cbind(sdf[[i]]$df, vals[, -1])

  # Keep track of the loop's progress
  if (i %% 10 == 0) {
    cat("Extracting for data frame", i, "of", length(sdf), "\n")
  }
}

rm(coords, vals, i, stack)


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# EXTRACT MODEL PERFORMANCE METRICS FOR LOGISTIC REGRESSION: ####
# ========================================================== #

# Create an empty data frame to store the model performance metrics:
model_performance <- data.frame(
  Model_type = character(),
  model_id = character(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  TSS = numeric(),
  AUC = numeric(),
  Kappa = numeric(),
  Accuracy = numeric(),
  AccuracyLower = numeric(),
  AccuracyUpper = numeric(),
  stringsAsFactors = FALSE
)


# Number of data frames/models to evaluate
n <- length(sdf) # Number of data frames/models to evaluate


system.time({
  {
    # Set seed once for reproducible, but iteration-varying, data splits
    set.seed(3791)

    # Train the models and calculate the performance metrics:
    for (i in 1:length(sdf)) {
      cat("Running model", i, "of", length(sdf), "\n") # Print the progress of the model fitting and evaluation

      # Split the data into training and testing sets:
      train_ind <- createDataPartition(sdf[[i]]$df$presence, p = 0.7, list = F)
      train <- sdf[[i]]$df[train_ind, ]
      test <- sdf[[i]]$df[-train_ind, ]
      rm(train_ind)

      # Fit the logistic regression model:
      mod <- glm(
        presence ~
          scale(ndvi) +
          scale(evi) +
          scale(lai) +
          scale(ptc) +
          scale(pnt) +
          scale(pnv) +
          scale(elev) +
          scale(water_dist) +
          scale(clay),
        data = train,
        family = binomial(link = "logit")
      )

      sdf[[i]]$lr_mod <- mod

      # Make predictions and create a confusion matrix:
      pred_lr <- predict(mod, newdata = test, type = "response")
      pred_class_lr <- ifelse(pred_lr >= 0.5, 1, 0) # Convert probabilities to class labels based on a threshold of 0.5
      pred_class_lr <- factor(pred_class_lr, levels = c(1, 0)) # Convert to factor with levels 1 and 0
      obs <- factor(test$presence, levels = c(1, 0)) # Convert observed values to factor with levels 1 and 0

      # Create the confusion matrix using the caret package:
      conf_mat <- caret::confusionMatrix(pred_class_lr, obs)

      # Store the model performance metrics in the data frame:
      model_performance[i, 1] <- "logistic_regression"
      model_performance[i, 2] <- sdf[[i]]$mod_id
      model_performance[i, 3] <- conf_mat$byClass["Sensitivity"]
      model_performance[i, 4] <- conf_mat$byClass["Specificity"]
      model_performance[i, 5] <- conf_mat$byClass["Sensitivity"] +
        conf_mat$byClass["Specificity"] -
        1
      model_performance[i, 6] <- auc(roc(test$presence, pred_lr, quiet = TRUE))
      model_performance[i, 7] <- conf_mat$overall["Kappa"]
      model_performance[i, 8] <- conf_mat$overall["Accuracy"]
      model_performance[i, 9] <- conf_mat$overall["AccuracyLower"]
      model_performance[i, 10] <- conf_mat$overall["AccuracyUpper"]
    }
  }

  rm(
    test,
    train,
    i,
    mod,
    pred_class_lr,
    pred_lr,
    conf_mat,
    obs
  )

  # --------------------------------------------------------------------------------------------------- #
  # --------------------------------------------------------------------------------------------------- #

  # RANDOM FOREST MODEL PERFORMANCE TRIAL: ####
  # ====================================== #

  # Tuning the model:
  # ----------------- #
  # Create a hyperparameter grid for tuning the random forest model:
  hyper_grid <- expand.grid(
    mtry = seq(2, 12, by = 1),
    # mtry = c(2, 3, 4, 5, 6, 7, 8, 10, 12),
    # node_size = seq(1, 15, by = 1),
    node_size = c(1, 2, 3, 4, 6, 8, 10, 15),
    # sample_size = seq(0.65, 0.95, by = 0.2),
    OOB_RMSE = 0,
    pred_error = 0
  )

  # Set up parallel processing:
  {
    num_cores <- detectCores() - 1 # Leave one core free to avoid overloading the system
    cl <- makeCluster(num_cores) # Create a cluster using the detected number of cores
    registerDoSNOW(cl) # Register the cluster for parallel processing
  }

  # Create progress bar
  {
    pb <- txtProgressBar(max = nrow(hyper_grid), style = 3)
    opts <- list(progress = function(n) {
      setTxtProgressBar(pb, n)
      if (n == nrow(hyper_grid)) cat("\n") # Add a newline after the progress bar is complete
    })
  }

  # system.time({
  {
    set.seed(3791)

    for (i in 1:length(sdf)) {
      # Print the progress of the model fitting and evaluation:
      cat("Completing model", i, "of", length(sdf), "\n")

      # ------------------------------------------------------------------------------------------------ #

      # Split the data into training and testing sets:
      {
        # Extract the data frame from the list and remove the geometry column:
        df <- as.data.frame(sdf[[i]]$df) |>
          dplyr::select(-geometry)
        train_ind <- createDataPartition(df$presence, p = 0.7, list = F)
        train <- df[train_ind, ]
        test <- df[-train_ind, ]
      }

      # ------------------------------------------------------------------------------------------------ #

      # Fit the random forest model with the current hyperparameters:
      results <- foreach(
        j = 1:nrow(hyper_grid), # Loop over each row of the hyperparameter grid
        .combine = rbind, # Combine results into a single data frame
        .packages = c("ranger"), # Load the ranger package in each parallel worker
        .export = c("train", "hyper_grid"), # Export the train data and hyper_grid to each worker
        .options.snow = opts, # Pass the progress bar options to the parallel workers
        .verbose = FALSE # Suppress warning messages
      ) %dopar%
        {
          # Fit the random forest model with the current hyperparameters:
          model <- ranger(
            formula = factor(presence, levels = c(0, 1)) ~ .,
            data = train,
            num.trees = 500,
            mtry = hyper_grid$mtry[j],
            min.node.size = hyper_grid$node_size[j],
            sample.fraction = 0.7,
            seed = 3791,
            num.threads = 1 # Ensure that each model runs on a single thread to prevent nested parallelism issues
          )

          # Store the OOB RMSE and prediction error in the hyperparameter grid:
          data.frame(
            OOB_RMSE = sqrt(model$prediction.error),
            pred_error = model$prediction.error * 100
          )
        }

      close(pb) # Close the progress bar after the loop is complete

      # ------------------------------------------------------------------------------------------------ #

      # Arrange the hyperparameter grid by prediction error
      hyper_params <- hyper_grid |>
        dplyr::arrange(pred_error)

      # Fit the final random forest model with the best hyperparameters:
      rf <- ranger(
        formula = factor(presence, levels = c(0, 1)) ~ .,
        data = train,
        num.trees = 500,
        mtry = hyper_params$mtry[1],
        min.node.size = hyper_params$node_size[1],
        sample.fraction = 0.7,
        seed = 3791,
        write.forest = T,
        importance = "impurity",
        probability = T
      )

      sdf[[i]]$rf_mod <- rf

      # ------------------------------------------------------------------------------------------------ #

      # EXTRACT THE MODEL PERFORMANCE METRICS: ####

      # Make predictions and create a confusion matrix:
      pred_rf <- predict(rf, data = test) # Get the predicted probabilities for the positive class (presence = 1)
      pred_rf <- pred_rf$predictions[, "1"] # Extract the predicted probabilities for the positive class (presence = 1)
      preds_rf <- ifelse(pred_rf >= 0.5, 1, 0) # Convert probabilities to class labels based on a threshold of 0.5
      pred_class_rf <- factor(preds_rf, levels = c(1, 0)) # Convert to factor with levels 1 and 0
      obs <- factor(test$presence, levels = c(1, 0)) # Convert observed values to factor with levels 1 and 0

      # Create the confusion matrix using the caret package:
      conf_mat <- confusionMatrix(pred_class_rf, obs)

      # Store the model performance metrics in the data frame:
      model_performance[i + n, 1] <- "Random_Forest"
      model_performance[i + n, 2] <- sdf[[i]]$mod_id
      model_performance[i + n, 3] <- conf_mat$byClass["Sensitivity"]
      model_performance[i + n, 4] <- conf_mat$byClass["Specificity"]
      model_performance[i + n, 5] <- conf_mat$byClass["Sensitivity"] +
        conf_mat$byClass["Specificity"] -
        1
      model_performance[i + n, 6] <- auc(roc(
        test$presence,
        pred_rf,
        quiet = TRUE
      ))
      model_performance[i + n, 7] <- conf_mat$overall["Kappa"]
      model_performance[i + n, 8] <- conf_mat$overall["Accuracy"]
      model_performance[i + n, 9] <- conf_mat$overall["AccuracyLower"]
      model_performance[i + n, 10] <- conf_mat$overall["AccuracyUpper"]
    } # End of loop

    # ------------------------------------------------------------------------------------------------ #

    stopCluster(cl)
  }
})


rm(
  df,
  hyper_grid,
  hyper_params,
  results,
  test,
  train,
  train_ind,
  cl,
  conf_mat,
  i,
  num_cores,
  obs,
  opts,
  pb,
  pred_class_rf,
  pred_rf,
  preds_rf,
  rf,
  n
)


# SAVE THE MODEL PERFORMANCE DATA AND THE LIST WITH THE MDOELS: ####
# ============================================================= #

write.csv(
  model_performance,
  "./Data/non_linear_200_data_model_performance.csv",
  row.names = FALSE
)

saveRDS(
  sdf,
  "./Data/non_linear_200_data_frames_with_models.rds"
)

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

model_perf <- model_performance |>
  dplyr::select(-model_id, -AccuracyLower, -AccuracyUpper) |>
  pivot_longer(cols = !Model_type, names_to = "type", values_to = "vals") |>
  dplyr::group_by(Model_type, type) |>
  dplyr::summarise(
    mean = mean(vals, na.rm = TRUE),
    sd = sd(vals, na.rm = TRUE)
  )

ggplot(model_perf, aes(x = type, y = mean, fill = Model_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(
    aes(ymax = mean + sd, ymin = mean - sd),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  # scale_fill_manual(
  # values = c("#535E31", "#A87B38"),
  # labels = c("Logistic Regression", "Random Forest")
  # ) +
  labs(fill = "", x = "Accuracy Measure", y = "Relative Value") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  # geom_point() +
  # geom_errorbar(aes(ymax = acc_up, ymin = acc_low), width = 0.1) +
  theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    axis.line = element_line(colour = "black", size = 0.8),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.text.x = element_text(
      size = 16,
      colour = "black",
      margin = ggplot2::margin(b = 15, t = 8)
    ),
    axis.title = element_text(size = 18),
    legend.key.height = unit(0.65, 'cm'),
    legend.key.width = unit(0.65, 'cm'),
    legend.key = element_rect(fill = "white", colour = "white"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.box = element_blank(),
    legend.position = "bottom"
  )


write.csv(
  model_performance,
  "./Data/non_linear_200_data_model_performance.csv",
  row.names = FALSE
)

mean(model_performance$Sensitivity[
  model_performance$Model_type == "logistic_regression"
])
mean(model_performance$Sensitivity[
  model_performance$Model_type == "Random_Forest"
])

mean(model_performance$Specificity[
  model_performance$Model_type == "logistic_regression"
])
mean(model_performance$Specificity[
  model_performance$Model_type == "Random_Forest"
])

mean(model_performance$TSS[
  model_performance$Model_type == "logistic_regression"
])
mean(model_performance$TSS[model_performance$Model_type == "Random_Forest"])

mean(model_performance$AUC[
  model_performance$Model_type == "logistic_regression"
])
mean(model_performance$AUC[model_performance$Model_type == "Random_Forest"])


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# Make a data frame with the equtions used to build the data frames:
equations <- data.frame(
  model_id = character(),
  beta_1 = numeric(),
  beta_2 = numeric(),
  beta_3 = numeric(),
  beta_4 = numeric(),
  beta_5 = numeric(),
  beta_6 = numeric(),
  c = numeric(),
  alpha = numeric()
)

for (i in 1:length(sdf)) {
  equations[i, 1] <- sdf[[i]]$mod_id
  equations[i, 2:9] <- sdf[[i]]$equation[1, ]
}

rm(i)


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

  # Convert hyper_grid to vectors for faster indexing:
  # beta_1_vec <- hyper_grid[, 1]
  # beta_2_vec <- hyper_grid[, 2]
  # beta_3_vec <- hyper_grid[, 3]
  # beta_4_vec <- hyper_grid[, 4]
  # beta_5_vec <- hyper_grid[, 5]
  # beta_6_vec <- hyper_grid[, 6]
  # c_vec <- hyper_grid[, 7]
  # alpha_vec <- hyper_grid[, 8]

  # b_1 <- beta_1_vec[i] / lai_sd
  #     b_2 <- beta_2_vec[i] / lai_sd
  #     b_3 <- beta_3_vec[i] / water_sd
  #     b_4 <- beta_4_vec[i] / ptc_sd
  #     b_5 <- beta_5_vec[i] / ndvi_sd
  #     b_6 <- beta_6_vec[i] / ndvi_sd

  test <- equations |>
    mutate(
      beta_1 = beta_1 / lai_sd,
      beta_2 = beta_2 / lai_sd,
      beta_3 = beta_3 / water_sd,
      beta_4 = beta_4 / ptc_sd,
      beta_5 = beta_5 / ndvi_sd,
      beta_6 = beta_6 / ndvi_sd
    )

  # Remove unnecessary objects
  rm(
    elev,
    evi,
    pnt,
    pnv,
    slope,
    soil,
    gnp
  )
  cat("Done!!!")
}


i <- 1

logit_p <- equations$alpha[i] +
  (equations$beta_1[i] * prob_rst$lai) +
  (equations$beta_2[i] * (prob_rst$lai^2)) +
  (equations$beta_3[i] * log(prob_rst$water_dist + 0.0001)) +
  ((equations$beta_4[i] * prob_rst$ptc) / (equations$c[i] + prob_rst$ptc)) +
  (equations$beta_5[i] * exp(equations$beta_6[i] * prob_rst$ndvi))

p <- 1 / (1 + exp(-logit_p))
p <- plogis(logit_p)
p <- terra::app(logit_p, fun = plogis)

plot(logit_p, col = inferno(150))
plot(p, col = inferno(150))

summary(logit_p)


((25002 / 60) / 60)


test <- read.csv("./Data/non_linear_200_data_model_performance.csv")


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# Make a data frame with the equtions used to build the data frames:
equations <- data.frame(
  model_id = character(),
  beta_1 = numeric(),
  beta_2 = numeric(),
  beta_3 = numeric(),
  beta_4 = numeric(),
  beta_5 = numeric(),
  beta_6 = numeric(),
  c = numeric(),
  alpha = numeric()
)

for (i in 1:length(sdf)) {
  equations[i, 1] <- sdf[[i]]$mod_id
  equations[i, 2:9] <- sdf[[i]]$equation[1, ]
}

rm(i)


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

  # Convert hyper_grid to vectors for faster indexing:
  # beta_1_vec <- hyper_grid[, 1]
  # beta_2_vec <- hyper_grid[, 2]
  # beta_3_vec <- hyper_grid[, 3]
  # beta_4_vec <- hyper_grid[, 4]
  # beta_5_vec <- hyper_grid[, 5]
  # beta_6_vec <- hyper_grid[, 6]
  # c_vec <- hyper_grid[, 7]
  # alpha_vec <- hyper_grid[, 8]

  # b_1 <- beta_1_vec[i] / lai_sd
  #     b_2 <- beta_2_vec[i] / lai_sd
  #     b_3 <- beta_3_vec[i] / water_sd
  #     b_4 <- beta_4_vec[i] / ptc_sd
  #     b_5 <- beta_5_vec[i] / ndvi_sd
  #     b_6 <- beta_6_vec[i] / ndvi_sd

  test <- equations |>
    mutate(
      beta_1 = beta_1 / lai_sd,
      beta_2 = beta_2 / lai_sd,
      beta_3 = beta_3 / water_sd,
      beta_4 = beta_4 / ptc_sd,
      beta_5 = beta_5 / ndvi_sd,
      beta_6 = beta_6 / ndvi_sd
    )

  # Remove unnecessary objects
  rm(
    elev,
    evi,
    pnt,
    pnv,
    slope,
    soil,
    gnp
  )
  cat("Done!!!")
}


i <- 1

logit_p <- equations$alpha[i] +
  (equations$beta_1[i] * prob_rst$lai) +
  (equations$beta_2[i] * (prob_rst$lai^2)) +
  (equations$beta_3[i] * log(prob_rst$water_dist + 0.0001)) +
  ((equations$beta_4[i] * prob_rst$ptc) / (equations$c[i] + prob_rst$ptc)) +
  (equations$beta_5[i] * exp(equations$beta_6[i] * prob_rst$ndvi))

p <- 1 / (1 + exp(-logit_p))
p <- plogis(logit_p)
p <- terra::app(logit_p, fun = plogis)

plot(logit_p, col = inferno(150))
plot(p, col = inferno(150))

summary(logit_p)


((25002 / 60) / 60)
