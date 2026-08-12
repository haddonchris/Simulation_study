rm(list = ls())
graphics.off()


# LOAD LIBRARIES: ####
# =============== #

{
  library(terra)
  library(sf)
  library(viridis)
  library(ranger)
  library(tidyverse)
}

# READ IN THE DATA: ####
# ================= #

sdf <- readRDS("./Data/Final_data/Non_linear_200_data_frames_with_models.rds")

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
  # ---------------------------------------- #

  cat("Resampling rasters...", "\n")
  rst <- terra::rast("./Satellite_images/Elevation/elevation_export.tif")

  ndvi <- resample(x = ndvi_mean, y = rst)
  lai <- resample(x = lai_mean, y = rst)
  ptc <- resample(x = ptc_mean, y = rst)
  water_dist <- resample(x = water_dist, y = rst)

  rm(
    rst,
    ndvi_mean,
    lai_mean,
    ptc_mean
  )

  # Load park boundary:
  # ------------------- #

  gnp <- st_read(
    "./Gorongosa_files/GNP_shp/WDPA_WDOECM_Sep2025_Public_801_shp-polygons.shp",
    quiet = T
  )

  # Edit/Clean the data
  gnp <- gnp[1, 31]
  gnp <- st_transform(gnp, "epsg:32736") # Set the CRS

  # Crop and mask the stack:
  # ------------------------ #

  cat("Masking and cropping rasters...", "\n")
  ndvi <- mask(crop(ndvi, gnp), gnp)
  lai <- mask(crop(lai, gnp), gnp)
  ptc <- mask(crop(ptc, gnp), gnp)
  water_dist <- mask(crop(water_dist, gnp), gnp)

  stack <- c(
    ndvi,
    lai,
    ptc,
    water_dist
  )

  cat("Calculating raster statistics...", "\n")
  ndvi_stats <- global(ndvi, fun = c("mean", "sd"), na.rm = TRUE)
  water_stats <- global(water_dist, fun = c("mean", "sd"), na.rm = TRUE)
  ptc_stats <- global(ptc, fun = c("mean", "sd"), na.rm = TRUE)
  lai_stats <- global(lai, fun = c("mean", "sd"), na.rm = TRUE)

  # Remove unnecessary objects
  rm(
    lai,
    ndvi,
    ptc,
    water_dist
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
  cellsize = 5000, # 5 km
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
  x = numeric(),
  y = numeric(),
  true_p = numeric(),
  lr_p = numeric(),
  rf_p = numeric()
)

# Make a vector of the sampling points to extract the probabilities from the rasters:
vect_p_true <- terra::vect(pts, geom = c("x", "y"), crs = "epsg:32736")

# Define a function to extract the probabilities from the random forest model:
rf_fun <- function(model, data, ...) {
  p <- predict(model, data = data, type = "response")$predictions

  if (is.matrix(p) || is.data.frame(p)) {
    p <- p[, "1"] # replace with the presence-class name if needed
  }
  as.numeric(p)
}

system.time({
  # Start time of the loop
  start_time <- Sys.time()

  for (i in 1:length(sdf)) {
    # Start time for the current iteration
    iteration_start <- Sys.time()

    # -------------------------------------------------------------------------- #

    # Make a vector containing the model ID for the current iteration:
    mod_id <- sdf[[i]]$mod_id

    # -------------------------------------------------------------------------- #

    # CALCULATE THE TRUE PROBABILITIES FOR THE CURRENT MODEL:
    # ------------------------------------------------------- #

    # Scale the coefficients by the standard deviation of each raster:
    b_ndvi <- sdf[[i]]$equation$beta_ndvi / ndvi_stats$sd
    b_water <- sdf[[i]]$equation$beta_water / water_stats$sd
    b_ptc <- sdf[[i]]$equation$beta_ptc / ptc_stats$sd
    b_lai <- sdf[[i]]$equation$beta_lai / lai_stats$sd

    logit_rst <- sdf[[i]]$equation$alpha +
      (b_ndvi *
        exp(
          -(stack$ndvi - sdf[[i]]$equation$mu_ndvi)^2 /
            (2 * sdf[[i]]$equation$sigma_ndvi^2)
        )) +
      (b_water *
        exp(
          -sdf[[i]]$equation$k_water * stack$water_dist
        )) +
      (b_ptc *
        (1 /
          (1 +
            exp(
              -sdf[[i]]$equation$k_ptc *
                (stack$ptc - sdf[[i]]$equation$c_ptc)
            )))) +
      (b_lai * (-((stack$lai - sdf[[i]]$equation$mu_lai)^2)))

    p_true_rst <- 1 / (1 + exp(-logit_rst))

    names(p_true_rst) <- "true_p"
    varnames(p_true_rst) <- "true_p"

    # Save True probability maps as PNG files:
    {
      png(
        filename = paste0(
          "./Data/Figures/Non_linear_predictive_maps/True/NL_True_surface_",
          i,
          ".png"
        ),
        width = 2000,
        height = 2200,
        res = 300
      )

      plot(
        p_true_rst,
        main = paste0(
          "Simulated Non-Linear Data Frame ",
          i,
          "\nTrue Probabilities"
        ),
        col = viridis(250)
      )

      dev.off()
    }

    var <- terra::extract(p_true_rst, vect_p_true, ID = FALSE)

    true_p <- var$true_p

    rm(
      logit_rst,
      p_true_rst
    )

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

    # Save Logistic regression predicted probability maps as PNG files:
    {
      png(
        filename = paste0(
          "./Data/Figures/Non_linear_predictive_maps/LR/NL_LR_prob_surface_",
          i,
          ".png"
        ),
        width = 2000,
        height = 2200,
        res = 300
      )

      plot(
        lr_rst,
        main = paste0(
          "Simulated Non-Linear Data Frame ",
          i,
          "\nLinear Regression Probabilities"
        ),
        col = viridis(250)
      )

      dev.off()
    }

    var <- terra::extract(lr_rst, vect_p_true, ID = FALSE)

    lr_p <- var$lr_p

    rm(
      lr_rst
    )

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

    na_mask <- is.na(stack[[1]]) # Create a mask for the NA values in the first layer of the stack
    rf_rst[na_mask] <- NA

    names(rf_rst) <- "rf_p"
    varnames(rf_rst) <- "rf_p"

    # Save Random forest predicted probability maps as PNG files:
    {
      png(
        filename = paste0(
          "./Data/Figures/Non_linear_predictive_maps/RF/NL_RF_prob_surface_",
          i,
          ".png"
        ),
        width = 2000,
        height = 2200,
        res = 300
      )

      plot(
        rf_rst,
        main = paste0("Simulated Non-Linear Data Frame ", i, "\nRandom Forest"),
        col = viridis(250)
      )

      dev.off()
    }

    var <- terra::extract(rf_rst, vect_p_true, ID = FALSE)

    rf_p <- var$rf_p

    rm(
      rf_rst,
      na_mask
    )

    # -------------------------------------------------------------------------- #

    # ADD THE PROBABILITIES TO THE DATA FRAME:
    # --------------------------------------- #

    temp <- data.frame(
      mod_id = mod_id,
      x = pts$x,
      y = pts$y,
      true_p = true_p,
      lr_p = lr_p,
      rf_p = rf_p
    )

    probs <- rbind(probs, temp)

    # -------------------------------------------------------------------------- #

    # GARAGE COLLECT AND REMOVE UNNECESSARY OBJECTS:
    # ---------------------------------------------- #

    if (i %% 5 == 0) {
      gc()
    }

    # -------------------------------------------------------------------------- #

    # TRACK PROGRESS OF MODELS:
    # ------------------------- #

    # Time taken for current iteration
    iteration_time <- as.numeric(difftime(
      Sys.time(),
      iteration_start,
      units = "mins"
    ))

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

    avg_time <- elapsed / i

    remaining <- avg_time * (length(sdf) - i)

    estimated_end <- Sys.time() + (remaining * 60)

    percent <- (i / length(sdf)) * 100

    # Print progress on the same line
    cat(
      "\r",
      sprintf(
        "Progress: %5.1f%% | Model: %d/%d | Iteration: %.1fm | Elapsed: %.1fm | Remaining: %.1fm | ETA: %s",
        percent,
        i,
        length(sdf),
        iteration_time,
        elapsed,
        remaining,
        format(estimated_end, "%H:%M:%S")
      ),
      sep = ""
    )

    # Flush the console to update the progress bar
    flush.console()
  } # End of the for loop

  # -------------------------------------------------------------------------- #

  # Clear objects from the environment to free up memory
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
    vect_p_true,
    i,
    true_p,
    lr_p,
    rf_p,
    mod_id,
    avg_time,
    elapsed,
    iteration_time,
    remaining,
    estimated_end,
    iteration_start
  )

  # -------------------------------------------------------------------------- #

  # SAVE THE RAW PROBABILITIES TO A CSV FILE:
  # ========================================= #

  write.csv(
    probs,
    "./Data/Final_data/Non_linear_200_probabilities.csv",
    row.names = FALSE
  )

  # -------------------------------------------------------------------------- #

  # CALCULATE THE METRICS FOR EACH MODEL: ####
  # ==================================== #

  # Calculate the metrics for the Random Forest models
  rf_metrics <- probs |>
    dplyr::select(-c(x, y, lr_p)) |>
    mutate(rf_error = true_p - rf_p) |>
    dplyr::relocate(mod_id, rf_error, true_p, rf_p) |>
    dplyr::group_by(mod_id) |>
    dplyr::summarise(
      sim = "NL",
      mod_type = "Random Forest",
      rmse = sqrt(mean(abs(rf_error)^2)),
      mae = mean(abs(rf_error)),
      bias = mean(rf_error),
      cor = cor(rf_p, true_p)
    ) |>
    dplyr::ungroup()

  # Calculate the metrics for the Logistic Regression models
  lr_metrics <- probs |>
    dplyr::select(-c(x, y, rf_p)) |>
    mutate(lr_error = true_p - lr_p) |>
    dplyr::relocate(mod_id, lr_error, true_p, lr_p) |>
    dplyr::group_by(mod_id) |>
    dplyr::summarise(
      sim = "NL",
      mod_type = "Logistic Regression",
      rmse = sqrt(mean(abs(lr_error)^2)),
      mae = mean(abs(lr_error)),
      bias = mean(lr_error),
      cor = cor(lr_p, true_p)
    ) |>
    dplyr::ungroup()

  # Combine the metrics for both models into a single data frame
  metrics <- rbind(rf_metrics, lr_metrics)

  # Save the metrics to a CSV file
  write.csv(
    metrics,
    "./Data/Final_data/Non_linear_200_model_metrics.csv",
    row.names = FALSE
  )

  # Remove the unnecessary objects
  rm(rf_metrics, lr_metrics)

  # -------------------------------------------------------------------------- #

  # TIME METRICS FOR THE LOOP:
  # ========================== #

  # End time
  end_time <- Sys.time()

  # Calculate the total time taken for the loop
  total_time <- difftime(
    end_time,
    start_time,
    units = "hours"
  )

  # To be printed in the console after the loop is complete
  cat("\n\n")
  cat("========================================\n")
  cat("Simulation complete\n")
  cat("========================================\n")
  cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("End time:  ", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total time:", round(as.numeric(total_time), 2), "hours\n")
  cat("========================================\n")

  rm(start_time, end_time, total_time)
})
