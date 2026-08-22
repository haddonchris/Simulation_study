rm(list = ls())
graphics.off()

# LOAD LIBRARIES: ####
# =============== #

{
  library(tidyverse)
  library(terra)
  library(sf)
  library(sp)
  library(raster)
  library(viridis)
  library(ranger)
  library(sdmTMB)
  library(adehabitatHR)
}

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# LOAD IN THE GEOGRAPHIC & SATELLITE DATA: ####
# ======================================== #

{
  # NDVI:
  {
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
  }

  # Leaf Area Index:
  {
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
  }

  # Percent Canopy Cover:
  {
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
  }

  # Distance to water:
  {
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
  }

  # Resample rasters to the same resolution:
  {
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
  }

  # Load park boundary:
  {
    gnp <- st_read(
      "./Gorongosa_files/GNP_shp/WDPA_WDOECM_Sep2025_Public_801_shp-polygons.shp",
      quiet = T
    )

    # Edit/Clean the data
    gnp <- gnp[1, 31]
    gnp <- st_transform(gnp, "epsg:32736") # Set the CRS
  }

  # Crop and mask the stack:
  {
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

    # Calculate the mean and standard deviation for the rasters (for scaling the coefficients)
    cat("Calculating raster statistics...", "\n")
    ndvi_stats <- global(
      ndvi,
      fun = c("mean", "sd", "min", "max"),
      na.rm = TRUE
    )
    ptc_stats <- global(ptc, fun = c("mean", "sd", "min", "max"), na.rm = TRUE)
    lai_stats <- global(lai, fun = c("mean", "sd", "min", "max"), na.rm = TRUE)
  }

  cat("Done!!!")
}


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# GENERATE SIMULATED DATA: ####
# ======================== #

# Create an empty list for all of the data frames to go into
data_frames <- list()

# Create the candidate points for the simulation:
{
  set.seed(3791)

  # Generate random points within GNP:
  cand_points <- st_sample(
    gnp,
    size = 20000,
    type = "random"
  )

  # Convert to sf object:
  cand_points <- st_sf(geometry = cand_points)

  # Extract coordinates:
  coords <- st_coordinates(cand_points)

  # Extract covariate values at candidate points:
  cand_points$NDVI <- terra::extract(ndvi, coords)[, 1]
  cand_points$ptc <- terra::extract(ptc, coords)[, 1]
  cand_points$lai <- terra::extract(lai, coords)[, 1]
  cand_points$dist_water <- terra::extract(water_dist, coords)[, 1]

  ndvi_q20_80 <- quantile(cand_points$NDVI, probs = c(0.20, 0.80), na.rm = TRUE)
  lai_q20_80 <- quantile(cand_points$lai, probs = c(0.20, 0.80), na.rm = TRUE)
}

# Create a data frame to store diagnostics for each simulation:
diagnostics <- data.frame(
  sim = seq(1, 2000000),
  eta_sd = NA_real_,
  p_mean = NA_real_,
  p_sd = NA_real_,
  p_min = NA_real_,
  p_max = NA_real_,
  p_low = NA_real_,
  p_high = NA_real_,
  n_pres = NA_integer_
)

# Create a vector to store the indices of the hyper_grid that were kept for the 200 data frames
diag_index <- NA

system.time({
  {
    set.seed(3791)

    # Create a hyper grid with all of the parameters for the simulated data frames
    hyper_grid <- data.frame(
      beta_ndvi = runif(2000000, -2, 2),
      beta_water = runif(2000000, -2, 2),
      beta_ptc = runif(2000000, -2, 2),
      beta_lai = runif(2000000, -2, 2),
      alpha = runif(2000000, -2, 2),
      mu_ndvi = runif(2000000, min = ndvi_q20_80[1], max = ndvi_q20_80[2]),
      sigma_ndvi = runif(2000000, 0.10, 0.50) *
        (ndvi_stats$max - ndvi_stats$min),
      mu_lai = runif(2000000, min = lai_q20_80[1], max = lai_q20_80[2]),
      r_lai = max(lai_stats$max - lai_stats$min),
      k_water = log(2) / runif(2000000, 100, 17000),
      k_ptc = 6 / runif(2000000, 5, 50),
      c_ptc = runif(2000000, ptc_stats$min, ptc_stats$max)
    )

    # Create an empty list for all of the data frames to go into
    data_frames <- list()
    kept <- 0

    # Convert hyper_grid to vectors for faster indexing:
    beta_ndvi_vec <- hyper_grid[, 1]
    beta_water_vec <- hyper_grid[, 2]
    beta_ptc_vec <- hyper_grid[, 3]
    beta_lai_vec <- hyper_grid[, 4]
    alpha_vec <- hyper_grid[, 5]
    mu_ndvi_vec <- hyper_grid[, 6]
    sigma_ndvi_vec <- hyper_grid[, 7]
    mu_lai_vec <- hyper_grid[, 8]
    r_lai_vec <- hyper_grid[, 9]
    k_water_vec <- hyper_grid[, 10]
    k_ptc_vec <- hyper_grid[, 11]
    c_ptc_vec <- hyper_grid[, 12]

    # Extract candidate points data as vectors for faster computation:
    cand_ndvi <- cand_points$NDVI
    cand_water <- cand_points$dist_water
    cand_ptc <- cand_points$ptc
    cand_lai <- cand_points$lai
    cand_geom <- cand_points$geometry
    n_cand <- length(cand_geom)

    # Create a progress bar to track the progress of the simulation
    pb <- txtProgressBar(
      min = 0,
      max = 200,
      style = 3
    )

    # Simulation For Loop (vectorized):
    for (i in 1:nrow(hyper_grid)) {
      # Calculate the non-linear transformations for each covariate (vectorized):
      g_ndvi <- exp(
        -(cand_ndvi - mu_ndvi_vec[i])^2 /
          (2 * sigma_ndvi_vec[i]^2)
      )
      g_lai <- -(((cand_lai - mu_lai_vec[i]) / r_lai_vec[i])^2)

      g_water <- exp(-k_water_vec[i] * cand_water)

      g_ptc <- 1 /
        (1 + exp(-k_ptc_vec[i] * (cand_ptc - c_ptc_vec[i])))

      # DEFINE THE EQUATION AND CALCULATE PROBABILITIES:
      logit_p <- alpha_vec[i] +
        beta_ndvi_vec[i] * g_ndvi +
        beta_water_vec[i] * g_water +
        beta_ptc_vec[i] * g_ptc +
        beta_lai_vec[i] * g_lai
      # Convert logit to probability:
      p_true <- plogis(logit_p)

      # Bernoulli trial to simulate presence/absence:
      presence <- rbinom(n = n_cand, size = 1, prob = p_true)

      # Get indices of presences (faster than filtering):
      presence_idx <- which(presence == 1)
      n_presences <- length(presence_idx)

      diagnostics$n_pres[i] <- n_presences
      diagnostics$eta_sd[i] <- sd(logit_p)
      diagnostics$p_mean[i] <- mean(p_true)
      diagnostics$p_sd[i] <- sd(p_true)
      diagnostics$p_min[i] <- min(p_true)
      diagnostics$p_max[i] <- max(p_true)
      diagnostics$p_low[i] <- mean(p_true < 0.10)
      diagnostics$p_high[i] <- mean(p_true > 0.90)

      # --------------------------------------------------------------------------------------------------- #

      # FILTERING CONDITIONS:
      # ===================== #

      if (diagnostics$eta_sd[i] < 0.61 || diagnostics$eta_sd[i] > 1.0) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_mean[i] < 0.2 || diagnostics$p_mean[i] > 0.4) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_sd[i] < 0.13 || diagnostics$p_sd[i] > 0.2) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_min[i] > 0.15) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_max[i] < 0.85) {
        next # Skip to the next iteration if the condition is met
      }

      # --------------------------------------------------------------------------------------------------- #

      # GENERATE ABSENCE POINTS & COMBINE DATA SETS:

      # Only proceed if presences are between 1000 and 10 000:
      if (n_presences >= 1000 && n_presences <= 10000) {
        # Keep track of how many data frames have been created
        kept <- kept + 1
        diag_index[kept] <- i

        # Extract presence geometries and create sf object:
        presences_geom <- cand_geom[presence_idx]
        presences_df <- st_sf(geometry = presences_geom, presence = 1)

        # Generate absence points:
        absence_points <- st_sample(gnp, size = n_presences, type = "random")
        absence_df <- st_sf(geometry = absence_points, presence = 0)

        # Combine presences and absences:
        sim_data <- rbind(presences_df, absence_df)

        data_frames[[kept]] <- list(
          df = sim_data,
          equation = hyper_grid[i, ]
        )

        # Create a raster of the true probabilities for visualization:
        # ------------------------------------------------------------ #

        f_ndvi <- exp(
          -(stack$ndvi - mu_ndvi_vec[i])^2 /
            (2 * sigma_ndvi_vec[i]^2)
        )
        f_lai <- -(((stack$lai - mu_lai_vec[i]) / r_lai_vec[i])^2)

        f_water <- exp(-k_water_vec[i] * stack$water_dist)

        f_ptc <- 1 /
          (1 + exp(-k_ptc_vec[i] * (stack$ptc - c_ptc_vec[i])))

        # Create probability raster:
        logit_rst <- alpha_vec[i] +
          (beta_ndvi_vec[i] * f_ndvi) +
          (beta_water_vec[i] * f_water) +
          (beta_ptc_vec[i] * f_ptc) +
          (beta_lai_vec[i] * f_lai)

        p_true_rst <- 1 / (1 + exp(-logit_rst))

        png(
          filename = paste0(
            "./Data/Cal_images/Non_linear_calibration/test3/Simulated_data_frame_",
            kept,
            ".png"
          ),
          width = 2000,
          height = 2200,
          res = 300
        )

        plot(
          p_true_rst,
          main = paste0("Simulated Data Frame ", kept, "\nIteration ", i),
          col = viridis(250)
        )

        dev.off()

        # Update the progress bar
        setTxtProgressBar(pb, kept)
      } else {
        next # Skip to the next iteration if the condition is not met
      }

      # Keep track of Loop process
      if (kept == 200) break
    }

    close(pb) # Close the progress bar

    # Remove unnecessary objects
    rm(
      b_1,
      b_2,
      b_3,
      b_4,
      b_5,
      b_6,
      beta_1_vec,
      beta_2_vec,
      beta_3_vec,
      beta_4_vec,
      beta_5_vec,
      beta_6_vec,
      c_vec,
      cand_points,
      coords,
      i,
      presence_idx,
      presence,
      p_true,
      logit_p,
      absence_df,
      presences_df,
      presences_geom,
      sim_data,
      absence_points,
      cand_ndvi,
      cand_water,
      cand_ptc,
      cand_geom,
      n_cand,
      n_presences,
      ndvi_mean,
      ndvi_sd,
      water_mean,
      water_sd,
      ptc_mean,
      ptc_sd,
      alpha_vec,
      hyper_grid,
      gnp,
      ndvi_stats,
      ptc_stats,
      water_stats,
      ndvi,
      ptc,
      water_dist,
      lai_stats,
      cand_lai,
      lai,
      lai_mean,
      lai_sd,
      kept
    )

    # FILTER THE LISTS AND EDIT: ####
    # ========================== #

    # Add a unique ID and extract all satellite variables for each data frame (vectorized):
    for (i in 1:length(data_frames)) {
      # Add unique ID:
      data_frames[[i]]$mod_id <- paste0("mod_", i)

      # Extract satellite imagery data (vectorized):
      coords <- st_coordinates(data_frames[[i]]$df)
      vals <- terra::extract(stack, terra::vect(coords), df = TRUE)
      data_frames[[i]]$df <- cbind(data_frames[[i]]$df, vals[, -1])

      if (i %% 10 == 0) {
        cat("Extracting for data frame", i, "of", length(data_frames), "\n")
      }
    }

    # Remove unnecessary objects:
    rm(coords, i, stack, vals)
  }
})

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# SAVE THE DATA FRAMES: ####
# ===================== #

saveRDS(
  data_frames,
  "./Data/Final_data/Non_linear_data.rds"
)
