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
    ndvi_stats <- global(ndvi, fun = c("mean", "sd"), na.rm = TRUE)
    water_stats <- global(water_dist, fun = c("mean", "sd"), na.rm = TRUE)
    ptc_stats <- global(ptc, fun = c("mean", "sd"), na.rm = TRUE)
    lai_stats <- global(lai, fun = c("mean", "sd"), na.rm = TRUE)

    ndvi_range <- terra::minmax(ndvi, compute = TRUE)
    water_range <- terra::minmax(water_dist, compute = TRUE)
    ptc_range <- terra::minmax(ptc, compute = TRUE)
    lai_range <- terra::minmax(lai, compute = TRUE)
  }

  cat("Done!!!")
}

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# GENERATE SIMULATED DATA: ####
# ======================== #

# Function to sample z values excluding a specific range (e.g., around 1):
sample_z_excluding_one <- function(
  n,
  lower = 0.1,
  upper = 4,
  excl_low = 0.5,
  excl_high = 1.5
) {
  # Two allowed segments: [lower, excl_low) and (excl_high, upper]
  w1 <- excl_low - lower
  w2 <- upper - excl_high
  total_w <- w1 + w2

  # Split n draws between the two segments proportional to their width,
  # so the result is still uniform density over the *allowed* region
  n1 <- rbinom(1, n, w1 / total_w)
  n2 <- n - n1

  z <- c(
    runif(n1, lower, excl_low),
    runif(n2, excl_high, upper)
  )
  sample(z) # shuffle so segment membership isn't tied to row position
}

# Create a hyper grid with all of the parameters for the simulated data frames
{
  set.seed(3791)
  hyper_grid <- data.frame(
    beta_ndvi = runif(1000000, -2, 2),
    beta_water = runif(1000000, -2, 2),
    beta_ptc = runif(1000000, -2, 2),
    beta_lai = runif(1000000, -2, 2),
    alpha = runif(1000000, -2, 2),
    z_ndvi = sample_z_excluding_one(1000000),
    z_water = sample_z_excluding_one(1000000),
    z_ptc = sample_z_excluding_one(1000000),
    z_lai = sample_z_excluding_one(1000000)
  )
}

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
}

# Create a diagnostics data frame to store the statistics for each simulated data frame:
diagnostics <- data.frame(
  sim = seq(1, nrow(hyper_grid)),
  eta_sd = NA_real_,
  p_mean = NA_real_,
  p_sd = NA_real_,
  p_min = NA_real_,
  p_max = NA_real_,
  n_pres = NA_integer_
)

diag_index <- NA

system.time({
  {
    set.seed(3791)

    # Randomise the hyper grid and sample 1,000,000 rows for the simulation:
    hyper_grid <- hyper_grid[sample(nrow(hyper_grid), 1000000), ]

    # Create an empty list for all of the data frames to go into
    data_frames <- list()
    kept <- 0

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
    beta_ndvi_vec <- hyper_grid[, 1]
    beta_water_vec <- hyper_grid[, 2]
    beta_ptc_vec <- hyper_grid[, 3]
    beta_lai_vec <- hyper_grid[, 4]
    alpha_vec <- hyper_grid[, 5]
    z_ndvi_vec <- hyper_grid[, 6]
    z_water_vec <- hyper_grid[, 7]
    z_ptc_vec <- hyper_grid[, 8]
    z_lai_vec <- hyper_grid[, 9]

    # Extract candidate points data as vectors for faster computation:
    cand_ndvi <- cand_points$NDVI
    cand_water <- cand_points$dist_water
    cand_ptc <- cand_points$ptc
    cand_lai <- cand_points$lai
    cand_geom <- cand_points$geometry
    n_cand <- length(cand_geom)

    cand_ndvi_norm <- (cand_ndvi - ndvi_range["min", 1]) /
      (ndvi_range["max", 1] - ndvi_range["min", 1])
    cand_water_norm <- (cand_water - water_range["min", 1]) /
      (water_range["max", 1] - water_range["min", 1])
    cand_ptc_norm <- (cand_ptc - ptc_range["min", 1]) /
      (ptc_range["max", 1] - ptc_range["min", 1])
    cand_lai_norm <- (cand_lai - lai_range["min", 1]) /
      (lai_range["max", 1] - lai_range["min", 1])

    # Create a progress bar to track the progress of the simulation
    pb <- txtProgressBar(
      min = 0,
      max = 200,
      style = 3
    )

    # Simulation For Loop (vectorized):
    for (i in 1:nrow(hyper_grid)) {
      # CALCULATE SCALED COEFFICIENTS (vectorized):
      ndvi_term <- cand_ndvi_norm^z_ndvi_vec[i]
      water_term <- cand_water_norm^z_water_vec[i]
      ptc_term <- cand_ptc_norm^z_ptc_vec[i]
      lai_term <- cand_lai_norm^z_lai_vec[i]

      b_ndvi <- beta_ndvi_vec[i] / sd(ndvi_term)
      b_water <- beta_water_vec[i] / sd(water_term)
      b_ptc <- beta_ptc_vec[i] / sd(ptc_term)
      b_lai <- beta_lai_vec[i] / sd(lai_term)

      logit_p <- alpha_vec[i] +
        b_ndvi * ndvi_term +
        b_water * water_term +
        b_ptc * ptc_term +
        b_lai * lai_term

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

      # --------------------------------------------------------------------------------------------------- #

      # FILTERING CONDITIONS:
      # ===================== #

      if (diagnostics$eta_sd[i] < 1.634 || diagnostics$eta_sd[i] > 2.5) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_mean[i] < 0.21 || diagnostics$p_mean[i] > 0.5) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_sd[i] < 0.24 || diagnostics$p_sd[i] > 0.33) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_min[i] > 0.10) {
        next # Skip to the next iteration if the condition is met
      }

      if (diagnostics$p_max[i] < 0.90) {
        next # Skip to the next iteration if the condition is met
      }

      # --------------------------------------------------------------------------------------------------- #

      # GENERATE ABSENCE POINTS & COMBINE DATA SETS:

      # Only proceed if presences are between 1000 and 10000:
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
          equation = data.frame(
            alpha = alpha_vec[i],
            b_ndvi = b_ndvi,
            b_water = b_water,
            b_ptc = b_ptc,
            b_lai = b_lai,
            z_ndvi = z_ndvi_vec[i],
            z_water = z_water_vec[i],
            z_ptc = z_ptc_vec[i],
            z_lai = z_lai_vec[i]
          )
        )

        # Create a raster of the true probabilities for visualization:
        # ------------------------------------------------------------ #

        stack_ndvi_norm <- (stack$ndvi - ndvi_range["min", 1]) /
          (ndvi_range["max", 1] - ndvi_range["min", 1])
        stack_water_norm <- (stack$water_dist - water_range["min", 1]) /
          (water_range["max", 1] - water_range["min", 1])
        stack_ptc_norm <- (stack$ptc - ptc_range["min", 1]) /
          (ptc_range["max", 1] - ptc_range["min", 1])
        stack_lai_norm <- (stack$lai - lai_range["min", 1]) /
          (lai_range["max", 1] - lai_range["min", 1])

        logit_rst <- alpha_vec[i] +
          (b_ndvi * (stack_ndvi_norm^z_ndvi_vec[i])) +
          (b_water * (stack_water_norm^z_water_vec[i])) +
          (b_ptc * (stack_ptc_norm^z_ptc_vec[i])) +
          (b_lai * (stack_lai_norm^z_lai_vec[i]))

        p_true_rst <- 1 / (1 + exp(-logit_rst))

        # Create a raster of the true probabilities for visualization
        png(
          filename = paste0(
            "./Data/Cal_images/Semi_non_linear_calibration/Simulated_data_frame_",
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

      # Break the loop if 200 data frames have been created
      if (kept == 200) break
    }

    close(pb) # Close the progress bar

    # Remove unnecessary objects
    rm(
      cand_points,
      coords,
      i,
      presence_idx,
      presence,
      p_true,
      logit_p,
      b_ndvi,
      b_water,
      b_ptc,
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
      beta_ndvi_vec,
      beta_water_vec,
      beta_ptc_vec,
      alpha_vec,
      hyper_grid,
      gnp,
      ndvi_stats,
      ptc_stats,
      water_stats,
      ndvi,
      ptc,
      water_dist,
      b_lai,
      lai_stats,
      beta_lai_vec,
      cand_lai,
      lai,
      lai_mean,
      lai_sd,
      z_lai_vec,
      z_ndvi_vec,
      z_water_vec,
      z_ptc_vec,
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

      # Keep track of the loop's progress
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
  "./Data/Final_data/Semi_non_linear_data.rds"
)
