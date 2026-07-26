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

  # EVI
  {
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

  # Percent non-tree vegetation:
  {
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
  }

  # Percent non-vegetation
  {
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
  }

  # Terrain variables
  {
    cat("Getting terrain data...", "\n")

    elev <- terra::rast("./Satellite_images/Elevation/elevation_export.tif")
    names(elev) <- "elev"
    slope <- terra::rast("./Satellite_images/Slope/slope_export.tif")
    names(slope) <- "slope"
  }

  # Soil Variables:
  {
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

    # Calculate the mean and standard deviation for the rasters (for scaling the coefficients)
    cat("Calculating raster statistics...", "\n")
    ndvi_stats <- global(ndvi, fun = c("mean", "sd"), na.rm = TRUE)
    water_stats <- global(water_dist, fun = c("mean", "sd"), na.rm = TRUE)
    ptc_stats <- global(ptc, fun = c("mean", "sd"), na.rm = TRUE)
    lai_stats <- global(lai, fun = c("mean", "sd"), na.rm = TRUE)
  }

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

# GENERATE SIMULATED DATA: ####
# ======================== #

# Create a hyper grid with all of the parameters for the simulated data frames
hyper_grid <- expand.grid(
  beta_ndvi = setdiff(seq(from = -2, to = 2, by = 0.5), 0),
  beta_water = setdiff(seq(from = -2, to = 2, by = 0.5), 0),
  beta_ptc = setdiff(seq(from = -2, to = 2, by = 0.5), 0),
  beta_lai = setdiff(seq(from = -2, to = 2, by = 0.5), 0),
  alpha = setdiff(seq(from = -2, to = 2, by = 0.5), 0)
)

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
}

system.time({
  {
    set.seed(3791)

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

    # Extract candidate points data as vectors for faster computation:
    cand_ndvi <- cand_points$NDVI
    cand_water <- cand_points$dist_water
    cand_ptc <- cand_points$ptc
    cand_lai <- cand_points$lai
    cand_geom <- cand_points$geometry
    n_cand <- length(cand_geom)

    # Simulation For Loop (vectorized):
    for (i in 1:nrow(hyper_grid)) {
      # CALCULATE SCALED COEFFICIENTS (vectorized):
      b_ndvi <- beta_ndvi_vec[i] / ndvi_sd
      b_water <- beta_water_vec[i] / water_sd
      b_ptc <- beta_ptc_vec[i] / ptc_sd
      b_lai <- beta_lai_vec[i] / lai_sd # Using LAI stats for scaling LAI

      # DEFINE THE EQUATION AND CALCULATE PROBABILITIES:
      logit_p <- alpha_vec[i] +
        b_ndvi * cand_ndvi +
        b_water * cand_water +
        b_ptc * cand_ptc +
        b_lai * cand_lai

      # Convert logit to probability:
      p_true <- plogis(logit_p)

      # Bernoulli trial to simulate presence/absence:
      presence <- rbinom(n = n_cand, size = 1, prob = p_true)

      # Get indices of presences (faster than filtering):
      presence_idx <- which(presence == 1)
      n_presences <- length(presence_idx)

      # --------------------------------------------------------------------------------------------------- #

      # GENERATE ABSENCE POINTS & COMBINE DATA SETS:

      # Only proceed if presences are between 1000 and 10000:
      if (n_presences >= 1000 && n_presences <= 10000) {
        # Extract presence geometries and create sf object:
        presences_geom <- cand_geom[presence_idx]
        presences_df <- st_sf(geometry = presences_geom, presence = 1)

        # Generate absence points:
        absence_points <- st_sample(gnp, size = n_presences, type = "random")
        absence_df <- st_sf(geometry = absence_points, presence = 0)

        # Combine presences and absences:
        sim_data <- rbind(presences_df, absence_df)

        data_frames[[i]] <- list(
          df = sim_data,
          equation = hyper_grid[i, ]
        )
      } else {
        data_frames[[i]] <- "Failed"
      }

      # Keep track of Loop process
      if (i %% 100 == 0) {
        cat("Creating dataframe ", i, " of ", nrow(hyper_grid), "\n")
      }
    }

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
      lai_sd
    )
  }
})


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# FILTER THE LISTS AND EDIT: ####
# ========================== #

# Keep only the data frames that have data:
data_frames <- keep(
  data_frames,
  ~ is.list(.x)
)

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

rm(coords, i, stack, vals)
# rm(hyper_grid, vals, coords, i, stack)

# Save the lists to avoid re-running the code:

{
  set.seed(3791)
  sdf2 <- data_frames[sample(length(data_frames), 200)]
  saveRDS(sdf2, "./Data/Simple_data_200_data_frames.rds")
  rm(sdf2)
}

#
{
  set.seed(3791)
  sdf3 <- data_frames[sample(length(data_frames), 50)]
  saveRDS(sdf3, "./Data/Simple_50_data_frames_for_testing.rds")
  rm(sdf3)
}


saveRDS(data_frames, file = "./Data/simple_data_frames.rds")
