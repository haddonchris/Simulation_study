rm(list = ls())
graphics.off()

# LOAD LIBRARIES: ####
# =============== #

{
  library(tidyverse)
  library(terra)
  library(sf)
  library(viridis)
  library(ranger)
}

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# LOAD IN THE GEOGRAPHIC & SATELLITE DATA: ####
# ======================================== #

{
  # NDVI:
  {
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
    rm(ndvi_files)
    # plot(ndvi_mean, col = plasma(150))
  }

  # Percent Canopy Cover:
  {
    ptc_files <- data.frame(
      files = list.files(
        # List all of the file names and add as a column in a data frame
        path = "./Satellite_images/percent_canopy_cover",
        pattern = "\\.tif$",
        full.names = TRUE
      )
    )

    # Create a Raster Stac:
    ptc_stack <- mean(terra::rast(ptc_files$files)) # Read in the files as a raster stack
    rm(ptc_files)
  }

  # Percent non-tree vegetation:
  {
    pnt_files <- data.frame(
      files = list.files(
        # List all of the file names and add as a column in a data frame
        path = "./Satellite_images/percent_non_tree",
        pattern = "\\.tif$",
        full.names = TRUE
      )
    ) |>
      mutate(
        date = as.POSIXct(as.factor(substr(files, start = 37, stop = 46))) # Add a column just for the month
      )

    # Create a Raster Stack:
    pnt_stack <- mean(terra::rast(pnt_files$files)) # Read in the files as a raster stack
    rm(pnt_files)
  }

  # Distance to water:
  {
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
    # plot(water_dist, col = inferno(150))
    rm(rivers, lake, rst, lake_rst, rivers_rst, water)
  }

  # Resample rasters to the same resolution:
  rst <- terra::rast("./Satellite_images/Elevation/elevation_export.tif")

  ndvi <- resample(x = ndvi_mean, y = rst)
  ptc <- resample(x = ptc_stack, y = rst)
  pnt <- resample(x = pnt_stack, y = rst)
  water <- resample(x = water_dist, y = rst)

  rm(rst, ndvi_mean, ptc_stack, pnt_stack, water_dist)

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
  ndvi <- mask(crop(ndvi, gnp), gnp)
  ptc <- mask(crop(ptc, gnp), gnp)
  pnt <- mask(crop(pnt, gnp), gnp)
  water <- mask(crop(water, gnp), gnp)

  stack <- c(ndvi, ptc, pnt, water)
  names(stack) <- c("ndvi", "ptc", "pnt", "water_dist")
  varnames(stack) <- c("ndvi", "ptc", "pnt", "water_dist")
}

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# GENERATE SIMULATED DATA: ####
# ======================== #

# Create a hyper grid with all of the parameters for the simulated data frames
hyper_grid <- expand.grid(
  beta_ndvi = seq(from = -2.5, to = 2.5, by = 0.5),
  beta_water = seq(from = -2.5, to = 2.5, by = 0.5),
  beta_ptc = seq(from = -2.5, to = 2.5, by = 0.5),
  alpha = seq(from = -5, to = 5, by = 0.5)
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
  # cand_points$pnt <- terra::extract(pnt, coords)[, 1]
  cand_points$dist_water <- terra::extract(water, coords)[, 1]
}

# Simulation For Loop:

for (i in 1:nrow(hyper_grid)) {
  # CALCULATE THE UNSCALED COEFFICIENTS:

  # NDVI coefficient:
  ndvi_stats <- global(ndvi, fun = c("mean", "sd"), na.rm = TRUE)
  b_ndvi <- hyper_grid[i, 1] / ndvi_stats[1, 2] # Beta for scaled NDVI = 0.5
  a_ndvi <- hyper_grid[i, 4] -
    hyper_grid[i, 1] * (ndvi_stats[1, 1] / ndvi_stats[1, 2]) # Alpha for scaled NDVI = -0.1

  # Distance to water coefficient:
  water_stats <- global(water, fun = c("mean", "sd"), na.rm = TRUE)
  b_water <- hyper_grid[i, 2] / water_stats[1, 2]
  a_water <- hyper_grid[i, 4] -
    hyper_grid[i, 2] * (water_stats[1, 1] / water_stats[1, 2])

  # Percent tree cover coefficient:
  ptc_stats <- global(ptc, fun = c("mean", "sd"), na.rm = TRUE)
  b_ptc <- hyper_grid[i, 3] / ptc_stats[1, 2]
  a_ptc <- hyper_grid[i, 4] -
    hyper_grid[i, 3] * (ptc_stats[1, 1] / ptc_stats[1, 2])

  # --------------------------------------------------------------------------------------------------- #

  # DEFINE THE EQUATION, CALCULATE PROBABILITIES, AND SIMULATE PRESENCE:

  # The Equation
  logit_p <- (hyper_grid[i, 4] + b_ndvi * cand_points$NDVI) +
    (b_water * cand_points$dist_water) +
    (b_ptc * cand_points$ptc)

  # Convert logit to probability:
  cand_points$p_true <- plogis(logit_p)

  # Simulate presences based on probability:
  presences <- cand_points |>
    filter(p_true > 0.5)

  # --------------------------------------------------------------------------------------------------- #

  # GENERATE ABSENCE POINTS & COMBINE DATA SETS:

  # Create a `if else` statement so that if presences are > 1000 & <10 000 its defined as "Failed"
  if (nrow(presences) >= 1000 && nrow(presences) <= 10000) {
    # Generate absence points:
    {
      set.seed(3791)
      absence_points <- st_sample(
        gnp,
        size = nrow(presences),
        type = "random"
      )
    }

    absence_points <- st_sf(geometry = absence_points)

    # Add new variable, presence
    absence_points$presence <- 0
    presences <- presences |>
      mutate(presence = 1) |>
      dplyr::select(c("geometry", "presence")) # keep only presence and geometry to add the two data frames

    # Combine the absences and presences
    sim_data <- rbind(presences, absence_points)

    # Extract the co-ordinates to extract the other variables
    coords <- st_coordinates(sim_data)

    # Extract covariate values for simulated data:
    sim_data$NDVI <- terra::extract(ndvi, coords)[, 1]
    sim_data$ptc <- terra::extract(ptc, coords)[, 1]
    sim_data$pnt <- terra::extract(pnt, coords)[, 1]
    sim_data$dist_water <- terra::extract(water, coords)[, 1]
    sim_data$presence <- as.factor(sim_data$presence)

    data_frames[[i]] <- list(
      df = sim_data,
      equation = hyper_grid[i, ]
    )
  } else {
    data_frames[[i]] <- "Failed"
  }

  cat("Creating dataframe ", i, " of ", dim(hyper_grid)[1], "\n")
}

rm(
  a_ndvi,
  a_ptc,
  a_water,
  b_ndvi,
  b_ptc,
  b_water,
  i,
  logit_p,
  absence_points,
  cand_points,
  coords,
  gnp,
  ndvi_stats,
  presences,
  ptc_stats,
  sim_data,
  water_stats
)

# Save the lists to avoid re-running the code:
saveRDS(data_frames, file = "simulated_data_frames.rds")
saveRDS(test, file = "filtered_simulated_data_frames.rds")

big_list <- readRDS("simulated_data_frames.rds")
smaller_list <- readRDS("filtered_simulated_data_frames.rds")


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

test <- keep(
  big_list,
  ~ is.list(.x)
)

# Give each model a unique ID
for (i in 1:2455) {
  test[[i]]$ID <- paste0("mod_", i)
}

plot()


plot(test[[6]])

test[[6]]$df

Sdf_test <- test[[6]]

plot(df_test)
plot(df_test[, geometry])

df_test <- df_test |>
  filter(presence == 1)

plot(df_test)
plot()

lol <- cat("Poes")
lol
lol <- "Poes"
