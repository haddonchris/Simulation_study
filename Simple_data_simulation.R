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

  # Animal variables
  {
    cat("Getting Animal interaction data...", "\n")
    gnp_species <- read.csv("./Data/GNP_game_counts.csv") |>
      filter(year >= 2000) |> # Only keep counts after the year 2000
      dplyr::select(-c(X, X.1, X.2, X.3, X.4)) |> # Remove unnecessary vars
      mutate(
        species = as.factor(species),
        species = recode_factor(
          # Fix naming issues
          species,
          "Blue wildebeest" = "Blue Wildebeest",
          "oribi" = "Oribi"
        ),
        species = str_replace_all(species, "\\s", "_") # Replace spaces with underscores
      ) |>
      filter(
        species %in%
          c(
            "Blue_Wildebeest",
            "Buffalo",
            "Eland",
            "Elephant",
            "Hartebeest",
            "Impala",
            "Waterbuck",
            "Zebra"
          )
      )

    # Get the counts of each species:
    animal_counts <- gnp_species |>
      group_by(species) |>
      summarise(obs = n(), counts = sum(number))

    # Add UTM columns to the data frame
    gnp_species <- suppressMessages(add_utm_columns(
      gnp_species,
      c("long", "lat"),
      units = "m"
    )) |>
      rename(x = X, y = Y)

    # Make a spatial points data frame:
    spdf <- gnp_species
    coordinates(spdf) <- c("x", "y")
    projection(spdf) <- CRS("+init=epsg:32736")

    # Make the Utilisation Distributions
    ud_select_species <- kernelUD(spdf[, "species"], h = 'href')

    rst <- rast("./satellite_images/Elevation/elevation_export.tif")

    m <- 1000000000

    gnu <- resample(terra::rast(ud_select_species[["Blue_Wildebeest"]]), rst) *
      m
    varnames(gnu) <- "Blue_wildebeest"
    buffalo <- resample(terra::rast(ud_select_species[["Buffalo"]]), rst) * m
    varnames(buffalo) <- "Buffalo"
    eland <- resample(terra::rast(ud_select_species[["Eland"]]), rst) * m
    varnames(eland) <- "Eland"
    elephant <- resample(terra::rast(ud_select_species[["Elephant"]]), rst) * m
    varnames(elephant) <- "Elephant"
    hart <- resample(terra::rast(ud_select_species[["Hartebeest"]]), rst) * m
    varnames(hart) <- "Hartebeest"
    impala <- resample(terra::rast(ud_select_species[["Impala"]]), rst) * m
    varnames(impala) <- "Impala"
    waterbuck <- resample(terra::rast(ud_select_species[["Waterbuck"]]), rst) *
      m
    varnames(waterbuck) <- "Waterbuck"
    zebra <- resample(terra::rast(ud_select_species[["Zebra"]]), rst) * m
    varnames(zebra) <- "Zebra"

    animals <- c(gnu, buffalo, eland, elephant, hart, impala, waterbuck, zebra)
    names(animals) <- c(
      "blue_wildebeest",
      "buffalo",
      "eland",
      "elephant",
      "hartebeest",
      "impala",
      "waterbuck",
      "zebra"
    )
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
    animals <- resample(x = animals, y = rst)

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
    animals <- mask(crop(animals, gnp), gnp)

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
      water_dist,
      animals
    )
  }

  # Remove unnecessary objects
  rm(
    animals,
    buffalo,
    eland,
    elephant,
    elev,
    evi,
    gnu,
    hart,
    impala,
    lai,
    m,
    pnt,
    pnv,
    slope,
    soil,
    spdf,
    ud_select_species,
    waterbuck,
    zebra,
    animal_counts,
    gnp_species
  )
  cat("Done!!!")
}

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# GENERATE SIMULATED DATA: ####
# ======================== #

# Create a hyper grid with all of the parameters for the simulated data frames
hyper_grid <- expand.grid(
  beta_ndvi = seq(from = -2.5, to = 2.5, by = 2.5),
  beta_water = seq(from = -2.5, to = 2.5, by = 2.5),
  beta_ptc = seq(from = -2.5, to = 2.5, by = 2.5),
  alpha = seq(from = -5, to = 5, by = 2.5)
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
  cand_points$dist_water <- terra::extract(water_dist, coords)[, 1]
}

set.seed(3791)

# Simulation For Loop:
for (i in 1:nrow(hyper_grid)) {
  # CALCULATE THE UNSCALED COEFFICIENTS:

  # NDVI coefficient:
  {
    ndvi_stats <- global(ndvi, fun = c("mean", "sd"), na.rm = TRUE)
    b_ndvi <- hyper_grid[i, 1] / ndvi_stats[1, 2] # Beta for scaled NDVI = 0.5
    a_ndvi <- hyper_grid[i, 4] -
      hyper_grid[i, 1] * (ndvi_stats[1, 1] / ndvi_stats[1, 2]) # Alpha for scaled NDVI = -0.1
  }

  # Distance to water coefficient:
  {
    water_stats <- global(water_dist, fun = c("mean", "sd"), na.rm = TRUE)
    b_water <- hyper_grid[i, 2] / water_stats[1, 2]
    a_water <- hyper_grid[i, 4] -
      hyper_grid[i, 2] * (water_stats[1, 1] / water_stats[1, 2])
  }

  # Percent tree cover coefficient:
  {
    ptc_stats <- global(ptc, fun = c("mean", "sd"), na.rm = TRUE)
    b_ptc <- hyper_grid[i, 3] / ptc_stats[1, 2]
    a_ptc <- hyper_grid[i, 4] -
      hyper_grid[i, 3] * (ptc_stats[1, 1] / ptc_stats[1, 2])
  }

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

    data_frames[[i]] <- list(
      df = sim_data,
      equation = hyper_grid[i, ]
    )
  } else {
    data_frames[[i]] <- "Failed"
  }

  # --------------------------------------------------------------------------------------------------- #

  # Keep track of Loop process
  cat("Creating dataframe ", i, " of ", dim(hyper_grid)[1], "\n")
}

# Remove unnessary objects
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
  water_stats,
  ndvi,
  ptc,
  water_dist
)

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# FILTER THE LISTS AND EDIT: ####
# ========================== #

# Keep only the data frames that have data:
data_frames <- keep(
  data_frames,
  ~ is.list(.x)
)

# Add a unique ID and extract all satellite variables for each data frame:
for (i in 1:length(data_frames)) {
  # Add unique ID:
  data_frames[[i]]$mod_id <- paste0("mod_", i)

  # Extract satellite imagery data:
  coords <- terra::vect(st_coordinates(data_frames[[i]]$df)) # Get coordinates
  vals <- terra::extract(stack, coords, df = T) # Extract the values into a new data frame
  data_frames[[i]]$df <- cbind(data_frames[[i]]$df, vals[, -1]) # Combine with data

  # Keep track of the loop's progress
  cat("Extracting for data frame", i, "of", length(data_frames), "\n")
}

rm(hyper_grid, vals, coords, i, stack)

# Save the lists to avoid re-running the code:
saveRDS(data_frames, file = "./Data/simple_data_frames.rds")

big_list <- readRDS("./Data/simple_data_frames.rds")


# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# RUN MODELS FOR EACH DATA FRAME IN

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# CREATE TRUE PROBABILITY MAP: ####
# ============================ #

prob_stack <- c(ndvi, ptc, water_dist)

i <- 120

# Create a raster with based off of the equation:
logit_rst <- (data_frames[[i]]$equation$alpha +
  (data_frames[[i]]$equation$beta_ndvi * prob_stack$ndvi) +
  (data_frames[[i]]$equation$beta_water * prob_stack$water_dist) +
  (data_frames[[i]]$equation$beta_ptc * prob_stack$ptc))

# Convert to probability
p_true_rst <- 1 / (1 + exp(-logit_rst))
plot(p_true_rst, col = inferno(150))
plot(logit_rst, col = inferno(150))
