rm(list = ls())
graphics.off()

# Load required libraries:
{
  library(terra)
  library(sf)
  library(gstat)
  library(spdep)
}

# ============================================================================================== #
# ============================================================================================== #

# LOAD RASTER DATA: ####
# ================= #

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
  lai <- mask(crop(lai, gnp), gnp)
  ptc <- mask(crop(ptc, gnp), gnp)
  water_dist <- mask(crop(water_dist, gnp), gnp)

  stack <- c(
    ndvi,
    lai,
    ptc,
    water_dist
  )

  cat("Done!!!")
}

# ============================================================================================== #
# ============================================================================================== #

# CREATE FUNCTION TO GET LOCAL RANGE OF SPATIAL AUTOCORRELATION: ####

get_local_range <- function(
  raster,
  name,
  sample_size = 5000,
  cutoff = 15000,
  width = 500,
  seed = 3791
) {
  set.seed(seed)

  # Sample raster:
  pts <- spatSample(
    raster,
    size = sample_size,
    method = "random",
    as.points = TRUE,
    values = TRUE,
    na.rm = TRUE
  )

  # Convert the sampled points to an sf object and rename the value column to the variable name
  pts <- st_as_sf(pts)
  names(pts)[1] <- name

  # Extract coordinates:
  xy <- st_coordinates(pts)
  pts$x <- xy[, 1]
  pts$y <- xy[, 2]

  # Fit quadratic trend surface:
  form <- as.formula(paste(name, "~ x + y + I(x^2) + I(y^2) + I(x*y)")) # Create a formula for the trend surface model

  trend <- lm(form, data = pts) # Fit the trend surface model using linear regression

  # Calculate the residuals from the trend surface model and add them as a new column in the pts data frame
  pts$residuals <- residuals(trend)

  # Origional variogram:
  original_form <- as.formula(
    paste(name, "~ 1")
  )

  vg_original <- variogram(
    original_form,
    data = pts,
    cutoff = cutoff,
    width = width
  )

  # Residual variogram:
  vg_residual <- variogram(
    residuals ~ 1,
    data = pts,
    cutoff = cutoff,
    width = width
  )

  # Candidate models:
  candidate_models <- c("Sph", "Exp", "Gau")
  fits <- vector("list", length(candidate_models))
  names(fits) <- candidate_models

  sse <- numeric(length(candidate_models))
  names(sse) <- candidate_models

  # Loop through the candidate models and fit each one to the residual variogram, storing the fitted model and its sum of
  # squared errors (SSE) in the fits and sse lists, respectively. Use vg_residual results to set the initial parameters for
  # the model fitting.
  for (i in seq_along(candidate_models)) {
    fits[[i]] <- fit.variogram(
      vg_residual,
      model = vgm(
        model = candidate_models[i], # Use the candidate model type
        psill = max(vg_residual$gamma), # Set the partial sill to the maximum semivariance of the residual variogram
        range = max(vg_residual$dist) / 3, # Set the range to one-third of the maximum distance of the residual variogram
        nugget = min(vg_residual$gamma) # Set the nugget to the minimum semivariance of the residual variogram
      )
    )

    # Store the sum of squared errors (SSE) for the fitted model in the sse list
    sse[i] <- attr(fits[[i]], "SSErr")
  }

  # Identify the best-fitting model by finding the model with the minimum SSE and storing its name in the best variable
  best <- names(which.min(sse))

  # Create a plot comparing the original and residual variograms, saving it as a PNG file in the specified directory. The plot includes points for both variograms, with different colors and symbols, and a legend indicating which is which.
  {
    png(
      filename = paste0(
        "./Data/Figures/Variogram_diagnostics/",
        name,
        "_comparison.png"
      ),
      width = 1400,
      height = 1000,
      res = 200
    )

    plot(
      vg_original$dist,
      vg_original$gamma,
      type = "b",
      pch = 15,
      col = "black",
      xlab = "Distance (m)",
      ylab = "Semivariance",
      main = paste(name, "- Original vs Residual")
    )

    lines(
      vg_residual$dist,
      vg_residual$gamma,
      type = "b",
      pch = 18,
      col = "red"
    )

    legend(
      "bottomright",
      legend = c("Original", "Residual"),
      pch = c(15, 18),
      col = c("black", "red"),
      bty = "n"
    )

    dev.off()
  }

  # Return a list containing the variable name, R-squared value of the trend surface model, best-fitting model type,
  # range, nugget, partial sill, SSE values for all candidate models, original variogram, fitted variogram for the
  # best model, and the trend surface model object.
  list(
    variable = name,
    r_squared = summary(trend)$r.squared,
    best_model = best,
    range = fits[[best]]$range[2],
    nugget = fits[[best]]$psill[1],
    partial_sill = fits[[best]]$psill[2],
    sse = sse,
    variogram = vg_residual,
    fit = fits[[best]],
    trend = trend
  )
} # End of get_local_range function

# ============================================================================================== #
# ============================================================================================== #

# RUN THE FUNCTION FOR EACH VARIABLE IN THE STACK: ####
# ================================================ #

# Create a list to store the results for each variable
results <- vector("list", nlyr(stack))

# Loop through each layer in the stack and apply the get_local_range function:
for (i in 1:nlyr(stack)) {
  # Get the name of the current variable
  variable_name <- names(stack)[i]

  # Print the name of the variable being processed
  cat("Processing variable:", variable_name, "\n")

  # Call the get_local_range function and store the result in the results list
  results[[i]] <- get_local_range(
    raster = stack[[i]],
    name = variable_name
  )
}

# Assign names to the results list based on the variable names in the stack
names(results) <- names(stack)

# Create a summary table to store the results for each variable
summary_table <- data.frame(
  Variable = character(), # Column for variable names
  R2 = numeric(), # Column for R-squared values of the trend surface models
  Model = character(), # Column for best-fitting model types
  Range_m = numeric(), # Column for ranges of the best-fitting models
  Nugget = numeric(), # Column for nugget values of the best-fitting models
  Partial_sill = numeric(), # Column for partial sill values of the best-fitting models
  stringsAsFactors = FALSE
)

# Loop through the results list and populate the summary table with the relevant information
for (i in seq_along(results)) {
  summary_table <- rbind(
    summary_table,
    data.frame(
      Variable = results[[i]]$variable,
      R2 = results[[i]]$r_squared,
      Model = results[[i]]$best_model,
      Range_m = results[[i]]$range,
      Nugget = results[[i]]$nugget,
      Partial_sill = results[[i]]$partial_sill
    )
  )
}

# Print the summary table to the console
summary_table

# Write the summary table to a CSV file for further analysis or reporting
write.csv(
  summary_table,
  file = "./Data/Spatial_autocorrelation_summary.csv",
  row.names = FALSE
)

# ============================================================================================== #
# ============================================================================================== #

# CREATE A PLOT TO SHOW THE EFFECT OF SAMPLING DISTANCE ON EVALUATION SAMPLE SIZE: ####
# ================================================================================ #

# Define a vector of distances (in meters) to be used for creating hexagonal grids
dist <- c(2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000)

# Create a data frame to store the number of evaluation points for each distance
n_points <- data.frame(
  n_rows = integer(length(dist)),
  dist = numeric(length(dist))
)

# Loop through each distance value to create a hexagonal grid and count the number of evaluation points
for (i in 1:length(dist)) {
  # Create a hexagonal grid
  hex <- st_make_grid(
    gnp,
    cellsize = dist[i], # 2.5 km
    square = FALSE
  )

  # Get the centroids of the hexagonal grid cells and convert them to a data frame
  pts <- st_centroid(hex)
  pts <- as.data.frame(st_coordinates(pts[gnp]))

  # Store the number of evaluation points and the corresponding distance in the n_points data frame
  n_points$n_rows[i] <- nrow(pts)
  n_points$dist[i] <- dist[i]
}

# Create a plot to visualize the relationship between hexagon spacing and the number of evaluation points
png(
  filename = "./Data/Figures/Variogram_diagnostics/Sampling_distance_vs_variogram_range.png",
  width = 1800,
  height = 1400,
  res = 300
)

{
  # Set the margins for the plot
  par(mar = c(5, 5, 4, 2))

  # Plot the number of evaluation points against hexagon spacing
  plot(
    n_points$dist,
    n_points$n_rows,
    type = "b",
    pch = 16,
    lwd = 2,
    xlab = "Hexagon spacing (m)",
    ylab = "Number of evaluation points",
    main = "Effect of sampling distance on evaluation sample size"
  )

  # Define colors for the vertical lines representing the ranges of spatial autocorrelation for each variable
  cols <- c("red", "blue", "forestgreen", "orange")

  # Loop through each variable in the summary table and add vertical lines to the plot at the corresponding range values
  for (i in 1:nrow(summary_table)) {
    abline(
      v = summary_table$Range_m[i],
      col = cols[i],
      lty = 2,
      lwd = 2
    )
  }

  # Add a legend to the plot indicating which color corresponds to which variable and its range of spatial autocorrelation
  legend(
    "topright",
    inset = c(0.12, 0),
    legend = paste0(
      summary_table$Variable,
      " (",
      round(summary_table$Range_m / 1000, 1),
      " km)"
    ),
    col = cols,
    lty = 2,
    lwd = 2,
    bty = "n",
    cex = 0.8
  )
}

# Close the graphics device
dev.off()
