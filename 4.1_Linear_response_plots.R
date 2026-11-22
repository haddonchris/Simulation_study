rm(list = ls())
graphics.off()

{
  library(terra)
  library(sf)
  library(viridis)
  library(ranger)
  library(tidyverse)
}

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

sdf <- readRDS("./Data/Final_data/Linear_data_with_models.rds")

# ======================================================= #
# CREATE RESPONSE PLOTS FOR EACH VARIABLE FOR EACH MODEL:
# ======================================================= #

system.time({
  # response_plots <- setNames(
  #   vector("list", length(levels(sable_data$sable_id))),
  #   levels(sable_data$sable_id)
  # )

  vars <- c("ndvi", "ptc", "lai", "water_dist")

  for (i in 1:10) {
    cat("Creating partial plots for", sdf[[i]]$mod_id, "\n")

    data <- as.data.frame(sdf[[i]]$df) |>
      dplyr::select(-"geometry")

    # response_plots[[i]] <- list()

    for (j in seq_along(vars)) {
      focal_var <- vars[j]

      base_vals <- data[vars][rep(1, 200), ]

      cat("- Variable:", focal_var, "\n")

      for (k in seq_along(vars)) {
        base_vals[[vars[k]]] <- mean(data[[vars[k]]], na.rm = TRUE)
      }

      base_vals[[focal_var]] <- seq(
        quantile(data[[focal_var]], 0.05, na.rm = TRUE),
        quantile(data[[focal_var]], 0.95, na.rm = TRUE),
        length.out = 200
      )

      {
        {
          b_ndvi <- sdf[[i]]$equation$beta_ndvi / ndvi_stats$sd
          b_water <- sdf[[i]]$equation$beta_water / water_stats$sd
          b_ptc <- sdf[[i]]$equation$beta_ptc / ptc_stats$sd
          b_lai <- sdf[[i]]$equation$beta_lai / lai_stats$sd # Using LAI stats for scaling LAI
        }

        true_logit <- (
          # Intercept
          sdf[[i]]$equation$alpha +

            # NDVI response
            (b_ndvi * (base_vals$ndvi)) +

            # Distance-to-water response
            (b_water * (base_vals$water_dist)) +

            # PTC response
            (b_ptc * (base_vals$ptc)) +

            # LAI response
            (b_lai * (base_vals$lai))
        )

        base_vals$True <- plogis(true_logit)
      }

      x_lab <- ifelse(
        focal_var == "ndvi",
        "NDVI",
        ifelse(
          focal_var == "ptc",
          "Canopy Cover (%)",
          ifelse(
            focal_var == "lai",
            "Leaf Area Index",
            ifelse(
              focal_var == "water_dist",
              "Distance to Water (m)",
              focal_var
            )
          )
        )
      )

      # response_plot <- data.frame(
      #   Predictor = base_vals[[focal_var]],
      #   True = base_vals$True,
      #   RF = predict(sdf[[i]]$rf_mod, data = base_vals)$predictions[, "1"],
      #   LR = predict(sdf[[i]]$lr_mod, newdata = base_vals, type = "response")
      # ) |>
      #   pivot_longer(
      #     cols = c(True, RF, LR),
      #     names_to = "Model",
      #     values_to = "Probability"
      #   ) |>
      #   ggplot(
      #     aes(
      #       x = Predictor,
      #       y = Probability,
      #       colour = Model
      #     )
      #   ) +
      #   geom_line(linewidth = 2) +
      #   labs(
      #     y = "Predicted Probability",
      #     x = x_lab
      #   ) +
      #   scale_colour_manual(
      #     values = c(
      #       "True" = "black",
      #       "LR" = "#A87B38",
      #       "RF" = "#535E31"
      #     ),
      #     labels = c(
      #       "True" = "True relationship",
      #       "LR" = "Logistic Regression",
      #       "RF" = "Random Forest"
      #     )
      #   ) +
      #   theme(
      #     panel.background = element_rect(
      #       fill = "white",
      #       colour = "white"
      #     ),
      #     axis.line = element_line(
      #       colour = "black",
      #       linewidth = 0.8
      #     ),
      #     axis.text = element_text(
      #       size = 16,
      #       colour = "black"
      #     ),
      #     axis.title = element_text(
      #       size = 20
      #     ),
      #     legend.key.height = unit(0.65, "cm"),
      #     legend.key.width = unit(0.65, "cm"),
      #     legend.key = element_rect(
      #       fill = "white",
      #       colour = "white"
      #     ),
      #     legend.title = element_blank(),
      #     legend.text = element_text(
      #       size = 18
      #     ),
      #     legend.box = element_blank(),
      #     legend.position = "bottom"
      #   )

      # response_plot

      test <- data.frame(
        Predictor = base_vals[[focal_var]],
        True_logit = true_logit
      ) |>
        ggplot(aes(x = Predictor, y = True_logit)) +
        geom_line()

      ggsave(
        paste0(
          "./Data/Figures/Response_test/Response_plot_",
          sdf[[i]]$mod_id,
          "_",
          vars[j],
          "_test.png"
        ),
        plot = test,
        units = "cm",
        width = 25,
        height = 20,
        dpi = 1500
      )

      # ggsave(
      #   paste0(
      #     "./Data/Figures/Response_test/Response_plot_",
      #     sdf[[i]]$mod_id,
      #     "_",
      #     vars[j],
      #     ".png"
      #   ),
      #   plot = response_plot,
      #   units = "cm",
      #   width = 25,
      #   height = 20,
      #   dpi = 1500
      # )

      # response_plots[[i]][[focal_var]] <- response_plot
    }
  }

  rm(
    i,
    j,
    k,
    vars,
    focal_var,
    x_lab,
    base_vals,
    response_plot,
    data,
    b_lai,
    b_ndvi,
    b_ptc,
    b_water,
    true_logit
  )
})
