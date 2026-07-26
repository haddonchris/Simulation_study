rm(list = ls())
graphics.off()

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

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

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# READ IN THE DATA: ####
# ================= #

sdf <- readRDS("./Data/Final_data/Simple_data_200_data_frames.rds")

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# EXTRACT MODEL PERFORMANCE METRICS FOR LOGISTIC REGRESSION: ####
# ========================================================== #

n <- length(sdf) # Number of data frames/models to evaluate

system.time({
  {
    # Set seed once for reproducible, but iteration-varying, data splits
    set.seed(3791)

    # Train the models and calculate the performance metrics:
    for (i in 1:n) {
      cat("Running logistic regression", i, "of", length(sdf), "\n") # Print the progress of the model fitting and evaluation

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
    node_size = c(1, 2, 3, 4, 6, 8, 10, 15),
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
      cat("Running random forest", i, "of", length(sdf), "\n")

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
    } # End of loop

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

saveRDS(
  sdf,
  "./Data/Final_data/Simple_200_data_frames_with_models.rds"
)
