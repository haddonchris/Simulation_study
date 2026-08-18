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

sdf <- readRDS("./Data/Final_data/Non_linear_data.rds")

# --------------------------------------------------------------------------------------------------- #
# --------------------------------------------------------------------------------------------------- #

# EXTRACT MODEL PERFORMANCE METRICS FOR LOGISTIC REGRESSION: ####
# ========================================================== #

# Create an empty data frame to store the model performance metrics:
model_performance <- data.frame(
  mod_id = character(), # Model type (e.g., logistic regression, random forest)
  sim = character(),
  mod_type = character(), # Unique identifier for each model/data frame
  Sensitivity = numeric(), # Sensitivity (true positive rate)
  Specificity = numeric(), # Specificity (true negative rate)
  TSS = numeric(), # True Skill Statistic (TSS = Sensitivity + Specificity - 1)
  AUC = numeric(), # Area Under the Curve (AUC) for the ROC curve
  stringsAsFactors = FALSE # Prevent automatic conversion of character columns to factors
)

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
          scale(lai) +
          scale(ptc) +
          scale(water_dist),
        data = train,
        family = binomial(link = "logit")
      )

      # Store the fitted model in the list for later use:
      sdf[[i]]$lr_mod <- mod

      # Make predictions and create a confusion matrix:
      pred_lr <- predict(mod, newdata = test, type = "response")
      pred_class_lr <- ifelse(pred_lr >= 0.5, 1, 0) # Convert probabilities to class labels based on a threshold of 0.5
      pred_class_lr <- factor(pred_class_lr, levels = c(1, 0)) # Convert to factor with levels 1 and 0
      obs <- factor(test$presence, levels = c(1, 0)) # Convert observed values to factor with levels 1 and 0
      conf_mat <- caret::confusionMatrix(pred_class_lr, obs)

      # Store the model performance metrics in the data frame:
      model_performance[i, 1] <- sdf[[i]]$mod_id
      model_performance[i, 2] <- "NL"
      model_performance[i, 3] <- "Logistic Regression"
      model_performance[i, 4] <- conf_mat$byClass["Sensitivity"]
      model_performance[i, 5] <- conf_mat$byClass["Specificity"]
      model_performance[i, 6] <- (conf_mat$byClass["Sensitivity"] +
        conf_mat$byClass["Specificity"] -
        1)
      model_performance[i, 7] <- auc(roc(test$presence, pred_lr, quiet = TRUE))
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
    mtry = c(1, 2, 3, 4),
    node_size = c(1, 2, 3, 4, 6, 8, 10, 15)
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
            num.trees = 750,
            mtry = hyper_grid$mtry[j],
            min.node.size = hyper_grid$node_size[j],
            sample.fraction = 0.7,
            seed = 3791,
            num.threads = 1, # Ensure that each model runs on a single thread to prevent nested parallelism issues
            probability = TRUE
          )

          # Store the OOB RMSE and prediction error in the hyperparameter grid:
          data.frame(
            Brier = model$prediction.error,
            OOB_RMSE = sqrt(model$prediction.error)
          )
        }

      close(pb) # Close the progress bar after the loop is complete

      # ------------------------------------------------------------------------------------------------ #

      results <- cbind(hyper_grid, results) # Combine the hyperparameter grid with the results

      # Arrange the hyperparameter grid by prediction error
      best <- results |>
        dplyr::slice_min(Brier, n = 1)

      # Fit the final random forest model with the best hyperparameters:
      rf <- ranger(
        formula = factor(presence, levels = c(0, 1)) ~ .,
        data = train,
        num.trees = 750,
        mtry = best$mtry,
        min.node.size = best$node_size,
        sample.fraction = 0.7,
        seed = 3791,
        write.forest = T,
        importance = "impurity",
        probability = T
      )

      sdf[[i]]$rf_mod <- rf

      # Make predictions and create a confusion matrix:
      pred_rf <- predict(rf, data = test) # Get the predicted probabilities for the positive class (presence = 1)
      pred_rf <- pred_rf$predictions[, "1"] # Extract the predicted probabilities for the positive class (presence = 1)
      preds_rf <- ifelse(pred_rf >= 0.5, 1, 0) # Convert probabilities to class labels based on a threshold of 0.5
      pred_class_rf <- factor(preds_rf, levels = c(1, 0)) # Convert to factor with levels 1 and 0
      obs <- factor(test$presence, levels = c(1, 0)) # Convert observed values to factor with levels 1 and 0

      # Create the confusion matrix using the caret package:
      conf_mat <- confusionMatrix(pred_class_rf, obs)

      # Store the model performance metrics in the data frame:
      model_performance[i + n, 1] <- sdf[[i]]$mod_id
      model_performance[i + n, 2] <- "NL"
      model_performance[i + n, 3] <- "Random Forest"
      model_performance[i + n, 4] <- conf_mat$byClass["Sensitivity"]
      model_performance[i + n, 5] <- conf_mat$byClass["Specificity"]
      model_performance[i + n, 6] <- conf_mat$byClass["Sensitivity"] +
        conf_mat$byClass["Specificity"] -
        1
      model_performance[i + n, 7] <- auc(roc(
        test$presence,
        pred_rf,
        quiet = TRUE
      ))
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
  "./Data/Final_data/Non_linear_data_with_models.rds"
)

write.csv(
  model_performance,
  "./Data/Final_data/Non_linear_model_performance.csv",
  row.names = F
)
