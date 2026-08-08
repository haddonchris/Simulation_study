# Simulation Study Chapter

This repository contains an R-based simulation study for Gorongosa National Park. The workflow builds environmental predictors from raster and vector data, simulates presence/absence datasets under different response shapes, fits logistic regression and random forest models, and compares predicted probabilities to the known generating surface.

## Project structure

The repository is organized into four main script groups:

- `1.x_*_data_generation.R` — generate simulated datasets
- `2.x_*_models.R` — fit logistic regression and random forest models
- `3.x_*_analysis.R` — evaluate model performance and produce diagnostic figures

## Script overview

### Data generation

### `1.1_Linear_data_generation.R`
Creates the linear simulation dataset. It loads raster predictors, computes distance to water, resamples and masks all covariates to the Gorongosa National Park boundary, calculates summary statistics for scaling, samples candidate locations, simulates presence/absence data under a linear response surface, and saves the final list of simulated data frames.

$$\Large{y_i \sim \text{Bernoulli}(p_i)}$$
$$\Large{\text{logit}(p_i) = \alpha + \beta_{1}X_{1i} + \beta_{2}X_{2i} + \beta_{3}X_{3i} + \beta_{4} X_{4i}}$$

### `1.2_Semi_non_linear_data_generation.R`
Creates the semi-non-linear simulation dataset. It follows the same general structure as the linear generator, but uses a semi-non-linear response surface and a large hyperparameter grid to produce accepted simulated data frames and calibration figures.

$$\Large{y_i \sim \text{Bernoulli}(p_i)}$$
$$\Large{\text{logit}(p_i) = \alpha + \beta_{1}X_{1i}^{z} + \beta_{2}X_{2i}^{z} + \beta_{3}X_{3i}^{z} + \beta_{4}X_{4i}^{z}}$$

### `1.3_Non_linear_data_generation.R`
Creates the non-linear simulation dataset. It uses the same environmental predictors and sampling framework, then simulates datasets under a more strongly non-linear response function and writes the final simulated data frames and calibration outputs.

$$\Large{y_i \sim \text{Bernoulli}(p_i)}$$
$$\Large{\text{logit}(p_i) = \alpha + {\beta_{1} \space e^{\left(-\frac{(x_{1i} - \mu_{x_1})^2}{2 \sigma_{x_1}^2}\right)}} + \beta_{2} (-(x_{2i}-\mu_{x_{2}})^2) + \beta_3 \space e^{-k_{x_3}x_{3i}} + \beta_4 \frac{1}{1 + e^{-k_{x_4} (x_{4i} - c)}}}$$

### Model fitting

### `2.1_Linear_models.R`
Fits logistic regression and random forest models to each linear simulated dataset. The script splits each dataset into training and testing sets, fits a logistic model with scaled predictors, tunes a random forest with a hyperparameter grid, stores the fitted models back into the simulation object, and records performance metrics such as sensitivity, specificity, TSS, and AUC.

### `2.2_Semi_non_linear_models.R`
Fits the same model types for the semi-non-linear simulation datasets. The workflow mirrors the linear model-fitting script, including data partitioning, logistic regression fitting, random forest tuning, and performance summarization.

### `2.3_Non_linear_models.R`
Fits logistic regression and random forest models for the non-linear simulation datasets. As with the other model-fitting scripts, it stores fitted models in the simulation object and calculates summary performance measures for each dataset.

### Analysis and evaluation

### `3.0_spatial_autocorrelation_test.R`
Estimates the local range of spatial autocorrelation for the environmental predictors. It samples raster values, fits quadratic trend surfaces, computes original and residual variograms, fits candidate variogram models, saves diagnostic plots, and builds a summary table of estimated ranges and variogram parameters. It also explores how hexagonal sampling distance affects the number of evaluation points.

### `3.1_Linear_analysis.R`
Evaluates the linear simulation scenario. It rebuilds the raster stack, computes true probability surfaces, predicts probabilities from the fitted logistic regression and random forest models, saves probability maps, extracts predicted values at sampling locations, and writes model comparison metrics.

### `3.2_Semi_non_linear_analysis.R`
Evaluates the semi-non-linear simulation scenario. It uses the same analysis workflow as the linear script, producing true, logistic regression, and random forest probability surfaces and exporting summary metrics and figures.

### `3.3_Non_linear_analysis.R`
Evaluates the non-linear simulation scenario. It loads the simulated datasets and fitted models, rebuilds the environmental rasters, predicts true and model-based probability surfaces, saves maps, extracts values at evaluation points, and writes summary metrics.

## Data and output folders

### `Satellite_images/`
Raster source data used to build the environmental predictor stack.

### `Gorongosa_files/`
Vector data for the park boundary, rivers, lake, and related spatial inputs.

### `Data/Final_data/`
Final RDS and CSV outputs from the simulation study, including simulated datasets, predicted probabilities, and performance summaries.

### `Data/Figures/`
Diagnostic figures, variogram plots, and predicted probability maps.

### `Data/Cal_images/`
Calibration plots created during data generation.


## Reproducibility notes

- The scripts assume the project root is the working directory.
- Most scripts begin by clearing the workspace and graphics device before loading packages.
- Several scripts rebuild the predictor stack from source rasters and shapefiles, so the raw input files must remain available.
- Some file names differ slightly across script versions, so if a script fails to find an input file, check the names in `Data/Final_data/`.

## Dependencies

Common packages used across the scripts include `terra`, `sf`, `tidyverse`, `ranger`, `caret`, `pROC`, `doSNOW`, `parallel`, `gstat`, `spdep`, `sdmTMB`, `viridis`, `sp`, `raster`, and `adehabitatHR`.
