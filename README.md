# Simulation Study Chapter

This repository contains an R-based simulation study for Gorongosa National Park. It builds environmental raster predictors, simulates presence/absence datasets under different response shapes, fits statistical and machine-learning models, and compares predicted probabilities to known truth.

## Project workflow

The project is organized into three main stages:

1. **Data generation** — scripts in the `1.x_*_data_generation.R` series create simulated datasets for linear, semi-non-linear, and non-linear response surfaces.
2. **Model fitting** — scripts in the `2.x_*_models.R` series fit the models used in the simulation study.
3. **Analysis and evaluation** — scripts in the `3.x_*_analysis.R` series compare predicted and true probabilities, calculate error metrics, and save diagnostic figures.

A separate set of `Testing_*.R` scripts contains exploratory and diagnostic work.

## Main scripts tracked in this repository

### `1.2_Semi_non_linear_data_generation.R`
Generates semi-non-linear simulated presence/absence datasets. It loads raster predictors from `Satellite_images/`, computes distance to water from vector data in `Gorongosa_files/`, resamples and masks all layers to the Gorongosa National Park boundary, builds a large hyperparameter grid, simulates candidate datasets, and saves the first valid set of simulations along with calibration PNGs.

### `3.1_Linear_analysis.R`
Evaluates the linear simulation scenario. It rebuilds the environmental raster stack, creates sampling points, computes true probabilities, predicts with logistic regression and random forest models, saves probability maps, and calculates summary error metrics.

### `3.2_Semi_non_linear_analysis.R`
Evaluates the semi-non-linear simulation scenario. It follows the same general workflow as the linear analysis script: rebuild the raster stack, generate sampling points, compute true probabilities, predict with fitted models, save figures, and summarize prediction error.

### `3.3_Non_linear_analysis.R`
Evaluates the non-linear simulation scenario. It loads the simulated datasets and fitted models, rebuilds the raster predictors, predicts true, logistic regression, and random forest probability surfaces, writes prediction maps, and exports probability and metric summaries.

### `Testing_11.R`
A diagnostic and exploratory script for the non-linear scenario. It computes model comparison metrics and produces quick visual summaries for checking distributions of error statistics.

## Key folders

- `Data/Final_data/` — final RDS and CSV outputs from the simulation study
- `Data/Figures/` — diagnostic figures, variogram plots, and predicted probability maps
- `Data/Cal_images/` — calibration plots created during data generation
- `Satellite_images/` — raster source data used to build predictors
- `Gorongosa_files/` — spatial vector data for the park boundary, rivers, and lake

## Notes

- The scripts assume the project root is the working directory.
- Some script and output file names do not match exactly across versions of the project. If a script fails to find an input file, check the names in `Data/Final_data/`.
- Several scripts create large raster objects, so memory usage can be high.

## Reproducibility

Most scripts begin by clearing the workspace and graphics device, loading the required packages, and rebuilding the predictor stack from source files. This means they are intended to be run from a clean session.
