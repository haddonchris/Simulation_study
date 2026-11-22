# Machine learning for habitat selection: a simulation study

**MSc thesis, Chapter 2 · University of the Witwatersrand**

**Author:** Chris Haddon Garrett

This repository holds the R code for a simulation study that compares **logistic regression (LR)** and **random forest (RF)** models for habitat selection. Presence/absence data are simulated across Gorongosa National Park, Mozambique, from real environmental predictors under three response surfaces that are progressively more non-linear. Because the generating ("true") probability surface is known, each model's predictions can be checked directly against the truth.

> **Hypothesis:** as the species–environment relationship becomes more non-linear, RF's predictive performance should improve relative to LR.

---

## Contents

- [Study design](#study-design)
- [Workflow](#workflow)
- [Simulation scenarios](#simulation-scenarios)
- [Script reference](#script-reference)
- [Repository structure](#repository-structure)
- [Getting started](#getting-started)
- [Input data](#input-data)
- [Reproducibility notes](#reproducibility-notes)
- [Dependencies](#dependencies)
- [Citation and contact](#citation-and-contact)
- [License](#license)

---

## Study design

| Component | Details |
|---|---|
| Study area | Gorongosa National Park, Mozambique |
| Predictors | NDVI, LAI, distance to water, percent tree cover (PTC) |
| Scenarios | Linear, semi-non-linear, non-linear |
| Replicates | 200 simulated datasets per scenario |
| Candidate locations | Pool of 20,000 locations, shared across scenarios (seed `3791`) |
| Models | Logistic regression (scaled predictors); random forest (`ranger`, tuned over a hyperparameter grid) |
| Primary metrics (vs. truth) | RMSE, bias, and correlation between predicted and true probabilities |
| Secondary metrics | AUC, TSS, sensitivity, specificity |
| Spatial checks | Variogram-based estimates of predictor autocorrelation range, used to set the spacing of evaluation points |
| Model comparison | Pairwise (per-dataset) comparison of RF and LR performance; bootstrapped overlap between metric distributions |

The primary metrics are the point of the study: they compare predictions with the known generating surface, which is only possible with simulated data. AUC and TSS are reported as well so that results can be compared with the thesis's empirical case-study chapter, where the true surface is unknown.

---

## Workflow

```mermaid
flowchart LR
    A["Raw inputs<br/>Satellite_images/<br/>Gorongosa_files/"] --> B["1.x Data generation<br/>simulate presence/absence"]
    B --> C["2.x Model fitting<br/>LR and RF per dataset"]
    C --> D["3.x Analysis<br/>predicted vs. true surfaces"]
    A --> E["3.0 Spatial autocorrelation<br/>variograms, evaluation spacing"]
    E --> D
    C --> G["4.x Response plots<br/>partial response curves"]
    D --> H["5 Cross-scenario analysis<br/>overlap, pairwise comparison"]
    D --> F["Data/Final_data/<br/>Data/Figures/"]
    G --> F
    H --> F
```

Each scenario (linear `x.1`, semi-non-linear `x.2`, non-linear `x.3`) runs through the same stages: **generate → fit → evaluate → plot responses**. `5_Analysis.R` then combines the three scenarios.

---

## Simulation scenarios

In all three scenarios the observed response at location $i$ is drawn as

$$\Large{y_i \sim \text{Bernoulli}(p_i)}$$

The scenarios differ only in how $\text{logit}(p_i)$ depends on the four predictors.

Predictors: **NDVI** = normalised difference vegetation index, **LAI** = leaf area index, **PTC** = percent tree cover, **DTW** = distance to water. $\alpha$ is the intercept and $\beta_1$–$\beta_4$ are the coefficients.

### 1. Linear

$$\Large{\text{logit}(p_i) = \alpha + \beta_1\text{NDVI}_i + \beta_2\text{LAI}_i + \beta_3\text{PTC}_i + \beta_4\text{DTW}_i}$$

Predictors are scaled.

### 2. Semi-non-linear

$$\Large{\text{logit}(p_i) = \alpha + \beta_1\text{NDVI}_i^{z_1} + \beta_2\text{LAI}_i^{z_2} + \beta_3\text{PTC}_i^{z_3} + \beta_4\text{DTW}_i^{z_4}}$$

$z_1$–$z_4$ are power-transformation exponents, one per predictor. They bend each relationship away from linearity while keeping it monotonic.

### 3. Non-linear

$$\Large{\text{logit}(p_i) = \alpha + \beta_1 e^{-\frac{(\text{NDVI}_i - \mu_1)^2}{2\sigma_1^2}} + \beta_2\left(-\left(\frac{\text{LAI}_i - \mu_2}{r_2}\right)^2\right) + \beta_3\frac{1}{1 + e^{-k_1(\text{PTC}_i - c_1)}} + \beta_4 e^{-k_2\text{DTW}_i}}$$

Each predictor has a different functional form:

| Predictor | Response shape | Parameters |
|---|---|---|
| NDVI | Gaussian (unimodal peak) | $\mu_1$ = optimum, $\sigma_1$ = width |
| LAI | Inverted quadratic (unimodal peak) | $\mu_2$ = optimum, $r_2$ = width |
| PTC | Sigmoid (logistic threshold) | $k_1$ = steepness, $c_1$ = inflection point |
| DTW | Exponential decay | $k_2$ = decay rate |

The optima $\mu_1$ and $\mu_2$ are constrained to fall inside the predictor ranges, which keeps the Gaussian and quadratic terms non-monotonic.

---

## Script reference

Run the scripts in numerical order within each scenario.

### 1.x Data generation

| Script | Description |
|---|---|
| `1.1_Linear_data_generation.R` | Loads raster predictors, computes distance to water, then resamples and masks all covariates to the park boundary. Calculates summary statistics for scaling, samples candidate locations, simulates presence/absence under the linear surface, and saves the list of simulated data frames. |
| `1.2_Semi_non_linear_data_generation.R` | Same structure as `1.1`, using the semi-non-linear surface. Draws parameter combinations from a large grid, keeps the datasets that pass the acceptance criteria, and writes calibration figures. |
| `1.3_Non_linear_data_generation.R` | Same predictors and sampling framework, using the non-linear surface. Writes the accepted simulated data frames and calibration outputs. |

### 2.x Model fitting

| Script | Description |
|---|---|
| `2.1_Linear_models.R` | Splits each dataset into training and testing sets, fits LR with scaled predictors, and tunes RF over a hyperparameter grid. Stores the fitted models in the simulation object and records sensitivity, specificity, TSS and AUC. |
| `2.2_Semi_non_linear_models.R` | Same workflow as `2.1`, for the semi-non-linear datasets. |
| `2.3_Non_linear_models.R` | Same workflow as `2.1`, for the non-linear datasets. |

### 3.x Analysis and evaluation

| Script | Description |
|---|---|
| `3.0_spatial_autocorrelation_test.R` | Estimates the range of spatial autocorrelation in the predictors. Samples raster values, fits quadratic trend surfaces, and computes original and residual variograms. Fits candidate variogram models, saves diagnostic plots, and tabulates estimated ranges and parameters. Also tests how hexagonal sampling distance affects the number of evaluation points. |
| `3.1_Linear_analysis.R` | Rebuilds the raster stack, computes the true probability surface, and predicts surfaces from the fitted LR and RF models. Saves probability maps, extracts values at evaluation points, and writes comparison metrics. |
| `3.2_Semi_non_linear_analysis.R` | Same workflow as `3.1`, for the semi-non-linear scenario. |
| `3.3_Non_linear_analysis.R` | Same workflow as `3.1`, for the non-linear scenario. |

### 4.x Response plots

| Script | Description |
|---|---|
| `4.1_Linear_response_plots.R` | Builds partial response curves for each predictor from the fitted models: the focal predictor varies across its 5th–95th percentile while the others are held at their means. Runs on a subset of the simulated datasets and saves the plots to `Data/Figures/Response_test/`. |
| `4.2_Semi_non_linear_response_plots.R` | Same workflow as `4.1`, for the semi-non-linear scenario. |
| `4.3_Non_linear_response_plots.R` | Same workflow as `4.1`, for the non-linear scenario. |

### 5 Cross-scenario analysis

| Script | Description |
|---|---|
| `5_Analysis.R` | Combines the truth-comparison metrics from all three scenarios. Compares RF and LR performance pairwise for each dataset and estimates bootstrapped overlap between the metric distributions. Saves density plots and boxplots to `Data/Figures/`. |

---

## Repository structure

```text
.
├── 1.1_Linear_data_generation.R
├── 1.2_Semi_non_linear_data_generation.R
├── 1.3_Non_linear_data_generation.R
├── 2.1_Linear_models.R
├── 2.2_Semi_non_linear_models.R
├── 2.3_Non_linear_models.R
├── 3.0_spatial_autocorrelation_test.R
├── 3.1_Linear_analysis.R
├── 3.2_Semi_non_linear_analysis.R
├── 3.3_Non_linear_analysis.R
├── 4.1_Linear_response_plots.R
├── 4.2_Semi_non_linear_response_plots.R
├── 4.3_Non_linear_response_plots.R
├── 5_Analysis.R
├── LICENSE
├── README.md
├── Satellite_images/        # Raster sources for the predictor stack (not tracked)
├── Gorongosa_files/         # Vector data: park boundary, rivers, lake (not tracked)
└── Data/                    # All outputs (not tracked)
    ├── Final_data/          # RDS/CSV outputs: simulated data, models, metrics
    ├── Figures/             # Diagnostics, variograms, maps, response plots, boxplots
    └── Cal_images/          # Calibration plots from data generation
```

The input and output folders are excluded from Git because of their size. See [Input data](#input-data).

---

## Getting started

1. **Clone the repository**
   ```bash
   git clone https://github.com/haddonchris/Simulation_study.git
   ```
2. **Open the project root as your working directory.** All paths in the scripts are relative to it.
3. **Install the dependencies** (see [Dependencies](#dependencies)).
4. **Add the raw inputs** to `Satellite_images/` and `Gorongosa_files/` (see [Input data](#input-data)), and create the `Data/` output folders.
5. **Run the pipeline in order:**
   ```text
   3.0  →  1.1 → 2.1 → 3.1 → 4.1    (linear)
           1.2 → 2.2 → 3.2 → 4.2    (semi-non-linear)
           1.3 → 2.3 → 3.3 → 4.3    (non-linear)
                             ↓
                      5_Analysis.R  (all scenarios)
   ```
   `3.0` informs the spacing of evaluation points used by the `3.x` analysis scripts. The three scenario pipelines are independent of each other. `4.x` needs only the fitted models from `2.x`. `5_Analysis.R` needs the `3.x` outputs from all three scenarios.

---

## Input data

The raw inputs are not tracked in this repository because of their size. The scripts expect them in the following locations:

| Location | Contents | Used for |
|---|---|---|
| `Satellite_images/NDVI/` | NDVI rasters (`.tif`), averaged into a mean NDVI layer | Predictor |
| `Satellite_images/LAI/` | Leaf area index rasters (`.tif`), averaged into a mean LAI layer | Predictor |
| `Satellite_images/percent_canopy_cover/` | Percent tree cover rasters (`.tif`), averaged into a mean layer | Predictor |
| `Satellite_images/Elevation/elevation_export.tif` | Elevation raster | Template grid for resampling all predictors to a common resolution |
| `Gorongosa_files/Rivers_pans/` | Permanent rivers and Lake Urema shapefiles | Distance-to-water predictor |
| `Gorongosa_files/GNP_shp/` | Gorongosa National Park boundary (World Database on Protected Areas, September 2025 release) | Cropping and masking |

The input files are available from the author on request.

---

## Reproducibility notes


- **Working directory:** the scripts assume the project root is the working directory.
- **Clean sessions:** most scripts clear the workspace and graphics device before loading packages, so each one runs standalone.
- **Random seeds:** all three data-generation scripts use seed `3791` for the candidate location pool, so every scenario draws on the same 20,000 locations and covariate values. Parameter draws and presence/absence outcomes still differ between scenarios, because each equation uses the random number stream differently and the acceptance filters differ. Datasets are therefore **not** paired across scenarios and should be treated as independent groups in any cross-scenario comparison. Within a scenario, LR and RF are fitted to the same 200 datasets, so model results are paired per dataset.
- **Raw inputs:** several scripts rebuild the predictor stack from the source rasters and shapefiles, so the raw input files must stay available.
- **File names:** some file names differ slightly between script versions. If a script can't find an input file, check the names in `Data/Final_data/`.
- **Parallel processing:** model fitting and tuning use `doSNOW` / `parallel`. Run times depend on the number of available cores.

---

## Dependencies

Developed in R. The main packages are:

| Purpose | Packages |
|---|---|
| Spatial data | `terra`, `sf`, `sp`, `raster`, `adehabitatHR` |
| Data handling | `tidyverse` |
| Modelling | `ranger`, `caret`, `sdmTMB` |
| Evaluation | `pROC` |
| Spatial statistics | `gstat`, `spdep` |
| Cross-scenario analysis | `effectsize`, `overlapping`, `DescTools`, `car` |
| Parallel computing | `doSNOW`, `parallel` (base R) |
| Visualisation | `viridis`, `cowplot`, `ggh4x` |

Install them with:

```r
install.packages(c(
  "terra", "sf", "sp", "raster", "adehabitatHR",
  "tidyverse", "ranger", "caret", "sdmTMB", "pROC",
  "gstat", "spdep", "doSNOW", "viridis",
  "effectsize", "overlapping", "DescTools", "car",
  "cowplot", "ggh4x"
))
```

---

## Citation and contact

This code accompanies an MSc thesis in Animal, Plant and Environmental Sciences at the University of the Witwatersrand, supervised by Francesca Parrini, Jason Marshal and Kiara Haylock.

If you use or adapt this code, please cite the thesis:

> Garrett, C. H. (in prep.). *Machine learning for habitat selection: a simulation study.* MSc thesis, University of the Witwatersrand, Johannesburg.

For questions, open an issue on this repository or contact the author through GitHub ([@haddonchris](https://github.com/haddonchris)).

---

## License

The code in this repository is released under the [MIT License](LICENSE).