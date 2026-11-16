# SML optimal sampling tool — Sampling design for soil monitoring under the Soil Monitoring Law

This repository contains an R workflow that selects the **sampling sites of a soil monitoring network**
according to the Soil Monitoring Law (SML). It is developed at ISPRA by Claudia Cagnarini and Giulio Barcaroli.

The included data are a **worked example on a small test area in Italy**: 6 soil units (soil unit × land use,
e.g. `BA_Cropland`, `60.4_Forest`) on a 1 km grid. The soil units are the estimation domains.

## What the script does

The script `SML_optimal_sampling_tool.R` runs in six phases:

1. **Set-up and domains**: reads the parameters and the sampling frame (1 km cells with the DSM predictions of
   the target variables, their prediction variance and the auxiliary variables). Cells of the LUCAS 2 km grid
   (`POINT_ID > 0`) are flagged.
2. **Spatial parameters**, according to the flag `spatial_component`:
   - `TRUE`: reads `range` (km) and `nugget` of every target variable for every soil unit from the legend;
   - `FALSE`: no spatial component, `range = 0` in every domain; the legend is not read.
3. **SamplingStrata frame**: cells with non-finite values are dropped, soil units with fewer than
   `min_point_per_domain` cells are dropped, precision constraints are set to `cv_target` in every domain.
4. **Initial solution** per domain with K-means (`KmeansSolutionSpatial` if `TRUE`, `KmeansSolution2` if
   `FALSE`).
5. **Optimal stratification and allocation** with [SamplingStrata](https://barcaroli.github.io/SamplingStrata/):
   genetic algorithm (`optimStrata`) with Bethel allocation, so that the coefficient of variation (CV) of every
   target variable is at most `cv_target` in every domain. Only if `spatial_component = TRUE`, **Phase 5b**
   (report only) re-allocates the strata found by the GA without the spatial component, to show the reduction
   of the sample size due to it. The GA result is the final allocation.
6. **Sample selection** (`selectSampleSpatial`) and outputs: points, summary tables, report, figures and the
   strata as polygons, ready to be used in the SML tool.

## Inputs (folder `Input/`)

| File | Content |
|---|---|
| `Frame_1km_6_soil_units_Italy.csv` | Sampling frame, one row per 1 km cell: `POINT_ID`, `X_LAEA`, `Y_LAEA` (EPSG:3035), `soil_unit`, target variables with their prediction variance (`PH`, `PH_VAR`, `SOC`, `SOC_VAR`), auxiliary variables (`CLAY`, `DEM`) |
| `Test_6_soil_units_Italy_legend.csv` | Only if `spatial_component = TRUE`: one row per soil unit with `soil_unit` and, for each target variable `<var>`, `range_<var>_km` (range parameter of the exponential variogram of the prediction errors, km) and `nugget_<var>` (relative nugget c0/(c0+c1); 0 if the column is missing) |
| `Test_6_soil_units_Italy.gpkg` | Soil-unit polygons: not read by the script, used by the SML tool |

The range and nugget values in the example legend are **illustrative**: they are not estimated on the example
maps. Replace them with your own estimates (see below) before using the spatial component for real results.

## How to estimate range and nugget

The variance in the strata is computed from the DSM predictions and their prediction-error variances (additive
correction, de Gruijter et al. 2016): for two cells i, j of a stratum

    D2_ij = (y_i - y_j)^2 + s2_i + s2_j - 2 (1 - nugget) s_i s_j exp(-d_ij / range)

The prediction-error variances already account for the levelling (smoothing) of the predictions, so no R² of the
model is used: dividing by R² as well would count the error variance twice. The two spatial parameters describe
the correlation of the prediction errors, and are estimated, separately for each soil unit and each target
variable, on the variogram of the **out-of-sample residuals** (observed minus predicted values, from
cross-validation or an independent validation set) of the soil observations (e.g. LUCAS or regional surveys) that
fall in the soil unit. At the calibration points the residuals of the model are close to zero and not
informative about the prediction error.

- **range**: range parameter *a* of an **exponential** variogram, in km. The correlation at distance *d* is
  exp(−*d* / *a*); do not enter the practical range (about 3*a*), because the script uses `kappa = 1`.
- **nugget**: relative nugget c0/(c0+c1), in [0, 1]: the share of the error variance (micro-scale variability and
  measurement error) that is not spatially correlated, even between nearby cells.

Example with `gstat`, for one soil unit and one variable (`obs` is an `sf` object of the observations with the
observed value `PH`, the out-of-sample prediction `PH_pred`, the prediction-error variance `PH_var` and coordinates
in metres):

```r
library(gstat)
obs$res <- (obs$PH - obs$PH_pred) / sqrt(obs$PH_var)   # standardised residuals
v       <- variogram(res ~ 1, obs)
vfit    <- fit.variogram(v, vgm(NA, "Exp", NA, NA))     # exponential model with nugget
range_km <- vfit$range[vfit$model == "Exp"] / 1000      # range parameter, in km
nugget   <- vfit$psill[vfit$model == "Nug"] / sum(vfit$psill)
```

Standardised residuals are used because the error variance varies from cell to cell, while the formula assumes a
common correlation structure.

Some practical advice:

- A variogram needs enough observations with a good spatial spread: with only a few dozen points per soil unit
  the estimate is unstable. In that case, estimate the parameters on a group of similar soil units (same soil
  region or land use) and assign them to each unit of the group.
- Check the range against the size of the soil unit: a range much larger than the soil unit, or much smaller
  than the 1 km grid, usually means that the variogram is poorly identified.
- When no reliable estimate is possible, set `range = 0` (or `nugget = 1`) for that soil unit: it is then
  treated as non-spatial.
- With this formula the spatial component can only reduce the variance in the strata; a large nugget makes the
  reduction correspondingly smaller.
- Old legends with `fitting_<var>` columns: values < 1 are rejected; remove those columns.

## Main parameters (top of the script)

| Parameter | Meaning |
|---|---|
| `inp_path`, `out_path` | Input and output folders |
| `vars_cont_y`, `variance_cont_y` | Target variables and their prediction variance |
| `vars_cont_x` | Auxiliary (stratification) variables |
| `min_point_per_domain` | Soil units with fewer cells are dropped |
| `cv_target` | Target CV of every target variable in every domain (default 0.05) |
| `min_sample_per_stratum` | Minimum number of sampled units per stratum (≥ 2) |
| `max_clusters` | Maximum number of K-means clusters per domain |
| `ga_iter`, `ga_pops` | Generations and population size of the GA: the defaults are for quick tests; use much higher values for final results |
| `spatial_component` | `TRUE` to use range and nugget from the legend (Phase 5b reports the reduction due to the spatial component) |
| `use_lucas_grid` | `TRUE` to flag the selected cells that belong to the LUCAS 2 km grid |
| `base_name` | Prefix of the output files |

## Outputs (folder `Output/`, prefix `base_name`)

| File | Content |
|---|---|
| `_summary.txt` | Report: spatial parameters, optimisation, reduction due to the spatial component (Phase 5b), final sample |
| `_units.csv` | One row per soil unit: cells, whether it was stratified, selected points, of which on the LUCAS grid |
| `_allpoints.csv` | One row per selected point: soil unit, stratum, cell id, coordinates (EPSG:3035), LUCAS grid flag |
| `_sample.csv` | SamplingStrata sample, with weights |
| `_figures.pdf` | Correlation matrix, CV boxplot, sample maps, points per soil unit |
| `.gpkg` | 1 km grid with soil unit and sample counts |
| `_strata_SML.gpkg` | Strata as polygons, for the SML tool |
| `.RData` | Full workspace |

## Requirements

R (≥ 4.1) with `ggplot2`, `dplyr`, `tidyr`, `sf` and `SamplingStrata` in a version that has the `nugget`
parameter of `optimStrata` (later than 1.6). Missing packages are installed automatically; SamplingStrata is
installed from GitHub:

```r
pak::pak("barcaroli/SamplingStrata")
```

## Quick start

1. Open `SML_sampling.Rproj` in RStudio.
2. Edit `inp_path` and `out_path` at the top of the script and create the output folder if needed.
3. Set the parameters, in particular `spatial_component` (`TRUE` requires the legend with range and nugget).
4. Run `SML_optimal_sampling_tool.R`.
5. Read `<out_path>/<base_name>_summary.txt` first.

## Using the result in the SML tool

Copy `<base_name>_strata_SML.gpkg` to the tool's `Files/Strata/` folder and apply the settings in
`parameters_values_for_SML_tool.txt` (soil-unit field `soil_unit`, strata field `stratum_id`).
