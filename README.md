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
   - `TRUE`: reads `fitting` (R²) and `range` (km) of every target variable for every soil unit from the legend;
   - `FALSE`: no spatial component, `fitting = 1` and `range = 0` in every domain; the legend is not read.
3. **SamplingStrata frame**: cells with non-finite values are dropped, soil units with fewer than
   `min_point_per_domain` cells are dropped, precision constraints are set to `cv_target` in every domain.
4. **Initial solution** per domain with K-means (`KmeansSolutionSpatial` if `TRUE`, `KmeansSolution2` if
   `FALSE`).
5. **Optimal stratification and allocation** with [SamplingStrata](https://barcaroli.github.io/SamplingStrata/):
   genetic algorithm (`optimStrata`) with Bethel allocation, so that the coefficient of variation (CV) of every
   target variable is at most `cv_target` in every domain. Only if `spatial_component = TRUE`:
   - **Phase 5b**: the strata found by the GA are re-allocated with and without the spatial component; the
     non-spatial allocation is used where it needs fewer units;
   - **Phase 5c**: a second, non-spatial GA is run (`fitting = 1`, `range = 0`); for each soil unit the whole
     solution (strata and allocation) with the smaller sample size is kept. The spatial component is therefore
     used only in the soil units where it reduces the sample size.
6. **Sample selection** (`selectSampleSpatial`) and outputs: points, summary tables, report, figures and the
   strata as polygons, ready to be used in the SML tool.

## Inputs (folder `Input/`)

| File | Content |
|---|---|
| `Frame_1km_6_soil_units_Italy.csv` | Sampling frame, one row per 1 km cell: `POINT_ID`, `X_LAEA`, `Y_LAEA` (EPSG:3035), `soil_unit`, target variables with their prediction variance (`PH`, `PH_VAR`, `SOC`, `SOC_VAR`), auxiliary variables (`CLAY`, `DEM`) |
| `Test_6_soil_units_Italy_legend.csv` | Only if `spatial_component = TRUE`: one row per soil unit with `soil_unit` and, for each target variable `<var>`, `fitting_<var>` (R²) and `range_<var>_km` |
| `Test_6_soil_units_Italy.gpkg` | Soil-unit polygons: not read by the script, used by the SML tool |

The fitting and range values in the example legend are **illustrative**: they are not estimated on the example
maps. Replace them with your own estimates (see below) before using the spatial component for real results.

## How to estimate fitting (R²) and range

The two parameters describe how well the DSM model predicts each target variable and how far its errors are
spatially correlated. They are best derived from the same model used to produce the maps, separately for each
soil unit and each target variable, using the soil observations (e.g. LUCAS or regional surveys) that fall in
the soil unit.

- **fitting**: R² between observed and predicted values. Out-of-sample values (cross-validation or an
  independent validation set) are more realistic than the R² of the calibration. Allowed values are in (0, 1].
- **range**: range of the variogram of the model residuals, in km. SamplingStrata models the correlation between
  two cells at distance *d* as exp(−*d* / range), so fit an **exponential** model and enter its range
  parameter, not the practical range (which is about 3 times larger).

Example with `gstat`, for one soil unit and one variable (`obs` is an `sf` object of the observations with the
observed value `PH`, the predicted value `PH_pred` and coordinates in metres):

```r
library(gstat)
obs$res <- obs$PH - obs$PH_pred
fitting <- cor(obs$PH, obs$PH_pred)^2              # better on validation data
v       <- variogram(res ~ 1, obs)
vfit    <- fit.variogram(v, vgm("Exp"))
range_km <- vfit$range[vfit$model == "Exp"] / 1000  # range parameter, in km
```

Some practical advice:

- A variogram needs enough observations with a good spatial spread: with only a few dozen points per soil unit
  the estimate is unstable. In that case, estimate the parameters on a group of similar soil units (same soil
  region or land use) and assign them to each unit of the group.
- Check the range against the size of the soil unit: a range much larger than the soil unit, or much smaller
  than the 1 km grid, usually means that the variogram is poorly identified.
- When no reliable estimate is possible, set `fitting = 1` and `range = 0` for that soil unit: it is then
  treated as non-spatial.
- Low fitting values increase the expected variance, and therefore the sample size: Phases 5b and 5c make sure
  that the spatial component is kept only where it actually reduces the sample.

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
| `spatial_component` | `TRUE` to use fitting and range from the legend (Phases 5b and 5c are performed) |
| `use_lucas_grid` | `TRUE` to flag the selected cells that belong to the LUCAS 2 km grid |
| `base_name` | Prefix of the output files |

## Outputs (folder `Output/`, prefix `base_name`)

| File | Content |
|---|---|
| `_summary.txt` | Report: spatial parameters, optimisation, check of the spatial component (Phases 5b and 5c), final sample |
| `_units.csv` | One row per soil unit: cells, whether it was stratified, selected points, of which on the LUCAS grid |
| `_allpoints.csv` | One row per selected point: soil unit, stratum, cell id, coordinates (EPSG:3035), LUCAS grid flag |
| `_sample.csv` | SamplingStrata sample, with weights |
| `_figures.pdf` | Correlation matrix, CV boxplot, sample maps, points per soil unit |
| `.gpkg` | 1 km grid with soil unit and sample counts |
| `_strata_SML.gpkg` | Strata as polygons, for the SML tool |
| `.RData` | Full workspace |

## Requirements

R (≥ 4.1) with `ggplot2`, `dplyr`, `tidyr`, `sf` and `SamplingStrata` (≥ 1.6). Missing packages are installed
automatically; SamplingStrata 1.6 is installed from GitHub:

```r
pak::pak("barcaroli/SamplingStrata")
```

## Quick start

1. Open `SML_sampling.Rproj` in RStudio.
2. Edit `inp_path` and `out_path` at the top of the script and create the output folder if needed.
3. Set the parameters, in particular `spatial_component` (`TRUE` requires the legend with fitting and range).
4. Run `SML_optimal_sampling_tool.R`.
5. Read `<out_path>/<base_name>_summary.txt` first.

## Using the result in the SML tool

Copy `<base_name>_strata_SML.gpkg` to the tool's `Files/Strata/` folder and apply the settings in
`parameters_values_for_SML_tool.txt` (soil-unit field `soil_unit`, strata field `stratum_id`).
