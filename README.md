# SML_sampling — Sampling design for a Soil Monitoring Level (SML) network

This repository contains an R workflow that selects the **sampling sites of a soil monitoring network**
(Soil Monitoring Level, SML) in Italy. It is developed at ISPRA by Claudia Cagnarini and Giulio Barcaroli.

The included scripts are a **worked example on a small test area**: 6 domains (soil unit × land use) in
2 zones (southern Tuscany; eastern Emilia / Po delta), on a 1 km grid. They are the scaled-down version of a
national script that combines FAO/JRC raster values with a kriging-based description of spatial structure.

The workflow is split into two independent scripts:

| Script | Role | When to run it |
|---|---|---|
| `SML_sampling_example_with_6_Italy_soil_units.R` | Sampling design: frame, LUCAS points, optimal stratification and allocation, sample selection, outputs | Always |
| `Kriging_example_with_6_Italy_soil_units.R` | Spatial parameters: variogram fitting on point survey data, `legend.csv` with fitting and range per soil unit, reliability report | Only to (re)estimate the spatial parameters, before a run with `spatial_component = TRUE` |

## What the sampling script does

1. **Defines the sampling frame**: 1 km grid cells, each with model-predicted soil pH (`PH`) and soil organic
   carbon (`SOC`), their prediction variances, clay content (`ARG`) and elevation (`DEM`).
2. **Forces in existing LUCAS points** (those with `SOIL == 1` that fall inside the study area) as a
   take-all (census) stratum, and removes their cells from the frame that is actually sampled.
3. **Sets the spatial parameters** according to the flag `spatial_component`:
   - `TRUE`: reads, for every soil unit, `fitting` (R²) and `range` (km) of PH and SOC from `Input/legend.csv`;
   - `FALSE` (default): no spatial component, `fitting = 1` and `range = 0` in every domain.
4. **Optimises the stratification** with [SamplingStrata](https://barcaroli.github.io/SamplingStrata/):
   K-means initial solution (`KmeansSolutionSpatial` if `TRUE`, `KmeansSolution2` if `FALSE`), then a genetic
   algorithm (`optimStrata`) with Bethel allocation, using the soil-unit fitting/range (`TRUE`) or
   `fitting = 1`, `range = 0` (`FALSE`), so that the coefficient of variation (CV) of the estimates of PH and
   SOC is at most **5 % in every domain**.
5. **Checks the spatial component per domain** (Phase 6b, only if `TRUE`): the GA strata are re-allocated with
   and without the spatial component, and the non-spatial allocation replaces the spatial one where it needs
   fewer units.
6. **Selects the sites** (`selectSampleSpatial`) and exports points, summary tables, figures and the strata as
   polygons, ready to be used in the SML tool.

## What the kriging script does

1. Fits regression-kriging variograms (spherical model, trend `x + y + DEM`) on the point-level survey,
   nationally and per zone (a zone with fewer than 30 valid points, or a degenerate fit, uses the national model).
2. Assigns to each soil unit the fitting (R² of the variogram fit) and range of its zone and writes
   `Input/legend.csv`.
3. Assesses the reliability of the parameters of each soil unit (**ADEQUATE / USE WITH CAUTION / NOT RELIABLE**)
   from the fitting, the use of the national model, a range capped at the data support, and the size of the
   range compared with the soil unit, and writes a report with warnings and a recommendation on
   `spatial_component`.

In the example data the parameters are classified as **NOT RELIABLE** for all soil units (fitting 0.41 for PH,
0.49 for SOC, from the national model; ranges of 109 and 80 km against mean distances of 11–72 km inside the
soil units), so the sampling script is delivered with `spatial_component <- FALSE`.

Note: use the last version (1.6) of SamplingStrata by installing it as follows:

```r
pak::pak("barcaroli/SamplingStrata")
```

## Repository contents

| Path | Description |
|---|---|
| `SML_sampling_example_with_6_Italy_soil_units.R` | Sampling workflow (Phases 1–7) |
| `Kriging_example_with_6_Italy_soil_units.R` | Variogram fitting, `legend.csv` and reliability report (Phases K1–K4) |
| `Input/` | Frame (`.csv`/`.gpkg`), soil-unit legend, LUCAS points, point survey, 1 km DEM; `legend.csv` once the kriging script has been run |
| `Output/` | Results of a reference run (summary, sample, figures, strata, workspace) |
| `parameters_values_for_SML_tool.txt` | Settings to enter in the SML tool to reuse the strata produced here |
| `docs/SML_sampling_user_manual.pdf` | User manual: how to customise and run the two scripts |
| `SML_sampling.Rproj` | RStudio project |

## Requirements

R (≥ 4.1). Sampling script: `sp`, `raster`, `ggplot2`, `dplyr`, `tidyr`, `SamplingStrata` (≥ 1.6), `clhs`,
`sf`, `terra`, `readxl`, `parallel`, `doParallel`. Kriging script: `sf`, `terra`, `gstat`.
Both scripts install missing packages automatically.

## Quick start

1. Open `SML_sampling.Rproj` in RStudio.
2. Edit the path variables (`inp_path`, `out_path`) at the top of each script, which still point to a different
   folder. Create the output folder if needed.
3. Non-spatial design: run `SML_sampling_example_with_6_Italy_soil_units.R` as delivered
   (`spatial_component <- FALSE`).
4. Spatial design: run `Kriging_example_with_6_Italy_soil_units.R`, read
   `<out_path>/Example_with_6_soil_units_kriging_report.txt`, then set `spatial_component <- TRUE` in the
   sampling script and run it.
5. Read `<out_path>/Example_with_6_soil_units_summary.txt` first.

See the [user manual](docs/SML_sampling_user_manual.pdf) for the full list of parameters, inputs, outputs and
troubleshooting.

## Using the result in the SML tool

Copy `strata_SML.gpkg` to the tool's `Files/Strata/` folder and apply the settings in
`parameters_values_for_SML_tool.txt` (soil-unit field `soil_unit`, strata field `stratum_id`).
