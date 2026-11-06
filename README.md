# SML_sampling — Sampling design for a Soil Monitoring Level (SML) network

This repository contains an R workflow that selects the **sampling sites of a soil monitoring network**
(Soil Monitoring Level, SML) in Italy. It is developed at ISPRA by Claudia Cagnarini and Giulio Barcaroli.

The included script is a **worked example on a small test area**: 6 domains (soil unit × land use) in
2 zones (southern Tuscany; eastern Emilia / Po delta), on a 1 km grid. It is the scaled-down version of a
national script that combines FAO/JRC raster values with a kriging-based description of spatial structure.

## What the workflow does

1. **Defines the sampling frame**: 1 km grid cells, each with model-predicted soil pH (`PH`) and soil organic
   carbon (`SOC`), their prediction variances, clay content (`ARG`) and elevation (`DEM`).
2. **Forces in existing LUCAS points** (those with `SOIL == 1` that fall inside the study area) as a
   take-all (census) stratum, and removes their cells from the frame that is actually sampled.
3. **Fits variograms** (regression kriging, spherical model, trend `x + y + DEM`) on point-level survey data,
   nationally and per zone, to quantify spatial autocorrelation (fitting R² and range).
4. **Optimises the stratification** with [SamplingStrata](https://barcaroli.github.io/SamplingStrata/):
   spatial K-means initial solution, then a genetic algorithm with Bethel allocation, so that the
   coefficient of variation (CV) of the estimates of PH and SOC is at most **5 % in every domain**.
5. **Decides per domain** whether taking spatial correlation into account (Phase 6b) reduces the sample size
   compared with the non-spatial baseline.
6. **Selects the sites** (`selectSampleSpatial`) and exports points, summary tables, figures and the strata as
   polygons, ready to be used in the SML tool.
   
Note: use the last version (1.6) of SamplingStrata by installing it as follows:

pak::pak("barcaroli/SamplingStrata")


## Repository contents

| Path | Description |
|---|---|
| `SML_sampling_example_with_6_Italy_soil_units.R` | The full workflow (7 phases, ~1,300 lines) |
| `Input/` | Frame (`.csv`/`.gpkg`), soil-unit legend, LUCAS points, point survey, 1 km DEM |
| `output/` | Results of a reference run (summary, sample, figures, strata, workspace) |
| `parameters_values_for_SML_tool.txt` | Settings to enter in the SML tool to reuse the strata produced here |
| `docs/SML_sampling_user_manual.pdf` | User manual: how to customise and run the script |
| `SML_sampling.Rproj` | RStudio project |

## Requirements

R (≥ 4.1) with the packages `sp`, `raster`, `ggplot2`, `dplyr`, `tidyr`, `SamplingStrata`, `bethel`, `clhs`,
`sf`, `terra`, `readxl`, `gstat`, `parallel`, `doParallel`.

## Quick start

1. Open `SML_sampling.Rproj` in RStudio.
2. Edit the three path variables (`inp_path`, `out_path`, `kriging_survey_path`) at the top of the script,
   which still point to a different folder.
3. Create the output folder if needed, then run the script (`source()` or step by step).
4. Read `<out_path>/Example_with_6_soil_units_summary.txt` first.

See the [user manual](docs/SML_sampling_user_manual.pdf) for the full list of parameters, inputs, outputs and
troubleshooting.

## Using the result in the SML tool

Copy `strata_SML.gpkg` to the tool's `Files/Strata/` folder and apply the settings in
`parameters_values_for_SML_tool.txt` (soil-unit field `soil_unit`, strata field `stratum_id`).
