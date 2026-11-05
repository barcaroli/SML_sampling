# =============================================================================
# Script for optimal soil monitoring site selection in Italy (SML)
# EXAMPLE: 6 soil units in 2 zones (Toscana meridionale; Emilia orientale -
# delta del Po), on a 1 km grid
#
# Authors : Claudia Cagnarini, Giulio Barcaroli
# Adapted : 21/10/2026, from the national hybrid script (FAO/JRC raster
#           values + kriging-based spatial structure, 15/09/2026)
# =============================================================================
#
# PURPOSE
# -------
# Select the sampling sites of a soil monitoring network for a small test
# area: 6 domains, each one a soil unit x land use (LULUCF) combination,
# grouped in 2 zones:
#   - Toscana meridionale          : 60.4_Cropland, 60.4_Forest, 60.7_Cropland
#   - Emilia orientale - delta Po  : BE_Cropland, BA_Cropland, DO_Cropland
# The sample is optimised with SamplingStrata (optimal stratification by
# genetic algorithm + Bethel allocation) so that the target variables, soil
# pH (PH) and soil organic carbon (SOC), reach a 5% CV in every domain. The
# soil units are the estimation domains; the auxiliary variables used to
# build the strata are clay content (ARG) and elevation (DEM).
#
#
# INPUT (folder Input/)
# ---------------------
#   frame_SML_1km_6_soil_units.csv   1 km cells of the 6 soil units (frame)
#   Test_6_soil_units_Italy.gpkg     legend: one polygon per soil unit, with
#                                    its zone (field "zona")
#   Soil_observed_points_2022_Italy.xlsx (sheet "Points")
#                                    LUCAS points (X_WGS84, Y_WGS84, SOIL)
#   SharewithMS20260825-IT.csv       point-level soil survey (pH_H2O, SOC)
#                                    used only for the variogram fitting
#   DEM_3035_1000m.tif               elevation raster, 1 km, EPSG:3035; used
#                                    only for the DEM trend of the variogram
#                                    fit (Phase 3). It is always read, but it
#                                    affects the result only if
#                                    use_spatial_fit = TRUE. The DEM used to
#                                    build the strata comes from the frame.
#
# PROCESSING PHASES
# -----------------------------------------------------------------------------
#
# Phase 1  - Set-up and domains. Parameters, frame and legend reading, construction of the
#            grid (grid_df: coordinates x/y, cell id, dom_code = zona*soil_unit,
#            LUCAS_GRID flag from ORIGINE), correlation matrix of PH, SOC, ARG
#            and DEM (PDF).
# Phase 2  - LUCAS points to be retained. Points with SOIL == 1 are projected
#            to EPSG:3035 and linked to their grid cell through POINT_ID;
#            those outside the 6 soil units, or with missing values, are
#            dropped. The retained ones (frame_cens) are forced into the
#            sample as a census stratum, and their cells are removed from the
#            frame that is actually sampled (grid_samp).
# Phase 3  - Variogram fitting on the point survey. For PH and SOC a
#            regression-kriging variogram (trend: x + y + DEM, spherical
#            model) is fitted on the raw survey points: nationally, and
#            separately for each zone (zone assigned to the survey points by
#            spatial join with the legend; zones with fewer than
#            min_obs_district valid points, or with a degenerate fit, fall
#            back to the national model). Output: variogram PDF, tables
#            fitting_by_nuts (R2) and range_by_nuts_km, one row per zone.
#            NOTE: only a small part of the survey (about 110 of 1400 valid points for pH)
#            lies inside the study area, and its empirical variograms are flat
#            (Emilia) or based on too few points (Toscana), so the fits are
#            not reliable at zone level.
#            With use_spatial_fit = FALSE the tables are overridden with
#            fitting = 1 and range = 0 in every zone (non-spatial); with TRUE
#            (default here) they are used as fitted.
# Phase 4  - SamplingStrata frame. Cells with non-finite values are dropped;
#            the frame is built with buildFrameSpatial (X = ARG, DEM;
#            Y = PH, SOC; variance = PH_VAR, SOC_VAR; domain = soil unit);
#            domains with fewer than min_point_per_domain cells are dropped
#            and domain values are re-indexed; each domain gets the
#            fitting/range of its zone; precision constraints are set to
#            CV = 5% on both PH and SOC in every domain.
# Phase 5  - Initial solution by spatial K-means, domain by domain
#            (KmeansSolutionSpatial, at most 5 clusters per domain), with a
#            single-stratum fallback if it fails. It uses the zone fitting/
#            range, so with use_spatial_fit = TRUE the starting number of
#            strata can differ from the FALSE case.
# Phase 6  - Genetic algorithm (optimStrata, method = "spatial") with
#            fitting = 1 and range = 0 everywhere, i.e. the non-spatial
#            baseline, starting from the K-means number of strata. The
#            forced LUCAS points enter as a census stratum.
# Phase 6b - Per-domain choice between non-spatial and spatial settings. If
#            use_spatial_fit = TRUE, the strata of each domain are
#            re-allocated with the zone's real fitting/range and the spatial
#            version is kept only where it needs no more units; elsewhere the
#            GA result is kept as is (not the recomputed one, to avoid
#            solver noise). If use_spatial_fit = FALSE the GA result is kept
#            in every domain and no comparison is made.
# Phase 7  - Sample selection (selectSampleSpatial) and output assembly:
#            sampled points, forced LUCAS points, summary table per domain.
#
# OUTPUT (folder Output/, prefix Example_with_6_soil_units)
# -------------------------------------------------------
#   _correlation_matrix.pdf   correlation matrix of PH, SOC, ARG, DEM
#   _variograms.pdf           empirical variograms and fits (national, zones)
#   _summary.txt              text report: variogram fits (national, zones),
#                             non-spatial optimisation (Phase 6), spatial
#                             comparison (Phase 6b), warnings (e.g. low
#                             fitting), final sample
#   _units.csv                summary per domain: cells, sampled, LUCAS
#   _allpoints.csv            all selected points (sampled + LUCAS)
#   _sample.csv               SamplingStrata sample (with weights)
#   _figures.pdf              CV boxplot, sample maps, points per domain
#   .gpkg                     1 km grid with domain and sample counts
#   .RData                    full workspace
#   strata_SML.gpkg           strata as polygons (copy it to the Files/Strata/
#                             folder of the SML tool to use it there)
#
# =============================================================================
library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(SamplingStrata)
library(bethel)
library(clhs)
library(sf)
library(terra)
library(readxl)
library(gstat)
library(parallel)
library(doParallel)


# Modify the following paths:
inp_path <- "D:\\Google Drive\\ISPRA\\SML_sampling_procedure_example\\Input\\"
out_path <- "D:\\Google Drive\\ISPRA\\SML_sampling_procedure_example\\Output\\"
kriging_survey_path <- "D:/Google Drive/ISPRA/SML_sampling_procedure_example\\Input\\SharewithMS20260825-IT.csv"


frame <- read.csv(paste0(inp_path,"frame_SML_1km_6_soil_units.csv"))

legend <- st_read(paste0(inp_path,"Test_6_soil_units_Italy.gpkg"))

# Grid of 1 km cells (the sampling frame): one row per cell,
# with the coordinates/cell id/domain/LUCAS-grid flag the rest of the script
# expects. Domain = soil unit, as "zona*soil_unit": the first component (zona)
# is the key into fitting_by_nuts / range_by_nuts_km further below.
grid_df <- frame
grid_df$x          <- grid_df$X_LAEA
grid_df$y          <- grid_df$Y_LAEA
grid_df$cell       <- seq_len(nrow(grid_df))
grid_df$dom_code   <- paste(grid_df$zona, grid_df$soil_unit, sep = "*")
grid_df$LUCAS_GRID <- as.integer(grid_df$ORIGINE == "LUCAS_2km")
#----------------------------------------------------------------------------------------------------
# Parameters & settings (to be modified by users)
#----------------------------------------------------------------------------------------------------
min_point_per_domain <- 5 # domains (soil units) with fewer cells than this are dropped

# TRUE  = use the kriging-based fitting/range (Phase 3) in the K-means step and
#         in the Phase 6b per-domain comparison. The GA itself always runs
#         non-spatial (fitting = 1, range = 0). In each domain the spatial
#         version is kept only if it needs no more units than the non-spatial
#         one; otherwise the GA result is kept untouched. The fitting enters
#         the spatial variance at the denominator, so a low fitting (e.g.
#         ~0.4-0.5 of the national model) inflates the variance and the
#         comparison is expected to select the non-spatial version.
# FALSE = non-spatial everywhere (fitting = 1, range = 0), no comparison.
#         Use it when the local variograms are not reliable (flat empirical
#         variograms in Emilia; too few points in Toscana) and the DEM raster
#         is then not needed for the result.
use_spatial_fit <- TRUE

# Fitting (R2) below this value triggers a warning in the summary report
# (_summary.txt): the fitting enters the spatial variance at the denominator,
# so a low fitting is expected to inflate the variance with the spatial component.
fitting_warn_threshold <- 0.9

vars_cont_y <- c("PH", "SOC")
variance_cont_y <- c("PH_VAR", "SOC_VAR")
vars_cont_x <- c("ARG", "DEM")

use_lucas_points <- TRUE # TRUE if some LUCAS points are imposed (forced into the sample)
use_lucas_grid <- TRUE # TRUE if lucas grid is taken into account in the frame
base_names <- if (use_lucas_grid) c("LUCAS_GRID", "zona", "soil_unit") else c("zona", "soil_unit")

# Point-level soil survey data, used ONLY for the Phase 3 kriging-based
# variogram fitting (spatial-autocorrelation source). vars_cont_y (frame
# names) are mapped to the survey's own column names.

kriging_source_col  <- c(PH = "pH_H2O", SOC = "SOC")   # vars_cont_y -> survey column

base_nome <- "Example_with_6_soil_units"

#--------------------------------------------------------------------------------------------------
# Correlation matrix of continuous variables in the grid
#--------------------------------------------------------------------------------------------------
cor_vars <- c("PH", "SOC", "ARG","DEM")
cor_mat  <- cor(frame[, cor_vars], use = "pairwise.complete.obs")
message("Correlation matrix (PH, SOC, ARG, DEM):")
print(round(cor_mat, 3))

# Correlation matrix plot (colour-coded with coefficients)
plot_cor_matrix <- function(m) {
  n   <- ncol(m)
  pal <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(201)
  op  <- par(mar = c(1, 5, 5, 1))
  on.exit(par(op))
  image(1:n, 1:n, t(m[n:1, ]), col = pal, zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "",
        main = "Correlation matrix - grid")
  axis(3, at = 1:n, labels = colnames(m), tick = FALSE, las = 1)
  axis(2, at = 1:n, labels = rev(rownames(m)), tick = FALSE, las = 1)
  for (i in 1:n) for (j in 1:n) {
    text(j, n - i + 1, sprintf("%.2f", m[i, j]),
         col = ifelse(abs(m[i, j]) > 0.6, "white", "black"))
  }
}
plot_cor_matrix(cor_mat)

# Save the correlation matrix plot to PDF
pdf(paste0(out_path,paste0(base_nome, "_correlation_matrix.pdf")), width = 6, height = 6)
plot_cor_matrix(cor_mat)
dev.off()


#--------------------------------------------------------------------------------------------------
# Read Lucas points to be retained
#--------------------------------------------------------------------------------------------------
if (use_lucas_points) {
  lucas_raw <- read_excel(paste0(inp_path,"Soil_observed_points_2022_Italy.xlsx"),
    sheet = "Points"
  )
  plot(lucas_raw$X_WGS84, lucas_raw$Y_WGS84)

  lucas_soil <- subset(lucas_raw, SOIL == 1)

  lucas_v <- project(
    vect(lucas_soil, geom = c("X_WGS84", "Y_WGS84"), crs = "EPSG:4326"),
    "EPSG:3035"
  )

  # LUCAS points sit exactly on grid nodes: link each to its grid cell via
  # POINT_ID (NA = point outside the study area, dropped by !is.na(dom_code)).
  lucas_grid_row <- match(lucas_soil$POINT_ID, grid_df$POINT_ID)

  lucas_pts <- data.frame(
    lucas_soil[, setdiff(names(lucas_soil), names(grid_df))],
    crds(lucas_v) |> as.data.frame() |> setNames(c("x", "y")),
    grid_df[lucas_grid_row, setdiff(names(grid_df), c("x", "y"))]
  )

  lucas_pts_fixed_val <- subset(lucas_pts,
#                                (GRASSLAND == 1 | LF == 1) &
!is.na(dom_code))

  ok <- complete.cases(lucas_pts_fixed_val[, c(vars_cont_y, variance_cont_y, vars_cont_x)])
  lucas_points      <- lucas_pts_fixed_val[ ok, ]
  lucas_points_filt <- lucas_pts_fixed_val[!ok, ]

  message(sprintf(
    "Lucas points summary:\n  - Outside valid domains (no dom_code): %d\n  - Excluded (missing continuous vars):  %d\n  - Total retained:                      %d\n  - Grand total (SOIL == 1):             %d",
    sum(is.na(lucas_pts$dom_code)),
    nrow(lucas_points_filt),
    nrow(lucas_points),
    nrow(lucas_pts)
  ))

  lucas_points$id_lucas     <- seq_len(nrow(lucas_points))
  lucas_points$dom_code_num <- as.integer(factor(lucas_points$dom_code,
                                                 levels = levels(factor(grid_df$dom_code))))

  # Check the auxiliary (X) variables for negative values. SamplingStrata
  # expects non-negative auxiliary variables, so any variable containing
  # negative values is shifted upward (v - min(v)) to make all values >= 0.
  neg_x_cens <- vars_cont_x[sapply(vars_cont_x, function(v) any(lucas_points[[v]] < 0, na.rm = TRUE))]
  if (length(neg_x_cens) > 0) {
    for (v in neg_x_cens) {
      shift <- min(lucas_points[[v]], na.rm = TRUE)
      warning(sprintf(
        "frame_cens: variable '%s' contains negative values (min = %.4f): shifting by %.4f to make all values positive.",
        v, shift, -shift))
      lucas_points[[v]] <- lucas_points[[v]] - shift
    }
    message("frame_cens: shifted auxiliary variable(s) with negative values: ",
            paste(neg_x_cens, collapse = ", "))
  } else {
    message("frame_cens: no negative values found in auxiliary variables (vars_cont_x).")
  }

  frame_cens <- buildFrameSpatial(
    df          = lucas_points,
    id          = "id_lucas",
    X           = vars_cont_x,
    Y           = vars_cont_y,
    variance    = variance_cont_y,
    lon         = "x",
    lat         = "y",
    domainvalue = "dom_code_num"
  )



  grid_samp <- subset(grid_df, !cell %in% unique(lucas_points$cell))
}


#---------------------------------------------------------------------------
# Frame creation
#---------------------------------------------------------------------------
# Convert the grid to an sf spatial object (not used for variogram fitting:
# see Phase 3 below, which fits on the survey point data).
dat <- st_as_sf(grid_df, coords = c("X_LAEA", "Y_LAEA"), crs = 3035, remove = FALSE)

# ---------------------------------------------------------------------------
# Phase 3 -- Regression-kriging variogram fitting on POINT survey data, per
# zone.
#
# The variogram is fitted on the RAW point-level soil survey observations
# (SharewithMS20260825-IT.csv), not on the already-smoothed predicted values
# of the frame, which would show artificially strong spatial autocorrelation.
#
# For each target variable (PH, SOC) a trend v ~ x + y + DEM (regression
# kriging / kriging with external drift) is removed and a spherical variogram
# is fitted on the residuals, first nationally (all survey points) and then
# separately for each zone of the legend ("Toscana meridionale", "Emilia
# orientale - delta del Po"). The zone of each survey point comes from a
# spatial join with the legend polygons. A zone with fewer than
# min_obs_district valid observations, or whose fit is degenerate (R2 <= 0),
# falls back to the national model. Only the variogram FIT
# (fitting_by_nuts / range_by_nuts_km, one row per zone) is used here: no
# kriging PREDICTION is performed, since Y and variance for the sampling
# frame are already in the frame (see Phase 4 below).
#
# CAVEAT: only about 110 of the ~1400 valid pH survey points fall inside the 6 soil
# units. In Emilia the empirical variograms are flat (no spatial structure
# between 2 and 76 km; the spherical fit does not converge and is rejected),
# in Toscana there are only 12 (PH) / 10 (SOC) valid points. Both zones
# therefore fall back to the national model, which is not confirmed by local
# data (national R2 ~0.41 for PH, ~0.49 for SOC, range ~109 / ~80 km). Since
# the fitting enters the spatial variance at the denominator, such low values
# are not expected to reduce the sample size: Phase 6b decides per domain,
# and with use_spatial_fit = FALSE the tables are overridden below with
# fitting = 1, range = 0 (non-spatial).
# ---------------------------------------------------------------------------
dem_rast <- rast(paste0(inp_path, "DEM_3035_1000m.tif"))

# Per-variable trend covariate set: plain DEM for both (trend: x + y + DEM).
trend_covars <- list(
  PH = c("DEM"),
  SOC = c("DEM")
)

soil_survey <- read.csv(kriging_survey_path)
survey_pts  <- vect(soil_survey, geom = c("X_LAEA", "Y_LAEA"), crs = "EPSG:3035")
survey_covs <- terra::extract(dem_rast, survey_pts, ID = FALSE)
names(survey_covs) <- c("DEM")
soil_survey <- cbind(soil_survey, survey_covs)
dat_survey  <- st_as_sf(soil_survey, coords = c("X_LAEA", "Y_LAEA"), crs = st_crs(legend), remove = FALSE)
# assign each survey point the "zona" of the legend polygon it falls in
dat_survey  <- st_join(dat_survey, legend[, "zona"], left = TRUE)

min_obs_district <- 30  # below this, a district falls back to the national model

fit_variogram_model <- function(v_col, sub, covars) {
  # Reduce() over per-column !is.na() masks, rather than complete.cases() on a
  # matrix subset: `sub` is an sf object, and indexing it with `sub[, covars]`
  # keeps the geometry list-column, which complete.cases() cannot handle.
  ok <- !is.na(sub[[v_col]])
  for (cc in covars) ok <- ok & !is.na(sub[[cc]])
  sub <- sub[ok, ]
  if (nrow(sub) < 15) return(NULL)
  form  <- as.formula(paste(v_col, "~ X_LAEA + Y_LAEA +", paste(covars, collapse = " + ")))
  emp_v <- tryCatch(variogram(form, data = sub), error = function(e) NULL)
  if (is.null(emp_v) || nrow(emp_v) < 3) return(NULL)
  m0 <- vgm(
    psill  = var(sub[[v_col]], na.rm = TRUE) * 0.9,
    model  = "Sph",
    range  = 300000,
    nugget = var(sub[[v_col]], na.rm = TRUE) * 0.1
  )
  fit <- tryCatch(fit.variogram(emp_v, model = m0), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  pred <- variogramLine(fit, dist_vector = emp_v$dist)
  SSE  <- sum((emp_v$gamma - pred$gamma)^2)
  SST  <- sum((emp_v$gamma - mean(emp_v$gamma))^2)
  R2   <- 1 - SSE / SST
  range_km <- fit$range[2] / 1000

  # Reject a degenerate fit rather than returning it: a non-positive R2 means
  # the spherical model fits the empirical semivariance no better than a flat
  # line. The caller falls back to the national model when this returns NULL.
  if (R2 <= 0) return(NULL)

  # A fitted range far beyond the empirical variogram's own data support (>
  # range_frac_max of its max lag) is an extrapolation artefact; cap it at
  # the data support (the conservative "at least this far" reading).
  range_frac_max <- 0.8
  range_cap_km   <- range_frac_max * max(emp_v$dist) / 1000
  if (range_km > range_cap_km) {
    fit$range[2] <- range_cap_km * 1000
    range_km <- range_cap_km
  }

  list(sub = sub, formula = form, emp_v = emp_v, fit = fit, R2 = R2,
       range_km = range_km, n = nrow(sub))
}

pdf(file.path(out_path, paste0(base_nome, "_variograms.pdf")), width = 8, height = 5)
national_vario <- list()
district_vario <- list()
zona_codes <- sort(unique(legend$zona[!is.na(legend$zona)]))
fit_log <- list()   # one row per scope x variable, for the summary report

for (v in names(kriging_source_col)) {
  src_col <- kriging_source_col[[v]]
  covars  <- trend_covars[[v]]

  nat <- fit_variogram_model(src_col, dat_survey, covars)
  if (is.null(nat)) {
    stop(sprintf(paste(
      "fit_variogram_model() rejected the NATIONAL fit for %s (survey column '%s'):",
      "R2 <= 0 on the full national sample. There is no further fallback below",
      "national scope -- inspect this variable's distribution/trend before retrying."
    ), v, src_col))
  }
  national_vario[[v]] <- nat
  fit_log[[length(fit_log) + 1]] <- data.frame(
    scope = "National", variable = v, R2 = nat$R2, range_km = nat$range_km,
    n_used = nat$n, n_zone_valid = NA_integer_, source = "own fit", stringsAsFactors = FALSE)
  cat(sprintf("NATIONAL | %s (from '%s', covars: %s) -> R2 = %.3f, range = %.1f km, n = %d\n",
              v, src_col, paste(covars, collapse=","), nat$R2, nat$range_km, nat$n))
  print(plot(nat$emp_v, nat$fit,
             main = sprintf("%s [kriging-based spatial fit, survey column '%s'] NATIONAL (fallback)\nrange=%.0f km, R2=%.3f",
                             v, src_col, nat$range_km, nat$R2)))

  district_vario[[v]] <- list()
  for (nc in zona_codes) {
    sub_nc  <- dat_survey[!is.na(dat_survey$zona) & dat_survey$zona == nc, ]
    ok_nc   <- !is.na(sub_nc[[src_col]])
    for (cc in covars) ok_nc <- ok_nc & !is.na(sub_nc[[cc]])
    n_valid <- sum(ok_nc)
    fit_d   <- if (n_valid >= min_obs_district) fit_variogram_model(src_col, sub_nc, covars) else NULL
    zone_fallback <- is.null(fit_d)
    if (is.null(fit_d)) {
      reason <- if (n_valid < min_obs_district)
        sprintf("n = %d < %d valid observations", n_valid, min_obs_district)
      else
        sprintf("n = %d but the district fit was degenerate (R2<=0)", n_valid)
      message(sprintf("District %s: falling back to the national model for %s (%s)", nc, v, reason))
      fit_d <- nat
    } else {
      cat(sprintf("DISTRICT %s | %s -> R2 = %.3f, range = %.1f km, n = %d\n",
                  nc, v, fit_d$R2, fit_d$range_km, fit_d$n))
      print(plot(fit_d$emp_v, fit_d$fit,
                 main = sprintf("%s [kriging-based spatial fit, survey column '%s'] DISTRICT %s\nrange=%.0f km, R2=%.3f",
                                 v, src_col, nc, fit_d$range_km, fit_d$R2)))
    }
    district_vario[[v]][[as.character(nc)]] <- fit_d
    fit_log[[length(fit_log) + 1]] <- data.frame(
      scope = as.character(nc), variable = v, R2 = fit_d$R2, range_km = fit_d$range_km,
      n_used = fit_d$n, n_zone_valid = n_valid,
      source = if (zone_fallback) "national model (fallback)" else "own fit",
      stringsAsFactors = FALSE)
  }
}
dev.off()

# fitting/range lookup tables, one row per soil district (rownames = zona
# code), one column per target variable (vars_cont_y names, i.e. PH /
# SOC) -- used both by K-means initialisation (Phase 5) and by the GA
# optimisation (fitting_by_domain / range_by_domain below).
fitting_by_nuts <- as.data.frame(sapply(vars_cont_y, function(v)
  sapply(district_vario[[v]][as.character(zona_codes)], `[[`, "R2")))
rownames(fitting_by_nuts) <- as.character(zona_codes)

range_by_nuts_km <- as.data.frame(sapply(vars_cont_y, function(v)
  sapply(district_vario[[v]][as.character(zona_codes)], `[[`, "range_km")))
rownames(range_by_nuts_km) <- as.character(zona_codes)

# Guard: R2<=0 would make optimStrata compute zero/NaN within-stratum
# variance, blowing up the GA fitness (Inf/NaN). Floor at a small value.
if (any(fitting_by_nuts <= 0)) {
  message("fitting_by_nuts has non-positive R2 values -- flooring at 0.01 to keep GA fitness finite.")
  fitting_by_nuts[fitting_by_nuts <= 0] <- 0.01
}

cat("\nfitting_by_nuts (R2 per soil district, from point-survey kriging fit):\n");   print(fitting_by_nuts)
cat("range_by_nuts_km (range per soil district, from point-survey kriging fit):\n"); print(range_by_nuts_km)

if (!use_spatial_fit) {
  message("use_spatial_fit = FALSE: overriding fitting/range with fitting = 1, range = 0 in every district.")
  fitting_by_nuts[]   <- 1
  range_by_nuts_km[]  <- 0
}


# ---------------------------------------------------------------------------
# Phase 4 -- SamplingStrata frame construction. Y and variance come from the
# frame itself (predicted PH/SOC and their variance); fitting/range come from
# the per-zone tables built above.
# ---------------------------------------------------------------------------
df2 <- st_drop_geometry(dat)
# Drop cells with missing/non-finite X, Y or variance: optimStrata(method =
# "spatial") refuses to run on them.
fin_cols <- c(vars_cont_x, vars_cont_y, variance_cont_y, "x", "y")
ok_fin   <- Reduce(`&`, lapply(fin_cols, function(cc) is.finite(df2[[cc]])))
if (any(!ok_fin)) message(sprintf("Dropping %d grid cell(s) with non-finite values in: %s",
                                  sum(!ok_fin), paste(fin_cols, collapse = ", ")))
df2 <- df2[ok_fin, ]
df2$id <- c(1:nrow(df2))
# df2$dom_code is already "zona*soil_unit" (built in grid_df above).
dom_code_orig <- as.character(df2$dom_code)   # preserve original "zona*soil_unit" strings before integer conversion
df2$dom_code  <- as.integer(factor(df2$dom_code))

# Check the auxiliary (X) variables for negative values. SamplingStrata
# expects non-negative auxiliary variables, so any variable containing
# negative values is shifted upward (v - min(v)) to make all values >= 0.
neg_x <- vars_cont_x[sapply(vars_cont_x, function(v) any(df2[[v]] < 0, na.rm = TRUE))]
if (length(neg_x) > 0) {
  for (v in neg_x) {
    shift <- min(df2[[v]], na.rm = TRUE)
    warning(sprintf(
      "Variable '%s' contains negative values (min = %.4f): shifting by %.4f to make all values positive.",
      v, shift, -shift))
    df2[[v]] <- df2[[v]] - shift
  }
  message("Shifted auxiliary variable(s) with negative values: ",
          paste(neg_x, collapse = ", "))
} else {
  message("No negative values found in auxiliary variables (vars_cont_x).")
}

# Build the SamplingStrata spatial frame object: X = auxiliary variables,
# Y = predicted PH/SOC values, variance = their prediction variance,
# coordinates = x/y in metres. domainvalue defines the stratification domain.
frame <- buildFrameSpatial(df=df2,
                           id="id",
                           X=vars_cont_x,
                           Y=vars_cont_y,
                           variance=variance_cont_y,
                           lon="x",
                           lat="y",
                           domainvalue = "dom_code")

# Filter out domains with fewer than min_point_per_domain cells: these are too small to
# support meaningful stratification or reliable K-means clustering
dom_code_to_num <- unique(data.frame(
  dom_code     = df2$dom_code,
  domainvalue  = frame$domainvalue,
  dom_code_str = dom_code_orig,       # original string (bridge for range_by_domain lookup)
  stringsAsFactors = FALSE
))

points_per_domain <- tapply(frame$domainvalue, frame$domainvalue, length)
valid_domains <- names(points_per_domain[points_per_domain >= min_point_per_domain])
frame <- frame[frame$domainvalue %in% valid_domains, ]

# Re-index domain values to consecutive integers after filtering, as
# SamplingStrata functions expect a contiguous sequence starting from 1
domain_map <- data.frame(
  old = sort(unique(frame$domainvalue)),
  new = seq_along(sort(unique(frame$domainvalue)))
)
frame$domainvalue <- domain_map$new[match(frame$domainvalue, domain_map$old)]
table(frame$domainvalue)

# Assign each re-indexed domain the fitting (R2) and range of its OWN soil
# district (from fitting_by_nuts / range_by_nuts_km, the kriging-based fit
# on point survey data computed in Phase 3). The soil-district code is the
# first component of dom_code_str (format zona*soil_unit).
new_to_str  <- dom_code_to_num$dom_code_str[match(domain_map$old, dom_code_to_num$domainvalue)]
new_to_nuts <- as.character(sapply(strsplit(new_to_str, "\\*"), `[`, 1))

get_domain_row <- function(nuts_key, tbl) {
  if (nuts_key %in% rownames(tbl)) tbl[nuts_key, , drop = FALSE]
  else as.data.frame(as.list(colMeans(tbl, na.rm = TRUE)))  # defensive fallback, should not trigger
}

fitting_by_domain <- do.call(rbind, lapply(new_to_nuts, get_domain_row, tbl = fitting_by_nuts))
rownames(fitting_by_domain) <- as.character(seq_len(nrow(fitting_by_domain)))

range_by_domain <- do.call(rbind, lapply(new_to_nuts, get_domain_row, tbl = range_by_nuts_km))
rownames(range_by_domain) <- as.character(seq_len(nrow(range_by_domain)))
colnames(range_by_domain) <- paste0("range_", colnames(range_by_domain))

cat("\nfitting_by_domain (per-district R2 from kriging fit, one row per re-indexed domain, head):\n")
print(head(fitting_by_domain))
cat("range_by_domain (per-district range km from kriging fit, one row per re-indexed domain, head):\n")
print(head(range_by_domain))

# Set up precision targets: 5% CV for both PH and SOC in every domain.
# The cv data frame is the input constraints table required by optimStrata.
ndom <- length(unique(frame$domainvalue))
cv_cols <- setNames(
  replicate(length(vars_cont_y), rep(0.05, ndom), simplify = FALSE),
  paste0("CV", seq_along(vars_cont_y))
)
cv <- as.data.frame(c(
  list(DOM = rep("DOM1", ndom)),
  cv_cols,
  list(domainvalue = c(1:ndom))
))
cv

# ---------------------------------------------------------------------------
# K-means spatial clustering to provide strata initial solutions
# Run KmeansSolutionSpatial domain by domain (lapply over domain IDs), using
# each domain's own soil-district kriging-based fitting/range throughout (no
# pre-GA spatial/non-spatial choice here -- the choice between spatial and
# non-spatial settings is made per domain in Phase 6b, after the GA).
# The maximum number of clusters per domain is capped at 5 and limited
# to (number of distinct Y-value combinations - 1) to avoid empty strata.
# A tryCatch wrapper (and an explicit row-count check, since
# KmeansSolutionSpatial can return a degenerate 0-strata result without
# raising an error) handles domains where clustering fails, assigning all
# units to a single stratum as a fallback.
# ---------------------------------------------------------------------------

# Determine the maximum feasible number of clusters per domain:
# at most 5, at least 1, and strictly less than the number of distinct
# combinations of target variable values in that domain
y_cols <- paste0("Y", seq_along(vars_cont_y))
distinct_per_dom <- tapply(seq_len(nrow(frame)), frame$domainvalue, function(idx)
  nrow(unique(frame[idx, y_cols])))
max_k <- pmin(pmax(distinct_per_dom - 1L, 1L), 5L)

# Helper: row-bind a list of data frames that may have different columns,
# filling missing columns with NA to ensure a consistent structure
bind_fill <- function(lst) {
  all_cols <- unique(unlist(lapply(lst, names)))
  aligned <- lapply(lst, function(x) {
    miss <- setdiff(all_cols, names(x))
    if (length(miss) > 0) x[miss] <- NA
    x[all_cols]
  })
  do.call(rbind, aligned)
}

# Sort frame by domain before splitting to ensure consistent subsetting
frame <- frame[order(frame$domainvalue), ]

# Run spatial K-means clustering independently for each domain.
# Zone-specific R² (fitting) and range are passed so that domains within
# the same zone share the same spatial autocorrelation parameters.
kmean_list <- lapply(sort(unique(frame$domainvalue)), function(dv) {
  f_sub <- frame[frame$domainvalue == dv, ]
  f_sub$domainvalue <- 1L         # reset to 1 for the single-domain call
  cv_sub <- cv[cv$domainvalue == dv, , drop = FALSE]
  cv_sub$domainvalue <- 1L
  k_max <- min(max_k[[as.character(dv)]], nrow(f_sub) - 1L)
  k_max <- max(k_max, 1L)
  nuts_dv    <- new_to_nuts[dv]
  fitting_dv <- unlist(get_domain_row(nuts_dv, fitting_by_nuts))
  range_dv   <- unlist(get_domain_row(nuts_dv, range_by_nuts_km)) * 1000   # km -> metres for KmeansSolutionSpatial
  res <- tryCatch(
    KmeansSolutionSpatial(frame       = f_sub,
                          fitting     = fitting_dv,
                          range       = range_dv,
                          kappa       = 1,
                          errors      = cv_sub,
                          nstrata     = NA,
                          minnumstrat = 2,
                          maxclusters = k_max,
                          showPlot    = FALSE,
                          # default maxsamp = 5000 subsamples larger domains, which
                          # would trip the length check below and wrongly trigger
                          # the single-stratum fallback
                          maxsamp     = nrow(f_sub)),
    error = function(e) NULL
  )
  # KmeansSolutionSpatial can return a degenerate result (0 strata / empty
  # suggestions) for some domains instead of raising an error -- validate
  # row count explicitly and fall back to a single stratum in that case too.
  if (is.null(res) || is.null(res$suggestions) || length(res$suggestions) != nrow(f_sub)) {
    cat("  [bypass] dv =", dv, " nrow =", nrow(f_sub), " k_max =", k_max, " -- KmeansSolutionSpatial failed or degenerate\n")
    id_col <- if ("id" %in% names(f_sub)) f_sub$id else seq_len(nrow(f_sub))
    res <- data.frame(id          = as.character(id_col),
                       suggestions = factor(rep("1", nrow(f_sub))),
                       domainvalue = dv,
                       stringsAsFactors = FALSE)
  }
  res$domainvalue <- dv   # restore the original domain identifier
  res
})


# Combine per-domain K-means results into a single data frame
kmean <- bind_fill(kmean_list)

# Derive the number of strata per domain from the K-means cluster suggestions.
# Enforce a minimum of 2 strata per domain (required by optimStrata).
nstrat <- tapply(kmean$suggestions, kmean$domainvalue,
                 FUN = function(x) length(unique(x)))
nstrat <- ifelse(nstrat < 2, 2, nstrat)

# Cap nstrat to the maximum feasible value per domain.
# optimStrata requires each stratum to hold at least minnumstr=2 units,
# so the hard upper bound is floor(n_units_in_domain / 2).
# Exceeding it triggers a "nStrata is not feasible" error in check_spatial_frame.
minnumstr_val <- 2L
units_per_dom <- tapply(frame$domainvalue, frame$domainvalue, length)
max_feasible  <- floor(units_per_dom[names(nstrat)] / minnumstr_val)
nstrat        <- pmin(nstrat, max_feasible)
nstrat        <- pmax(nstrat, minnumstr_val)   # keep minimum of 2

# When use_lucas_points is TRUE, SamplingStrata's internal optimizeStrataSpatial
# appends the forced (census) LUCAS units to framenew/aggr_strata AFTER the GA
# has run, tagging them with STRATO = LABEL = max(nStrata) + 1 (see package
# source: framecens$STRATO <- max(nStrata) + 1) and forcing their SOLUZ to the
# full census N (censtot$SOLUZ <- censtot$N) -- i.e. NOT run through Bethel.
# Phase 6b (below) must exclude these rows before recomputing strata per
# domain, otherwise bethel_cpp() would "optimize" the census stratum instead
# of respecting its forced full inclusion, corrupting the spatial/non-spatial
# comparison for every domain that contains LUCAS points. This mirrors
# exactly how the GA's own reported "Sample size"/"Number of strata" is
# computed (dimens <- sum(round(outstrata$SOLUZ)) is printed BEFORE the
# census rows are appended in the package source).
census_label <- if (use_lucas_points) max(nstrat) + 1L else NA_integer_
# ---------------------------------------------------------------------------
# Genetic algorithm (GA) optimisation of the stratification
#
# The GA runs with fitting = 1 and range = 0 EVERYWHERE (i.e. the non-spatial
# baseline), not with the kriging-based fitting/range. Reasoning: a GA that
# searches under fitting = 1 / range = 0 finds its own optimal non-spatial
# stratification, which is expected to do better than taking a clustering
# optimised for something else (spatial, or K-means only) and merely evaluating
# fitting = 1 / range = 0 on it after the fact. The non-spatial baseline is
# therefore built directly, and Phase 6b (below) only asks, per domain, whether
# the real kriging-based fitting/range does even BETTER on top of that baseline.
# ---------------------------------------------------------------------------
frame <- frame[order(frame$domainvalue),]

fitting_nospatial <- as.data.frame(matrix(1, nrow = nrow(fitting_by_domain), ncol = ncol(fitting_by_domain)))
rownames(fitting_nospatial) <- rownames(fitting_by_domain); colnames(fitting_nospatial) <- colnames(fitting_by_domain)
range_nospatial <- as.data.frame(matrix(0, nrow = nrow(range_by_domain), ncol = ncol(range_by_domain)))
rownames(range_nospatial) <- rownames(range_by_domain); colnames(range_nospatial) <- colnames(range_by_domain)

set.seed(1234)
args <- list(
  method = "spatial",
  errors=cv,
  framesamp=frame,
  iter = 15,           # number of GA generations
  pops = 10,           # population size per generation
  nStrata = nstrat,    # initial number of strata (from K-means)
  fitting = fitting_nospatial,  # fitting=1 everywhere -- see the note above
  range = range_nospatial,      # range=0 everywhere -- see the note above
  kappa = 1,
  writeFiles = FALSE,
  showPlot = FALSE,
  parallel = TRUE)

if (use_lucas_points) args$framecens <- frame_cens
solution <- do.call(optimStrata, args)


#---------------------------------------------------------------------------
# Output and sampling
#---------------------------------------------------------------------------
framenew  <- solution$framenew
outstrata <- solution$aggr_strata

# NOTE on the two totals below when use_lucas_points is TRUE: `outstrata`
# includes the forced-census (LUCAS) rows at their full N (added by the
# package AFTER its own internal "Sample size" print), so sum(outstrata$SOLUZ)
# is NOT directly comparable to the optimized-only totals reported by Phase
# 6b further below (which correctly exclude the census stratum -- see
# census_label). The "optimized-only" total here reproduces exactly what
# optimStrata itself printed during the GA run (dimens/"*** Sample size").
optimized_only <- if (use_lucas_points) outstrata$STRATO != census_label else rep(TRUE, nrow(outstrata))
cat(sprintf("\nGA result (non-spatial baseline, fitting=1/range=0 everywhere): total sample size = %.1f, %d strata (optimized-only, excludes forced census/LUCAS units)\n",
            sum(outstrata$SOLUZ[optimized_only]), sum(optimized_only)))
if (use_lucas_points) {
  cat(sprintf("  (+ %.1f forced census/LUCAS units in %d additional take-all strata; total incl. census = %.1f, %d strata)\n",
              sum(outstrata$SOLUZ[!optimized_only]), sum(!optimized_only), sum(outstrata$SOLUZ), nrow(outstrata)))
}

# ---------------------------------------------------------------------------
# Phase 6b -- POST-HOC per-domain selection: does real spatial correlation
# improve on the non-spatial GA baseline?
#
# For each domain, its OWN final GA-refined stratification (framenew$STRATO,
# already fixed at this point, optimised under fitting=1/range=0) is scored
# under two settings -- fitting=1/range=0 (should reproduce the GA's own
# result, a sanity check) and the real kriging-based fitting/range for that
# domain -- via the exact same fitness formula the GA itself uses internally
# (buildStrataDFSpatial + bethel_cpp, both non-exported/internal functions of
# SamplingStrata, accessed via :::). Whichever setting gives the smaller
# implied allocation on that SAME clustering is kept, replacing that
# domain's rows in the aggregated strata table (outstrata) ONLY where the
# real spatial fitting/range genuinely improves on the non-spatial GA
# result; otherwise that domain's non-spatial GA result is left untouched.
# Because both settings are compared on an identical, already GA-refined
# clustering, this cannot do worse than the non-spatial GA baseline alone --
# Only this post-hoc comparison, on the final GA-refined clustering, is used.
# ---------------------------------------------------------------------------
buildStrataDFSpatial <- SamplingStrata:::buildStrataDFSpatial
bethel_cpp <- SamplingStrata:::bethel_cpp

domain_recompute_strata <- function(dv, fit_vec, range_vec_m) {
  sub <- framenew[framenew$DOMAINVALUE == dv, ]
  # Exclude the forced-census (LUCAS) stratum, if any -- see the note above
  # census_label: those units are never Bethel-optimized (SOLUZ is forced to
  # their full N by the package itself) and are re-added to the final sample
  # independently via the LABEL == max(framenew$LABEL) block further below,
  # so they must not be fed to bethel_cpp() here.
  if (use_lucas_points) sub <- sub[sub$STRATO != census_label, ]
  strato_levels <- sort(unique(sub$STRATO))
  sub$X1 <- match(sub$STRATO, strato_levels)  # re-index to 1..k for buildStrataDFSpatial
  sub$DOMAINVALUE <- 1L
  colnames(sub) <- toupper(colnames(sub))
  strata <- buildStrataDFSpatial(dataset = sub[, c("ID","X1","LON","LAT","Y1","Y2","VAR1","VAR2","DOMAINVALUE")],
                                  fitting = fit_vec, range = range_vec_m, kappa = 1,
                                  progress = FALSE, verbose = FALSE)
  strata$DOM1 <- 1L
  cv_dv <- cv[cv$domainvalue == dv, , drop = FALSE]
  cv_dv$domainvalue <- 1L
  soluz <- tryCatch(bethel_cpp(strata, cv_dv, minnumstr = 2, realAllocation = TRUE),
                     error = function(e) rep(NA_real_, nrow(strata)))
  strata$SOLUZ  <- as.numeric(soluz)
  strata$DOM1   <- dv                        # restore real domain id, for merging into outstrata
  strata$STRATO <- strato_levels[strata$X1]  # map the 1..k index back to the GA's own STRATO labels
  strata
}

# Keep the ORIGINAL (pre-Phase-6b) outstrata around so each domain's forced-
# census row (SOLUZ already correctly = N, set by the package itself) can be
# re-attached below. It is excluded from the bethel_cpp() comparison inside
# domain_recompute_strata() (see census_label), but summaryStrata() further
# down matches framenew's strata to outstrata BY POSITION per domain, so
# outstrata_final must keep the same per-domain row count as framenew --
# dropping the census row here would silently misalign that positional match.
outstrata_raw <- outstrata

domain_choice <- character(0)
outstrata_final <- outstrata[0, ]
# Accumulates the fresh per-domain nospatial recompute (n_nospatial), summed
# across all domains, regardless of which option Phase 6b actually picks.
# This isolates the true spatial benefit from the small (solver-convergence)
# noise between the GA's OWN internal nospatial Bethel allocation
# (outstrata_raw) and this standalone bethel_cpp recompute of the exact same
# fitting=1/range=0 configuration on the exact same clustering: the two
# should coincide almost exactly, but small floating-point/convergence
# differences between independent Bethel calls are expected. Comparing
# nospatial_recomputed_total to outstrata_raw's regular-strata total shows
# that noise; comparing it to outstrata_final's total shows the real,
# noise-free benefit (if any) of using real spatial fitting/range in the
# domains where Phase 6b picked it.
nospatial_recomputed_total <- 0
domain_log <- list()   # one row per domain, for the summary report
for (dv in sort(unique(framenew$DOMAINVALUE))) {
  cat("\nDomain: ",dv)
  fit_real   <- as.numeric(fitting_by_domain[as.character(dv), ])
  rng_real_m <- as.numeric(range_by_domain[as.character(dv), ]) * 1000

  strata_spatial   <- if (use_spatial_fit) domain_recompute_strata(dv, fit_real, rng_real_m) else NULL
  strata_nospatial <- domain_recompute_strata(dv, c(1, 1), c(0, 0))

  # n_spatial/n_nospatial compare only the regular (Bethel-optimized) strata --
  # the census contribution, if any, is identical on both sides (same forced
  # N regardless of fitting/range) so adding it here would not change which
  # option is smaller; it is reattached below, to whichever option wins, so
  # that outstrata_final's structure matches framenew.
  n_spatial   <- if (use_spatial_fit) sum(strata_spatial$SOLUZ, na.rm = TRUE) else NA_real_
  n_nospatial <- sum(strata_nospatial$SOLUZ, na.rm = TRUE)
  use_spatial <- use_spatial_fit && is.finite(n_spatial) && (n_spatial <= n_nospatial || !is.finite(n_nospatial))
  nospatial_recomputed_total <- nospatial_recomputed_total +
    if (is.finite(n_nospatial)) n_nospatial else 0

  # Whenever the spatial version does not win (use_spatial_fit = FALSE, or
  # TRUE but no gain in this domain) keep the GA's own result untouched:
  # strata_nospatial is used only for the comparison above, never as the
  # result, so no solver-convergence noise from the standalone recompute
  # enters the final allocation and TRUE/FALSE give identical output where
  # the spatial version brings no gain.
  winner <- if (use_spatial) strata_spatial
            else outstrata_raw[outstrata_raw$DOM1 == dv & !(outstrata_raw$STRATO %in% census_label), ]
  ga_reg <- outstrata_raw[outstrata_raw$DOM1 == dv & !(outstrata_raw$STRATO %in% census_label), ]
  domain_log[[length(domain_log) + 1]] <- data.frame(
    domain = dv, dom_code = new_to_str[dv],
    fit_PH = fit_real[1], fit_SOC = fit_real[2],
    range_PH_km = rng_real_m[1] / 1000, range_SOC_km = rng_real_m[2] / 1000,
    n_strata_ga = nrow(ga_reg), n_ga = sum(ga_reg$SOLUZ, na.rm = TRUE),
    n_nospatial = n_nospatial, n_spatial = n_spatial,
    n_strata_final = nrow(winner), n_final = sum(winner$SOLUZ, na.rm = TRUE),
    choice = if (use_spatial) "spatial" else "nospatial", stringsAsFactors = FALSE)
  if (use_lucas_points) {
    cens_row <- outstrata_raw[outstrata_raw$DOM1 == dv & outstrata_raw$STRATO == census_label, ]
    if (nrow(cens_row) > 0) winner <- rbind(winner, cens_row[, names(winner)])
  }
  outstrata_final <- rbind(outstrata_final, winner[, names(outstrata_final)])
  domain_choice[as.character(dv)] <- if (use_spatial) "spatial" else "nospatial"
}

cat("\n=== Phase 6b: does real spatial fitting/range improve on the non-spatial GA baseline, per domain? ===\n")
print(table(domain_choice))
# Both totals below include the forced-census/LUCAS units on an equal footing
# (outstrata_final carries the same untouched census rows as outstrata_raw --
# see the reattachment above), so this comparison is now apples-to-apples;
# any difference between the two totals comes only from the regular strata.
cat(sprintf("Total sample size, GA result (non-spatial baseline, fitting=1/range=0 everywhere): %.1f\n", sum(outstrata_raw$SOLUZ)))
cat(sprintf("Total sample size, after Phase 6b (real spatial kept only where it improves):       %.1f\n", sum(outstrata_final$SOLUZ)))
# Diagnostic: separates solver-convergence noise (between the GA's own
# internal nospatial Bethel allocation and this standalone recompute of the
# identical config) from the real, noise-free benefit of using real spatial
# fitting/range where Phase 6b picked it.
census_total <- if (use_lucas_points) sum(outstrata_raw$SOLUZ[outstrata_raw$STRATO == census_label]) else 0
cat(sprintf("Total sample size, fresh recompute of nospatial EVERYWHERE (same clustering, all domains): %.1f\n",
            nospatial_recomputed_total + census_total))
cat(sprintf("  -> vs GA result above: %.1f = solver-convergence noise between the two independent Bethel calls (expect small, not a real effect)\n",
            (nospatial_recomputed_total + census_total) - sum(outstrata_raw$SOLUZ)))
cat(sprintf("  -> vs Phase 6b result above: %.1f = real, noise-free benefit of using spatial fitting/range in the %d domain(s) where Phase 6b picked it\n",
            (nospatial_recomputed_total + census_total) - sum(outstrata_final$SOLUZ),
            sum(domain_choice == "spatial")))

# From here on, `outstrata` is the Phase-6b-corrected table: every domain
# uses whichever setting gave it the smaller allocation on its own final
# clustering. All downstream code (summaryStrata, sampling, outputs) is
# unchanged and simply uses this corrected
# `outstrata`.
outstrata <- outstrata_final

strataStructure <- summaryStrata(framenew, outstrata, progress = FALSE)
sum(strataStructure$Allocation)

cv_out <- as.data.frame(expected_CV(outstrata))
colnames(cv_out) <- paste0("cv_", vars_cont_y)
boxplot(cv_out)

stack_input <- if (use_lucas_points) grid_samp else grid_df

# framenew_samp: link framenew to stack_input through the cell id
framenew$cell <- df2$cell[match(as.integer(framenew$ID), df2$id)]
framenew_samp <- merge(
  subset(framenew, cell %in% stack_input$cell & !duplicated(cell)),
  stack_input[, c("cell", "dom_code", "x", "y",
                  if (use_lucas_grid) "LUCAS_GRID" else NULL)],
  by = "cell"
)

key          <- unique(framenew_samp[, c("DOMAINVALUE", "STRATO")])
outstrata_samp <- merge(outstrata, key,
                        by.x = c("DOM1", "STRATO"),
                        by.y = c("DOMAINVALUE", "STRATO"))

# Spatial sampling
sample1 <- selectSampleSpatial(framenew_samp, outstrata_samp, coord_names = c("x", "y"))

# Sampled points
points1 <- merge(
  sample1[, c("ID", "DOMAINVALUE", "LABEL", "STRATO")],
  framenew_samp[, c("ID", "cell", "dom_code", "x", "y",
                    if (use_lucas_grid) "LUCAS_GRID" else NULL)],
  by = "ID"
) %>%
  mutate(
    LUCAS_GRID = if (use_lucas_grid) LUCAS_GRID else 0L,
    source     = if_else(LUCAS_GRID > 0, "sampled_lucas", "sampled")
  )

# Forced Lucas points
lucas_out <- NULL
if (use_lucas_points) {
  lucas_out <- merge(
    subset(framenew, LABEL == max(framenew$LABEL))[, c("ID", "DOMAINVALUE", "LABEL", "STRATO")],
    lucas_points[, c("id_lucas", "cell", "x", "y", "dom_code", "LUCAS_GRID")],
    by.x = "ID", by.y = "id_lucas"
  )
  lucas_out$source <- "lucas"
}

# Final aggregation
points_all <- bind_rows(points1, lucas_out) %>%
  arrange(desc(source == "lucas")) %>%
  filter(!duplicated(cell))

message("Nrow points_all: ",      nrow(points_all))
message("Lucas forced: ",         sum(points_all$source == "lucas"))
message("Sampled lucas grid: ",   sum(points_all$source == "sampled_lucas"))
message("Sampled: ",              sum(points_all$source == "sampled"))

# Final tab: one row per domain (soil unit); dom_val flags domains with at
# least min_point_per_domain cells (those actually stratified)
tab2 <- grid_df %>%
  count(dom_code, zona, soil_unit, name = "n_point") %>%
  mutate(dom_val = n_point >= min_point_per_domain)

tab3 <- tab2 %>%
  rename(n_cell = n_point) %>%
  left_join(points_all %>% count(dom_code, name = "n_tot"),                                          by = "dom_code") %>%
  left_join(points_all %>% filter(source == "sampled_lucas") %>% count(dom_code, name = "n_lucas_grid"),   by = "dom_code") %>%
  left_join(points_all %>% filter(source == "lucas")         %>% count(dom_code, name = "n_lucas_forced"), by = "dom_code") %>%
  mutate(
    n_tot          = coalesce(n_tot, 0L),
    n_lucas_grid   = coalesce(n_lucas_grid, 0L),
    n_lucas_forced = coalesce(n_lucas_forced, 0L)
  ) %>%
  dplyr::select(
    dom_code, zona, soil_unit,
    n_cell, dom_val,
    n_tot, n_lucas_grid, n_lucas_forced
  )



#--------------------------------------------------------------
#Save files and figures
#--------------------------------------------------------------
#--------------------------------------------------------------
# Summary report (text file, English): main results of the run
#--------------------------------------------------------------
fit_tab <- do.call(rbind, fit_log)
dom_tab <- do.call(rbind, domain_log)
dom_tab$zone      <- sub("\\*.*", "", dom_tab$dom_code)
dom_tab$soil_unit <- sub("^[^*]*\\*", "", dom_tab$dom_code)

rep_lines <- character(0)
add  <- function(...) rep_lines <<- c(rep_lines, paste0(...))
addt <- function(df) {
  op <- options(width = 250)   # keep wide tables on one line
  on.exit(options(op))
  rep_lines <<- c(rep_lines, capture.output(print(df, row.names = FALSE)))
}
hdr  <- function(txt) add("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78))

add("SOIL MONITORING LEVEL (SML) - ITALY SAMPLING SITE EXAMPLE: RESULTS SUMMARY")
add("Run date            : ", format(Sys.time(), "%Y-%m-%d %H:%M"))
add("Scenario            : ", base_nome)
add("Domains             : ", nrow(dom_tab), " soil units in ", length(unique(dom_tab$zone)), " zones")
add("Target variables    : ", paste(vars_cont_y, collapse = ", "), " (CV target = 5% in every domain)")
add("Stratification vars : ", paste(vars_cont_x, collapse = ", "))
add("Forced LUCAS points : ", if (use_lucas_points) sum(points_all$source == "lucas") else 0)
add("use_spatial_fit     : ", use_spatial_fit)

hdr("1. VARIOGRAM FITTING (point-level soil survey, regression kriging, spherical model)")
add("R2 = goodness of fit of the spherical variogram (this is the 'fitting' used in the")
add("spatial variance); range = variogram range in km; trend = x + y + DEM.")
add("A zone with fewer than ", min_obs_district, " valid survey points, or with a degenerate fit,")
add("uses the national model.\n")
fit_out <- fit_tab
fit_out$R2       <- round(fit_out$R2, 3)
fit_out$range_km <- round(fit_out$range_km, 1)
names(fit_out) <- c("scope", "variable", "R2 (fitting)", "range (km)", "n of model", "n valid in zone", "source")
addt(fit_out)

hdr("2. OPTIMISATION WITHOUT SPATIAL COMPONENT (Phase 6: GA, fitting = 1, range = 0)")
add("Sample size per domain (optimised strata only; forced LUCAS census units excluded):\n")
ga_out <- dom_tab[, c("domain", "zone", "soil_unit", "n_strata_ga", "n_ga")]
ga_out$n_ga <- round(ga_out$n_ga, 1)
names(ga_out) <- c("domain", "zone", "soil unit", "strata", "sample size")
addt(ga_out)
census_n <- sum(outstrata_raw$SOLUZ[outstrata_raw$STRATO %in% census_label])
add("\nOptimised strata total : ", round(sum(dom_tab$n_ga), 1), " units in ", sum(dom_tab$n_strata_ga), " strata")
add("Forced LUCAS (census) : ", round(census_n, 1), " units")
add("TOTAL (GA result)     : ", round(sum(outstrata_raw$SOLUZ), 1), " units")

hdr("3. OPTIMISATION WITH SPATIAL COMPONENT (Phase 6b)")
if (!use_spatial_fit) {
  add("Not performed (use_spatial_fit = FALSE): the GA result above is the final allocation.")
} else {
  add("For each domain the strata found by the GA are re-allocated (Bethel) with the zone's")
  add("fitting and range; the spatial version is kept only if it needs no more units than the")
  add("non-spatial one recomputed on the same strata.\n")
  sp_out <- dom_tab[, c("domain", "zone", "soil_unit", "fit_PH", "fit_SOC",
                        "range_PH_km", "range_SOC_km", "n_nospatial", "n_spatial", "choice")]
  sp_out$difference <- round(sp_out$n_spatial - sp_out$n_nospatial, 1)
  sp_out[, c("fit_PH", "fit_SOC")] <- round(sp_out[, c("fit_PH", "fit_SOC")], 3)
  sp_out[, c("range_PH_km", "range_SOC_km")] <- round(sp_out[, c("range_PH_km", "range_SOC_km")], 1)
  sp_out[, c("n_nospatial", "n_spatial")] <- round(sp_out[, c("n_nospatial", "n_spatial")], 1)
  names(sp_out) <- c("domain", "zone", "soil unit", "fit PH", "fit SOC", "range PH km",
                     "range SOC km", "n non-spatial", "n spatial", "choice", "n spatial - n non-sp.")
  addt(sp_out)
  n_ns_all <- nospatial_recomputed_total + census_n
  add("\nDomains where the spatial version was kept : ", sum(dom_tab$choice == "spatial"), " of ", nrow(dom_tab))
  add("TOTAL, GA result (non-spatial)              : ", round(sum(outstrata_raw$SOLUZ), 1))
  add("TOTAL, non-spatial recomputed (same strata) : ", round(n_ns_all, 1),
      "   (difference from GA = ", round(n_ns_all - sum(outstrata_raw$SOLUZ), 1),
      ": solver noise between independent Bethel calls)")
  add("TOTAL, final allocation after Phase 6b      : ", round(sum(outstrata_final$SOLUZ), 1))
  add("Gain of the spatial component               : ", round(n_ns_all - sum(outstrata_final$SOLUZ), 1),
      " units (", round(100 * (n_ns_all - sum(outstrata_final$SOLUZ)) / n_ns_all, 1), "% of the non-spatial total)")
}

hdr("4. WARNINGS")
warn <- character(0)
low_zone <- fit_tab[fit_tab$scope != "National" & fit_tab$R2 < fitting_warn_threshold, ]
if (nrow(low_zone) > 0) {
  warn <- c(warn, sprintf(paste0(
    "LOW FITTING: the variogram fitting (R2) is below %.2f in %d of %d zone x variable fits ",
    "(range %.2f-%.2f). The fitting enters the spatial variance at the denominator, so with a ",
    "low fitting the spatial component is expected to INCREASE the variance, i.e. to give a ",
    "sample size not smaller than the non-spatial one. Treat any spatial gain with caution."),
    fitting_warn_threshold, nrow(low_zone), sum(fit_tab$scope != "National"),
    min(low_zone$R2), max(low_zone$R2)))
  if (use_spatial_fit && any(dom_tab$choice == "spatial")) {
    warn <- c(warn, sprintf(paste0(
      "SPATIAL VERSION SELECTED DESPITE LOW FITTING: in %d of %d domains Phase 6b kept the ",
      "spatial allocation (gain %.1f units). This is not the expected behaviour for fitting < %.2f; ",
      "check how fitting and range enter the spatial variance (the gain may come from the large ",
      "range, ~%.0f-%.0f km on 1 km cells, rather than from the fitting) before relying on it."),
      sum(dom_tab$choice == "spatial"), nrow(dom_tab),
      (nospatial_recomputed_total + census_n) - sum(outstrata_final$SOLUZ), fitting_warn_threshold,
      min(fit_tab$range_km), max(fit_tab$range_km)))
  }
  if (!use_spatial_fit) {
    warn <- c(warn, "The spatial component was NOT used (use_spatial_fit = FALSE), consistent with the low fitting.")
  }
}
fb <- fit_tab[fit_tab$source != "own fit", ]
if (nrow(fb) > 0) {
  warn <- c(warn, sprintf(paste0(
    "NATIONAL MODEL USED AT ZONE LEVEL: for %d of %d zone x variable fits the zone model could not be ",
    "estimated (fewer than %d valid survey points, or degenerate fit) and the national model was used. ",
    "Zones affected: %s. The national spatial structure is not confirmed by local data."),
    nrow(fb), sum(fit_tab$scope != "National"), min_obs_district, paste(unique(fb$scope), collapse = "; ")))
}
if (length(warn) == 0) add("No warnings.") else for (w in warn) add("* ", w, "\n")

hdr("5. FINAL SAMPLE")
add("Total sites selected : ", nrow(points_all))
add("  sampled           : ", sum(points_all$source == "sampled"))
add("  sampled (LUCAS grid): ", sum(points_all$source == "sampled_lucas"))
add("  LUCAS forced      : ", sum(points_all$source == "lucas"))
add("\nExpected CV per domain (target 5%), maximum over domains:")
for (cn in colnames(cv_out)) add("  ", cn, " : ", sprintf("%.4f", max(cv_out[[cn]], na.rm = TRUE)))
add("\nSites per domain (domains with at least ", min_point_per_domain, " cells are stratified):\n")
addt(as.data.frame(tab3))

hdr("6. OTHER NOTES")
add("* Survey points inside the study area are few (see Phase 3): the zone variograms are not reliable.")
add("* The DEM raster is used only for the trend of the variogram fit; the DEM used to build the strata comes from the frame.")
add("* The CV of the target variables is computed on the model-predicted values and variances of the frame.")
add("\nOutput files: ", base_nome, "_units.csv, _allpoints.csv, _sample.csv, _figures.pdf, _variograms.pdf, .gpkg, .RData")
writeLines(rep_lines, paste0(out_path, base_nome, "_summary.txt"))

write.csv(tab3, paste0(out_path, paste0(base_nome, "_units.csv")), row.names=FALSE)
write.csv(points_all, paste0(out_path, paste0(base_nome, "_allpoints.csv")), row.names=FALSE)
write.table(sample1,file.path(out_path, paste0(base_nome, "_sample.csv")), sep=",",row.names=F,quote=F)
save.image(file=paste0(out_path, paste0(base_nome, ".RData")))

# Save domains as geopackage
grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = 3035)
grid_sf <- grid_sf %>%
  left_join(tab3[, c("dom_code", "dom_val", "n_tot", "n_lucas_grid", "n_lucas_forced")],
            by = "dom_code")
grid_sf$dom_code <- grid_sf$soil_unit

st_write(grid_sf, paste0(out_path, paste0(base_nome, ".gpkg")),
         delete_dsn = TRUE)
table(grid_sf$soil_unit)


#################################################################################
# Figures
pdf(paste0(out_path, paste0(base_nome, "_figures.pdf")), width = 8.5, height = 8)

# CV boxplot
boxplot(cv_out,
        col  = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
        main = paste0("Coefficient of Variation Distribution (sample = ", nrow(points_all), ")"),
        ylab = "CV",
        xlab = "")
abline(a = 0.05, b = 0, col = "red")

# Retrieve zona and soil_unit from grid_df through the cell id
points_all$zona      <- grid_df$zona[match(points_all$cell, grid_df$cell)]
points_all$soil_unit <- grid_df$soil_unit[match(points_all$cell, grid_df$cell)]

# Convert to sf object
s_sf1 <- st_as_sf(points_all, coords = c("x", "y"), crs = 3035) %>%
  st_transform(crs = 4326) %>%
  mutate(zona      = as.factor(zona),
         soil_unit = as.factor(soil_unit))

# Map by soil unit
print(
  ggplot(s_sf1) +
    geom_sf(aes(color = soil_unit), size = 0.5, alpha = 0.5) +
    labs(title = paste0("Soil units (sample = ", nrow(points_all), ")"),
         color = "Soil unit") +
    theme_minimal()
)

# Barplot points per domain (stacked)
tab3_long <- tab3 %>%
  filter(dom_val) %>%
  mutate(
    n_sampled_only    = n_tot - n_lucas_grid - n_lucas_forced,
    n_lucas_grid_only = n_lucas_grid
  ) %>%
  dplyr::select(dom_code, n_sampled_only, n_lucas_grid_only, n_lucas_forced) %>%
  tidyr::pivot_longer(cols      = c(n_sampled_only, n_lucas_grid_only, n_lucas_forced),
                      names_to  = "type",
                      values_to = "count") %>%
  mutate(type = factor(type,
                       levels = c("n_sampled_only", "n_lucas_grid_only", "n_lucas_forced"),
                       labels = c("Sampled", "LUCAS Grid", "LUCAS Forced")))
print(
  ggplot(tab3_long, aes(x = dom_code, y = count, fill = type)) +
    geom_col() +
    scale_fill_manual(values = c("Sampled"      = "steelblue",
                                 "LUCAS Grid"   = "darkorange",
                                 "LUCAS Forced" = "red")) +
    labs(title = paste0("Number of Points per Domain (total = ", sum(tab3$n_tot), ")"),
         x     = "Domain",
         y     = "Count",
         fill  = "Point Type") +
    theme_minimal() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank())
)

# Geographic distribution
pts_all_plot <- points_all %>%
  mutate(type = case_when(
    source == "lucas" ~ "LUCAS Forced",
    LUCAS_GRID > 0    ~ "LUCAS Grid",
    TRUE              ~ "Sampled"
  )) %>%
  dplyr::select(x, y, type)

if (use_lucas_points && nrow(lucas_points_filt) > 0) {
  # Convert WGS84 -> EPSG:3035 coordinates outside mutate()
  lf_v <- project(
    vect(lucas_points_filt, geom = c("X_WGS84", "Y_WGS84"), crs = "EPSG:4326"),
    "EPSG:3035"
  )
  lf_coords <- as.data.frame(crds(lf_v))
  colnames(lf_coords) <- c("x", "y")
  lf_coords$type <- "LUCAS Discarded"
  pts_all_plot <- bind_rows(pts_all_plot, lf_coords)
}

pts_all_sf <- st_as_sf(pts_all_plot, coords = c("x", "y"), crs = 3035) %>%
  st_transform(crs = 4326)

print(
  ggplot(pts_all_sf) +
    geom_sf(aes(color = type, size = type, alpha = type)) +
    scale_color_manual(values = c("Sampled"         = "steelblue",
                                  "LUCAS Grid"      = "darkorange",
                                  "LUCAS Forced"    = "red",
                                  "LUCAS Discarded" = "grey40")) +
    scale_size_manual(values  = c("Sampled"         = 0.5,
                                  "LUCAS Grid"      = 0.8,
                                  "LUCAS Forced"    = 1.2,
                                  "LUCAS Discarded" = 0.8)) +
    scale_alpha_manual(values = c("Sampled"         = 0.4,
                                  "LUCAS Grid"      = 0.6,
                                  "LUCAS Forced"    = 1.0,
                                  "LUCAS Discarded" = 0.4)) +
    labs(title = paste0("Geographic Distribution of Points (total = ", nrow(pts_all_plot), ")"),
         color = "Point Type", size = "Point Type", alpha = "Point Type") +
    theme_minimal()
)

dev.off()


# ---------------------------------------------------------------------------
# Export optimized strata as polygons (to be copied to Files/Strata/ of the SML tool)
# ---------------------------------------------------------------------------
res_m <- 1000   # grid resolution (m)

dv_to_str <- dom_code_to_num %>%
  left_join(domain_map, by = c("dom_code" = "old")) %>%
  # dom_code_str is "zona*soil_unit": keep only the soil unit part for the export
  mutate(dom_code_str = sub("^[^*]*\\*", "", dom_code_str)) %>%
  dplyr::select(domainvalue = new, dom_code_str) %>%
  distinct()

strata_polygons <- framenew_samp %>%
  left_join(dv_to_str, by = c("DOMAINVALUE" = "domainvalue")) %>%
  st_as_sf(coords = c("x", "y"), crs = 3035) %>%
  st_buffer(dist = res_m / 2, endCapStyle = "SQUARE", joinStyle = "MITRE") %>%
  group_by(dom_code_str, STRATO) %>%
  summarise(.groups = "drop") %>%
  mutate(stratum_id = paste(dom_code_str, STRATO, sep = "_"))
colnames(strata_polygons)[1] <- "soil_unit"
strata_polygons$stratum_id <- paste(strata_polygons$soil_unit,strata_polygons$STRATO,sep="*")
st_write(strata_polygons,
         file.path(out_path, "strata_SML.gpkg"),
         delete_dsn = TRUE)
