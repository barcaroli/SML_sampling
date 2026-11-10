# =============================================================================
# Script for optimal soil monitoring site selection in Italy (SML)
# EXAMPLE: 6 soil units in 2 zones (Toscana meridionale; Emilia orientale -
# delta del Po), on a 1 km grid
#
# Authors : Claudia Cagnarini, Giulio Barcaroli
# Adapted : 21/10/2026, from the national hybrid script (FAO/JRC raster
#           values + kriging-based spatial structure, 15/09/2026)
# Revised : 10/11/2026, the variogram fitting (kriging) was moved to the
#           separate script Kriging_example_with_6_Italy_soil_units.R; this
#           script reads its results from legend.csv
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
# SPATIAL COMPONENT (flag spatial_component, parameters block)
# ------------------------------------------------------------
#   TRUE  : for every soil unit the spatial parameters fitting (R2) and range
#           (km) of PH and SOC are read from legend.csv, produced by
#           Kriging_example_with_6_Italy_soil_units.R (or edited by hand).
#           They are used in the initial solution (KmeansSolutionSpatial) and
#           in the genetic algorithm (optimStrata); Phase 6b is performed.
#   FALSE : no spatial component. The initial solution is computed with
#           KmeansSolution2, optimStrata runs with fitting = 1 and range = 0
#           in every domain; Phase 6b is not performed and legend.csv is not
#           needed.
# Before setting TRUE, read the report of the kriging script
# (<base_nome>_kriging_report.txt): it says whether the spatial parameters
# are reliable.
#
# INPUT (folder Input/)
# ---------------------
#   frame_SML_1km_6_soil_units.csv   1 km cells of the 6 soil units (frame,
#                                    with fields zona and soil_unit)
#   Soil_observed_points_2022_Italy.xlsx (sheet "Points")
#                                    LUCAS points (X_WGS84, Y_WGS84, SOIL)
#   legend.csv                       only if spatial_component = TRUE: one row
#                                    per soil unit with zona, soil_unit,
#                                    fitting_PH, range_PH_km, fitting_SOC,
#                                    range_SOC_km (output of the kriging
#                                    script; other columns are ignored)
# (Test_6_soil_units_Italy.gpkg, the soil-unit polygons, is not read here: it
# is used by the kriging script and by the SML tool.)
#
# PROCESSING PHASES
# -----------------------------------------------------------------------------
#
# Phase 1  - Set-up and domains. Parameters, frame reading, construction of the
#            grid (grid_df: coordinates x/y, cell id, dom_code = zona*soil_unit,
#            LUCAS_GRID flag from ORIGINE), correlation matrix of PH, SOC, ARG
#            and DEM (PDF).
# Phase 2  - LUCAS points to be retained. Points with SOIL == 1 are projected
#            to EPSG:3035 and linked to their grid cell through POINT_ID;
#            those outside the 6 soil units, or with missing values, are
#            dropped. The retained ones (frame_cens) are forced into the
#            sample as a census stratum, and their cells are removed from the
#            frame that is actually sampled (grid_samp).
# Phase 3  - Spatial parameters. If spatial_component = TRUE, fitting and
#            range of each target variable are read from legend.csv for every
#            soil unit and checked (one row per domain, fitting in (0, 1],
#            range >= 0). If FALSE, fitting = 1 and range = 0 everywhere.
# Phase 4  - SamplingStrata frame. Cells with non-finite values are dropped;
#            the frame is built with buildFrameSpatial (X = ARG, DEM;
#            Y = PH, SOC; variance = PH_VAR, SOC_VAR; domain = soil unit);
#            domains with fewer than min_point_per_domain cells are dropped
#            and domain values are re-indexed; each domain gets the
#            fitting/range of its soil unit; precision constraints are set to
#            CV = 5% on both PH and SOC in every domain.
# Phase 5  - Initial solution, domain by domain (at most 5 clusters per
#            domain), with a single-stratum fallback if it fails:
#            KmeansSolutionSpatial with the soil-unit fitting/range if
#            spatial_component = TRUE, KmeansSolution2 if FALSE.
# Phase 6  - Genetic algorithm (optimStrata, method = "spatial") starting
#            from the K-means number of strata, with the soil-unit
#            fitting/range (TRUE) or with fitting = 1 and range = 0 (FALSE).
#            The forced LUCAS points enter as a census stratum.
# Phase 6b - Only if spatial_component = TRUE: per-domain check of the
#            spatial component. The strata of each domain found by the GA are
#            re-allocated (Bethel) with fitting = 1 / range = 0 and with the
#            soil-unit fitting/range; the GA (spatial) result is kept where
#            the spatial version needs no more units, otherwise the
#            non-spatial allocation of the same strata is used.
# Phase 7  - Sample selection (selectSampleSpatial) and output assembly:
#            sampled points, forced LUCAS points, summary table per domain.
#
# OUTPUT (folder Output/, prefix Example_with_6_soil_units)
# -------------------------------------------------------
#   _correlation_matrix.pdf   correlation matrix of PH, SOC, ARG, DEM
#   _summary.txt              text report: spatial parameters used,
#                             optimisation (Phase 6), spatial check (Phase
#                             6b, if performed), warnings, final sample
#   _units.csv                summary per domain: cells, sampled, LUCAS
#   _allpoints.csv            all selected points (sampled + LUCAS)
#   _sample.csv               SamplingStrata sample (with weights)
#   _figures.pdf              CV boxplot, sample maps, points per domain
#   .gpkg                     1 km grid with domain and sample counts
#   .RData                    full workspace
#   strata_SML.gpkg           strata as polygons (copy it to the Files/Strata/
#                             folder of the SML tool to use it there)
# The variogram plots, legend.csv and the kriging report are produced by
# Kriging_example_with_6_Italy_soil_units.R, not by this script.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Install the required packages, if not yet installed
# -----------------------------------------------------------------------------
# CRAN packages ("parallel" is part of base R and needs no installation)
cran_pkgs <- c("sp", "raster", "ggplot2", "dplyr", "tidyr", "clhs",
               "sf", "terra", "readxl", "doParallel", "pak")
to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

# SamplingStrata: version 1.6 or later is required. If it is not available
# (or older), install the current version from GitHub with pak.
if (!requireNamespace("SamplingStrata", quietly = TRUE) ||
    packageVersion("SamplingStrata") < "1.6") {
  pak::pak("barcaroli/SamplingStrata")
}

library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(tidyr)
library(SamplingStrata)
library(clhs)
library(sf)
library(terra)
library(readxl)
library(parallel)
library(doParallel)


# Modify the following paths:
inp_path <- "D:\\Google Drive\\ISPRA\\SML_sampling\\Input\\"
out_path <- "D:\\Google Drive\\ISPRA\\SML_sampling\\Output\\"
# Spatial parameters per soil unit, written by Kriging_example_with_6_Italy_soil_units.R
# (same path as legend_csv_path there). Read only if spatial_component = TRUE.
legend_csv_path <- paste0(inp_path, "legend.csv")


frame <- read.csv(paste0(inp_path,"frame_SML_1km_6_soil_units.csv"))

# Grid of 1 km cells (the sampling frame): one row per cell,
# with the coordinates/cell id/domain/LUCAS-grid flag the rest of the script
# expects. Domain = soil unit, as "zona*soil_unit": the same key is built
# from the zona and soil_unit columns of legend.csv (Phase 3).
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

# TRUE  = spatial component: fitting/range of each soil unit are read from
#         legend.csv (output of Kriging_example_with_6_Italy_soil_units.R)
#         and used in KmeansSolutionSpatial (Phase 5) and optimStrata
#         (Phase 6); Phase 6b checks, per domain, that the spatial version
#         needs no more units than the non-spatial one on the same strata.
# FALSE = no spatial component: KmeansSolution2 (Phase 5), optimStrata with
#         fitting = 1 and range = 0 (Phase 6), no Phase 6b, legend.csv not read.
#         Use it when the kriging report classifies the spatial parameters as
#         not reliable (in this example: fitting ~0.4-0.5 from the national
#         model, ranges of 80-109 km against mean distances of 11-72 km
#         inside the soil units).
spatial_component <- FALSE

# Fitting (R2) below this value triggers a warning in the summary report
# (_summary.txt) when spatial_component = TRUE: the fitting divides the
# variance of the predicted values in the spatial variance, so a low fitting
# inflates it.
fitting_warn_threshold <- 0.9

vars_cont_y <- c("PH", "SOC")
variance_cont_y <- c("PH_VAR", "SOC_VAR")
vars_cont_x <- c("ARG", "DEM")

use_lucas_points <- TRUE # TRUE if some LUCAS points are imposed (forced into the sample)
use_lucas_grid <- TRUE # TRUE if lucas grid is taken into account in the frame
base_names <- if (use_lucas_grid) c("LUCAS_GRID", "zona", "soil_unit") else c("zona", "soil_unit")

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
# Convert the grid to an sf spatial object
dat <- st_as_sf(grid_df, coords = c("X_LAEA", "Y_LAEA"), crs = 3035, remove = FALSE)

# ---------------------------------------------------------------------------
# Phase 3 -- Spatial parameters per soil unit
#
# spatial_component = TRUE: legend.csv (written by
# Kriging_example_with_6_Italy_soil_units.R, or edited by hand) gives, for
# each soil unit (zona, soil_unit), the columns fitting_<var> (R2, in (0, 1])
# and range_<var>_km (range in km, >= 0) of every target variable in
# vars_cont_y. Other columns (diagnostics of the kriging script) are ignored,
# except "reliability", which is copied to the summary report if present.
# A row with fitting = 1 and range = 0 makes that soil unit non-spatial.
# The tables fitting_by_su / range_by_su_km (one row per dom_code
# "zona*soil_unit", one column per target variable) are used in Phases 4-6b.
#
# spatial_component = FALSE: legend.csv is not read; fitting = 1 and
# range = 0 are used for every domain (Phase 4).
# ---------------------------------------------------------------------------
if (spatial_component) {
  if (!file.exists(legend_csv_path))
    stop("spatial_component = TRUE but legend.csv was not found at: ", legend_csv_path,
         "\nRun Kriging_example_with_6_Italy_soil_units.R first, or set spatial_component <- FALSE.")
  legend_par <- read.csv(legend_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  fit_cols <- paste0("fitting_", vars_cont_y)
  rng_cols <- paste0("range_", vars_cont_y, "_km")
  miss_cols <- setdiff(c("zona", "soil_unit", fit_cols, rng_cols), names(legend_par))
  if (length(miss_cols) > 0)
    stop("legend.csv: missing column(s): ", paste(miss_cols, collapse = ", "))
  legend_par$dom_code <- paste(legend_par$zona, legend_par$soil_unit, sep = "*")
  if (any(duplicated(legend_par$dom_code)))
    stop("legend.csv: duplicated soil unit(s): ",
         paste(unique(legend_par$dom_code[duplicated(legend_par$dom_code)]), collapse = ", "))
  miss_dom <- setdiff(unique(grid_df$dom_code), legend_par$dom_code)
  if (length(miss_dom) > 0)
    stop("legend.csv: no row for the soil unit(s) of the frame: ", paste(miss_dom, collapse = ", "),
         "\n(zona and soil_unit must match the frame exactly)")
  fit_vals <- as.matrix(legend_par[, fit_cols])
  rng_vals <- as.matrix(legend_par[, rng_cols])
  # fitting = 0 is rejected by SamplingStrata (division by zero in the spatial
  # variance); fitting > 1 or negative ranges have no meaning.
  if (any(!is.finite(fit_vals)) || any(fit_vals <= 0 | fit_vals > 1))
    stop("legend.csv: every fitting value must be in (0, 1]")
  if (any(!is.finite(rng_vals)) || any(rng_vals < 0))
    stop("legend.csv: every range value must be >= 0 (km)")

  fitting_by_su <- as.data.frame(fit_vals, row.names = legend_par$dom_code)
  colnames(fitting_by_su) <- vars_cont_y
  range_by_su_km <- as.data.frame(rng_vals, row.names = legend_par$dom_code)
  colnames(range_by_su_km) <- vars_cont_y

  cat("\nSpatial parameters read from legend.csv (fitting = R2, range in km):\n")
  print(cbind(fitting_by_su, setNames(range_by_su_km, paste0("range_", vars_cont_y, "_km"))))
  if ("reliability" %in% names(legend_par) && any(legend_par$reliability != "ADEQUATE"))
    message("WARNING: legend.csv classifies the spatial parameters of some soil units as ",
            paste(unique(legend_par$reliability[legend_par$reliability != "ADEQUATE"]), collapse = " / "),
            ": see the kriging report before relying on the spatial component.")
} else {
  message("spatial_component = FALSE: fitting = 1 and range = 0 in every domain; legend.csv not read.")
}


# ---------------------------------------------------------------------------
# Phase 4 -- SamplingStrata frame construction. Y and variance come from the
# frame itself (predicted PH/SOC and their variance); fitting/range come from
# legend.csv (Phase 3) if spatial_component = TRUE, else fitting = 1, range = 0.
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
  dom_code_str = dom_code_orig,       # original string (bridge for the legend.csv lookup)
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
# unit: from legend.csv (Phase 3) if spatial_component = TRUE, otherwise
# fitting = 1 and range = 0. The key is dom_code_str ("zona*soil_unit").
# fitting_by_domain / range_by_domain_km: one row per re-indexed domain
# (rownames "1", "2", ...), one column per target variable; the GA receives
# the range in metres (range_by_domain_m), like the grid coordinates.
new_to_str <- dom_code_to_num$dom_code_str[match(domain_map$old, dom_code_to_num$domainvalue)]

ny <- length(vars_cont_y)
if (spatial_component) {
  fitting_by_domain  <- fitting_by_su[new_to_str, , drop = FALSE]
  range_by_domain_km <- range_by_su_km[new_to_str, , drop = FALSE]
} else {
  fitting_by_domain  <- as.data.frame(matrix(1, nrow = length(new_to_str), ncol = ny,
                                             dimnames = list(NULL, vars_cont_y)))
  range_by_domain_km <- as.data.frame(matrix(0, nrow = length(new_to_str), ncol = ny,
                                             dimnames = list(NULL, vars_cont_y)))
}
rownames(fitting_by_domain)  <- as.character(seq_len(nrow(fitting_by_domain)))
rownames(range_by_domain_km) <- as.character(seq_len(nrow(range_by_domain_km)))
range_by_domain_m <- range_by_domain_km * 1000

cat("\nfitting_by_domain (one row per re-indexed domain):\n")
print(cbind(dom_code = new_to_str, fitting_by_domain))
cat("range_by_domain_km (one row per re-indexed domain):\n")
print(cbind(dom_code = new_to_str, range_by_domain_km))

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
# Phase 5 -- K-means clustering to provide strata initial solutions
# Run domain by domain (lapply over domain IDs):
#   spatial_component = TRUE : KmeansSolutionSpatial, with the domain's own
#                              soil-unit fitting/range from legend.csv;
#   spatial_component = FALSE: KmeansSolution2 (non-spatial; it clusters on
#                              the Y values and evaluates each clustering
#                              with the non-spatial Bethel allocation).
# Only the resulting number of strata per domain is passed to the GA.
# The maximum number of clusters per domain is capped at 5 and limited
# to (number of distinct Y-value combinations - 1) to avoid empty strata.
# A tryCatch wrapper (and an explicit row-count check, since the K-means
# functions can return a degenerate 0-strata result without raising an
# error) handles domains where clustering fails, assigning all units to a
# single stratum as a fallback.
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

# Run K-means clustering independently for each domain.
kmean_list <- lapply(sort(unique(frame$domainvalue)), function(dv) {
  f_sub <- frame[frame$domainvalue == dv, ]
  f_sub$domainvalue <- 1L         # reset to 1 for the single-domain call
  cv_sub <- cv[cv$domainvalue == dv, , drop = FALSE]
  cv_sub$domainvalue <- 1L
  k_max <- min(max_k[[as.character(dv)]], nrow(f_sub) - 1L)
  k_max <- max(k_max, 1L)
  if (spatial_component) {
    fitting_dv <- as.numeric(fitting_by_domain[as.character(dv), ])
    range_dv   <- as.numeric(range_by_domain_m[as.character(dv), ])   # metres, like the coordinates
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
  } else {
    # Non-spatial frame for KmeansSolution2: id, X, Y and domain only
    # (variances and coordinates are not used by the non-spatial functions).
    x_cols  <- grep("^X[0-9]+$", names(f_sub), value = TRUE)
    f_sub_ns <- f_sub[, c("id", x_cols, y_cols, "domainvalue")]
    res <- tryCatch(
      KmeansSolution2(frame       = f_sub_ns,
                      errors      = cv_sub,
                      nstrata     = NA,
                      minnumstrat = 2,
                      maxclusters = k_max,
                      showPlot    = FALSE),
      error = function(e) NULL
    )
  }
  # The K-means functions can return a degenerate result (0 strata / empty
  # suggestions) for some domains instead of raising an error -- validate
  # row count explicitly and fall back to a single stratum in that case too.
  if (is.null(res) || is.null(res$suggestions) || length(res$suggestions) != nrow(f_sub)) {
    cat("  [bypass] dv =", dv, " nrow =", nrow(f_sub), " k_max =", k_max, " -- K-means failed or degenerate\n")
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
# Phase 6 -- Genetic algorithm (GA) optimisation of the stratification
#
# spatial_component = TRUE : the GA runs with the soil-unit fitting/range read
#                            from legend.csv (fitting_by_domain /
#                            range_by_domain_m), i.e. it searches the optimal
#                            stratification under the spatial variance.
# spatial_component = FALSE: the GA runs with fitting = 1 and range = 0 in
#                            every domain (non-spatial).
# ---------------------------------------------------------------------------
frame <- frame[order(frame$domainvalue),]

set.seed(1234)
args <- list(
  method = "spatial",
  errors=cv,
  framesamp=frame,
  iter = 15,           # number of GA generations
  pops = 10,           # population size per generation
  nStrata = nstrat,    # initial number of strata (from K-means)
  fitting = fitting_by_domain,   # legend.csv values, or 1 everywhere (Phase 4)
  range = range_by_domain_m,     # legend.csv values in metres, or 0 everywhere (Phase 4)
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
ga_label  <- if (spatial_component) "spatial, fitting/range from legend.csv" else "non-spatial, fitting=1/range=0 everywhere"

# NOTE on the two totals below when use_lucas_points is TRUE: `outstrata`
# includes the forced-census (LUCAS) rows at their full N (added by the
# package AFTER its own internal "Sample size" print), so sum(outstrata$SOLUZ)
# is NOT directly comparable to the optimized-only totals (which exclude the
# census stratum -- see census_label). The "optimized-only" total here
# reproduces exactly what optimStrata itself printed during the GA run
# (dimens/"*** Sample size").
optimized_only <- if (use_lucas_points) outstrata$STRATO != census_label else rep(TRUE, nrow(outstrata))
cat(sprintf("\nGA result (%s): total sample size = %.1f, %d strata (optimized-only, excludes forced census/LUCAS units)\n",
            ga_label, sum(outstrata$SOLUZ[optimized_only]), sum(optimized_only)))
if (use_lucas_points) {
  cat(sprintf("  (+ %.1f forced census/LUCAS units in %d additional take-all strata; total incl. census = %.1f, %d strata)\n",
              sum(outstrata$SOLUZ[!optimized_only]), sum(!optimized_only), sum(outstrata$SOLUZ), nrow(outstrata)))
}

# Keep the ORIGINAL (GA) outstrata: it is the reference of the report and,
# with spatial_component = TRUE, the starting point of Phase 6b.
outstrata_raw <- outstrata
census_n <- sum(outstrata_raw$SOLUZ[outstrata_raw$STRATO %in% census_label])

# ---------------------------------------------------------------------------
# Phase 6b -- (only if spatial_component = TRUE) POST-HOC per-domain check:
# does the spatial component reduce the sample size on the GA strata?
#
# For each domain, its OWN final GA stratification (framenew$STRATO, already
# fixed at this point, optimised under the spatial setting) is re-allocated
# under two settings -- the soil-unit fitting/range from legend.csv (should
# reproduce the GA's own result, a sanity check) and fitting = 1 / range = 0
# (non-spatial) -- via the exact same fitness formula the GA itself uses
# internally (buildStrataDFSpatial + bethel_cpp, both non-exported/internal
# functions of SamplingStrata, accessed via :::). Both recomputed allocations
# are compared on the SAME strata:
#   - spatial needs no more units  -> the GA result of the domain is kept
#                                     untouched (no solver noise from the
#                                     standalone recompute enters the result);
#   - non-spatial needs fewer units -> the domain's rows of outstrata are
#                                     replaced by the non-spatial allocation of
#                                     the same strata.
# With spatial_component = FALSE this phase is skipped: the GA result is the
# final allocation.
# ---------------------------------------------------------------------------
domain_choice <- character(0)
domain_log <- list()   # one row per domain, for the summary report

if (spatial_component) {
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
    strata <- buildStrataDFSpatial(dataset = sub[, c("ID", "X1", "LON", "LAT",
                                                     paste0("Y", seq_len(ny)), paste0("VAR", seq_len(ny)),
                                                     "DOMAINVALUE")],
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

  # Each domain's forced-census row (SOLUZ already correctly = N, set by the
  # package itself) is re-attached below. It is excluded from the bethel_cpp()
  # comparison inside domain_recompute_strata() (see census_label), but
  # summaryStrata() further down matches framenew's strata to outstrata BY
  # POSITION per domain, so outstrata_final must keep the same per-domain row
  # count as framenew -- dropping the census row would silently misalign that
  # positional match.
  outstrata_final <- outstrata[0, ]
  # spatial_recomputed_total accumulates the standalone spatial recompute of
  # every domain: compared with the GA total it shows the solver-convergence
  # noise between the GA's own internal Bethel allocation and an independent
  # bethel_cpp call on the identical configuration (expected to be small).
  spatial_recomputed_total   <- 0
  nospatial_recomputed_total <- 0
  for (dv in sort(unique(framenew$DOMAINVALUE))) {
    cat("\nDomain: ", dv)
    fit_real   <- as.numeric(fitting_by_domain[as.character(dv), ])
    rng_real_m <- as.numeric(range_by_domain_m[as.character(dv), ])

    strata_spatial   <- domain_recompute_strata(dv, fit_real, rng_real_m)
    strata_nospatial <- domain_recompute_strata(dv, rep(1, ny), rep(0, ny))

    # n_spatial/n_nospatial compare only the regular (Bethel-optimized) strata --
    # the census contribution, if any, is identical on both sides (same forced
    # N regardless of fitting/range); it is reattached below to the winner.
    n_spatial   <- sum(strata_spatial$SOLUZ, na.rm = TRUE)
    n_nospatial <- sum(strata_nospatial$SOLUZ, na.rm = TRUE)
    use_spatial <- is.finite(n_spatial) && (n_spatial <= n_nospatial || !is.finite(n_nospatial))
    spatial_recomputed_total   <- spatial_recomputed_total   + if (is.finite(n_spatial)) n_spatial else 0
    nospatial_recomputed_total <- nospatial_recomputed_total + if (is.finite(n_nospatial)) n_nospatial else 0

    ga_reg <- outstrata_raw[outstrata_raw$DOM1 == dv & !(outstrata_raw$STRATO %in% census_label), ]
    winner <- if (use_spatial) ga_reg else strata_nospatial
    domain_log[[length(domain_log) + 1]] <- data.frame(
      domain = dv, dom_code = new_to_str[dv],
      setNames(as.list(fit_real), paste0("fit_", vars_cont_y)),
      setNames(as.list(rng_real_m / 1000), paste0("range_", vars_cont_y, "_km")),
      n_strata_ga = nrow(ga_reg), n_ga = sum(ga_reg$SOLUZ, na.rm = TRUE),
      n_nospatial = n_nospatial, n_spatial = n_spatial,
      n_strata_final = nrow(winner), n_final = sum(winner$SOLUZ, na.rm = TRUE),
      choice = if (use_spatial) "spatial" else "nospatial", stringsAsFactors = FALSE)
    if (use_lucas_points) {
      cens_row <- outstrata_raw[outstrata_raw$DOM1 == dv & outstrata_raw$STRATO == census_label, ]
      if (nrow(cens_row) > 0) winner <- rbind(winner[, names(outstrata_final)], cens_row[, names(outstrata_final)])
    }
    outstrata_final <- rbind(outstrata_final, winner[, names(outstrata_final)])
    domain_choice[as.character(dv)] <- if (use_spatial) "spatial" else "nospatial"
  }

  cat("\n=== Phase 6b: does the spatial component reduce the sample size on the GA strata, per domain? ===\n")
  print(table(domain_choice))
  # All totals below include the forced-census/LUCAS units on an equal footing.
  cat(sprintf("Total sample size, GA result (spatial):                              %.1f\n", sum(outstrata_raw$SOLUZ)))
  cat(sprintf("Total sample size, spatial recomputed on the same strata:            %.1f   (difference = solver noise)\n",
              spatial_recomputed_total + census_n))
  cat(sprintf("Total sample size, non-spatial (fitting=1/range=0) on the same strata: %.1f\n",
              nospatial_recomputed_total + census_n))
  cat(sprintf("Total sample size, after Phase 6b:                                   %.1f\n", sum(outstrata_final$SOLUZ)))

  # From here on, `outstrata` is the Phase-6b table: every domain uses
  # whichever setting gave it the smaller allocation on its own GA strata.
  # All downstream code (summaryStrata, sampling, outputs) simply uses it.
  outstrata <- outstrata_final
} else {
  # No Phase 6b: the GA (non-spatial) result is the final allocation.
  for (dv in sort(unique(framenew$DOMAINVALUE))) {
    ga_reg <- outstrata_raw[outstrata_raw$DOM1 == dv & !(outstrata_raw$STRATO %in% census_label), ]
    domain_log[[length(domain_log) + 1]] <- data.frame(
      domain = dv, dom_code = new_to_str[dv],
      n_strata_ga = nrow(ga_reg), n_ga = sum(ga_reg$SOLUZ, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  outstrata_final <- outstrata_raw
}

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
add("spatial_component   : ", spatial_component)

hdr("1. SPATIAL PARAMETERS (Phase 3)")
if (spatial_component) {
  add("Read from ", legend_csv_path, " (produced by Kriging_example_with_6_Italy_soil_units.R:")
  add("see its report for the fit and the reliability of these values). fitting = R2;")
  add("range in km. Used in KmeansSolutionSpatial (Phase 5) and optimStrata (Phase 6).\n")
  par_out <- data.frame(domain = seq_along(new_to_str), zone = sub("\\*.*", "", new_to_str),
                        soil_unit = sub("^[^*]*\\*", "", new_to_str),
                        round(setNames(fitting_by_domain, paste0("fit_", vars_cont_y)), 3),
                        round(setNames(range_by_domain_km, paste0("range_", vars_cont_y, "_km")), 1),
                        check.names = FALSE)
  if ("reliability" %in% names(legend_par))
    par_out$reliability <- legend_par$reliability[match(new_to_str, legend_par$dom_code)]
  addt(par_out)
} else {
  add("Not used: fitting = 1 and range = 0 in every domain (non-spatial). Initial solution by")
  add("KmeansSolution2 (Phase 5); legend.csv not read.")
}

hdr(paste0("2. OPTIMISATION (Phase 6: GA, ", ga_label, ")"))
add("Sample size per domain (optimised strata only; forced LUCAS census units excluded):\n")
ga_out <- dom_tab[, c("domain", "zone", "soil_unit", "n_strata_ga", "n_ga")]
ga_out$n_ga <- round(ga_out$n_ga, 1)
names(ga_out) <- c("domain", "zone", "soil unit", "strata", "sample size")
addt(ga_out)
add("\nOptimised strata total : ", round(sum(dom_tab$n_ga), 1), " units in ", sum(dom_tab$n_strata_ga), " strata")
add("Forced LUCAS (census) : ", round(census_n, 1), " units")
add("TOTAL (GA result)     : ", round(sum(outstrata_raw$SOLUZ), 1), " units")

hdr("3. CHECK OF THE SPATIAL COMPONENT (Phase 6b)")
if (!spatial_component) {
  add("Not performed (spatial_component = FALSE): the GA result above is the final allocation.")
} else {
  add("For each domain the strata found by the GA are re-allocated (Bethel) with the soil-unit")
  add("fitting and range and with fitting = 1 / range = 0; the GA (spatial) result is kept where")
  add("the spatial version needs no more units, otherwise the non-spatial allocation of the same")
  add("strata is used.\n")
  sp_cols <- c("domain", "zone", "soil_unit", "n_nospatial", "n_spatial", "choice")
  sp_out <- dom_tab[, sp_cols]
  sp_out$difference <- round(sp_out$n_spatial - sp_out$n_nospatial, 1)
  sp_out[, c("n_nospatial", "n_spatial")] <- round(sp_out[, c("n_nospatial", "n_spatial")], 1)
  names(sp_out) <- c("domain", "zone", "soil unit", "n non-spatial", "n spatial", "choice", "n spatial - n non-sp.")
  addt(sp_out)
  n_ns_all <- nospatial_recomputed_total + census_n
  n_sp_all <- spatial_recomputed_total + census_n
  add("\nDomains where the spatial version was kept : ", sum(dom_tab$choice == "spatial"), " of ", nrow(dom_tab))
  add("TOTAL, GA result (spatial)                  : ", round(sum(outstrata_raw$SOLUZ), 1))
  add("TOTAL, spatial recomputed (same strata)     : ", round(n_sp_all, 1),
      "   (difference from GA = ", round(n_sp_all - sum(outstrata_raw$SOLUZ), 1),
      ": solver noise between independent Bethel calls)")
  add("TOTAL, non-spatial on the same strata       : ", round(n_ns_all, 1))
  add("TOTAL, final allocation after Phase 6b      : ", round(sum(outstrata_final$SOLUZ), 1))
  add("Reduction due to the spatial component      : ", round(n_ns_all - sum(outstrata_final$SOLUZ), 1),
      " units (", round(100 * (n_ns_all - sum(outstrata_final$SOLUZ)) / n_ns_all, 1),
      "% of the non-spatial total on the same strata)")
}

hdr("4. WARNINGS")
warn <- character(0)
if (spatial_component) {
  fit_mat  <- as.matrix(fitting_by_domain)
  low_dom  <- which(apply(fit_mat < fitting_warn_threshold, 1, any))
  if (length(low_dom) > 0) {
    warn <- c(warn, sprintf(paste0(
      "LOW FITTING: the fitting (R2) read from legend.csv is below %.2f in %d of %d domains ",
      "(values %.2f-%.2f). The fitting divides the variance of the predicted values in the ",
      "spatial variance, so a low fitting inflates it. Treat any spatial gain with caution."),
      fitting_warn_threshold, length(low_dom), nrow(fit_mat),
      min(fit_mat[low_dom, ]), max(fit_mat[low_dom, ])))
    if (any(dom_tab$choice == "spatial")) {
      warn <- c(warn, sprintf(paste0(
        "SPATIAL VERSION KEPT DESPITE LOW FITTING: in %d of %d domains Phase 6b kept the spatial ",
        "allocation (reduction %.1f units). With a range of %.0f-%.0f km on 1 km cells the ",
        "prediction errors of the cells of a stratum are treated as strongly correlated and the ",
        "prediction-variance term of the stratum variance nearly cancels: the reduction comes from ",
        "the range, not from the fitting. Rely on it only if the range is supported by local data ",
        "(see the kriging report)."),
        sum(dom_tab$choice == "spatial"), nrow(dom_tab),
        (nospatial_recomputed_total + census_n) - sum(outstrata_final$SOLUZ),
        min(range_by_domain_km), max(range_by_domain_km)))
    }
  }
  if ("reliability" %in% names(legend_par)) {
    rel_dom <- legend_par$reliability[match(new_to_str, legend_par$dom_code)]
    bad <- rel_dom != "ADEQUATE"
    if (any(bad)) {
      warn <- c(warn, sprintf(paste0(
        "SPATIAL PARAMETERS NOT VALIDATED: legend.csv classifies %d of %d domains as %s (%s). ",
        "The kriging report recommends spatial_component = FALSE when any soil unit is NOT RELIABLE."),
        sum(bad), length(bad), paste(unique(rel_dom[bad]), collapse = " / "),
        paste(sub("^[^*]*\\*", "", new_to_str[bad]), collapse = ", ")))
    }
  }
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
if (spatial_component) {
  add("* The spatial parameters come from legend.csv; the variogram fits, their reliability and the")
  add("  variogram plots are in the outputs of Kriging_example_with_6_Italy_soil_units.R.")
}
add("* The CV of the target variables is computed on the model-predicted values and variances of the frame.")
add("\nOutput files: ", base_nome, "_units.csv, _allpoints.csv, _sample.csv, _figures.pdf, _correlation_matrix.pdf, .gpkg, .RData, strata_SML.gpkg")
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
