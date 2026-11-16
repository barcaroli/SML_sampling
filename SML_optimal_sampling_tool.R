# =============================================================================
# Script for optimal soil monitoring site selection according to the Soil Monitoring Law
# 
#
# Authors : Claudia Cagnarini, Giulio Barcaroli
# Latest revision : 11/11/2026
# =============================================================================
#
# PURPOSE
#
# The sample is optimised with SamplingStrata (optimal stratification by
# genetic algorithm + Bethel allocation) so that the target variables reach 
# a target CV (cv_target) in every domain. The soil units are the estimation domains;
# the target and auxiliary variables are provided as DSM. 
# Optionally, the spatial component can be activated to further optimize
# the sampling.
#
#
# SPATIAL COMPONENT (flag spatial_component, parameters block)
# ------------------------------------------------------------
#   TRUE  : for every soil unit the spatial parameters range (km) and nugget
#           of each target variable are read from the legend.
#           They are used in the initial solution (KmeansSolutionSpatial) and
#           in the genetic algorithm (optimStrata); Phase 5b reports the
#           reduction of the sample size due to the spatial component.
#   FALSE : no spatial component. The initial solution is computed with
#           KmeansSolution2, optimStrata runs with range = 0 in every domain;
#           Phase 5b is not performed and the legend is not needed.
#
# Both parameters describe the spatial correlation of the prediction errors
# of the DSM, and are estimated, for each soil unit (or group of similar soil
# units) and each target variable, on the variogram of the out-of-sample
# residuals (observed - predicted), fitted with an exponential model:
#   range  : range parameter a of the exponential model, in km (>= 0)
#            (not the practical range, about 3a: the script uses kappa = 1)
#   nugget : relative nugget c0 / (c0 + c1), in [0, 1]
# The variance in the strata is the additive one (de Gruijter et al. 2016):
# squared differences of the predictions plus prediction-error variances
# minus their covariance (1 - nugget) s_i s_j exp(-d / range). No R2
# ("fitting") is used: dividing by R2 as well would count the error variance
# twice. Requires SamplingStrata with the 'nugget' parameter (version > 1.6).
#
# PROCESSING PHASES
# -----------------------------------------------------------------------------
#
# Phase 1  - Set-up and domains. Parameters, frame reading, construction of the
#            grid.
# Phase 2  - Spatial parameters. If spatial_component = TRUE, range and
#            nugget of each target variable are read from the legend for every
#            soil unit and checked (one row per domain, range >= 0, nugget in
#            [0, 1]). If FALSE, range = 0 everywhere.
# Phase 3  - SamplingStrata frame. Cells with non-finite values are dropped;
#            the frame is built with buildFrameSpatial;
#            domains with fewer than min_point_per_domain cells are dropped
#            and domain values are re-indexed; each domain gets the
#            range/nugget of its soil unit; precision constraints are set to
#            CV = cv_target in every domain.
# Phase 4  - Initial solution, domain by domain (at most max_clusters
#            clusters per domain), with a single-stratum fallback if it fails:
#            KmeansSolutionSpatial with the soil-unit range/nugget if
#            spatial_component = TRUE, KmeansSolution2 if FALSE.
# Phase 5  - Genetic algorithm (optimStrata, method = "spatial") starting
#            from the K-means number of strata, with the soil-unit
#            range/nugget (TRUE) or with range = 0 (FALSE).
# Phase 5b - Only if spatial_component = TRUE, report only: the strata found
#            by the GA are re-allocated (Bethel) without the spatial
#            component (range = 0), to show the reduction of the sample size
#            due to the spatial component. The GA result is the final
#            allocation. (The former selection between spatial and
#            non-spatial solutions is no longer needed: with the additive
#            variance the spatial component can only reduce the variance.)
# Phase 6  - Sample selection (selectSampleSpatial) and output assembly:
#            sampled points, summary table per domain, report, figures.
#
#
# INPUT (folder Input/)
# ---------------------
#   Frame_1km_6_soil_units_Italy.csv     sampling frame: one row per 1 km cell,
#                                        with the coordinates (X_LAEA, Y_LAEA,
#                                        EPSG:3035), the soil unit (soil_unit),
#                                        the DSM predictions of the target
#                                        variables with their prediction 
#                                        variance, and the auxiliary variables. 
#                                        The cells of the LUCAS 2 km grid have 
#                                        POINT_ID > 0
#   Test_6_soil_units_Italy_legend.csv   legend, only if spatial_component = TRUE:
#                                        one row per soil unit with columns
#                                        soil_unit and, for each target variable
#                                        <var>, range_<var>_km (range parameter
#                                        of the exponential variogram, km) and
#                                        nugget_<var> (relative nugget)
#   Test_6_soil_units_Italy.gpkg         soil-unit polygons: not read here, used
#                                        by the SML tool
# 
# OUTPUT (folder Output/, prefix Example_with_6_soil_units)
# -------------------------------------------------------
#   _summary.txt              text report: spatial parameters used,
#                             optimisation (Phase 5), reduction due to the
#                             spatial component (Phase 5b, if performed),
#                             final sample
#   _units.csv                one row per soil unit: number of cells, whether
#                             it was stratified (at least min_point_per_domain
#                             cells), number of selected points, of which on
#                             the LUCAS grid
#   _allpoints.csv            one row per selected point: domain, stratum, cell
#                             id, coordinates (EPSG:3035), LUCAS grid flag
#   _sample.csv               SamplingStrata sample (with weights)
#   _figures.pdf              correlation matrix, CV boxplot, sample maps,
#                             points per domain
#   .gpkg                     1 km grid with domain and sample counts
#   .RData                    full workspace
#   _strata_SML.gpkg          strata as polygons (copy it to the Files/Strata/
#                             folder of the SML tool to use it there)
#
#
# =============================================================================
# -----------------------------------------------------------------------------
# Install the required packages, if not yet installed
# -----------------------------------------------------------------------------
cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "sf", "pak")
to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

# SamplingStrata: a version with the 'nugget' parameter of optimStrata is
# required. If it is not available, install the current version from GitHub
# with pak.
if (!requireNamespace("SamplingStrata", quietly = TRUE) ||
    !"nugget" %in% names(formals(SamplingStrata::optimStrata))) {
  pak::pak("barcaroli/SamplingStrata")
}


library(ggplot2)
library(dplyr)
library(tidyr)
library(SamplingStrata)
library(sf)


# Modify the following paths:
inp_path <- "insert path\\Input\\"
out_path <- "insert path\\Output\\"
inp_path <- "D:\\Google Drive\\ISPRA\\SML_optimal_sampling_tool_IT\\Input\\"
out_path <- "D:\\Google Drive\\ISPRA\\SML_optimal_sampling_tool_IT\\Output\\"

frame_in <- read.csv(paste0(inp_path, "Frame_1km_6_soil_units_Italy.csv"))
legend_csv_path <- paste0(inp_path, "Test_6_soil_units_Italy_legend.csv")


#----------------------------------------------------------------------------------------------------
# Parameters & settings (to be modified by users)
#----------------------------------------------------------------------------------------------------
vars_cont_y     <- c("PH", "SOC")          # target variables
variance_cont_y <- c("PH_VAR", "SOC_VAR")  # prediction variance of the target variables
vars_cont_x     <- c("CLAY", "DEM")         # auxiliary variables

min_point_per_domain <- 2   # soil units with fewer cells are dropped

cv_target              <- 0.05   # target CV of every target variable in every domain
min_sample_per_stratum <- 2      # minimum number of sampled units per stratum (>= 2)
max_clusters           <- 5      # maximum number of K-means clusters per domain (Phase 4)

# Genetic algorithm (Phase 5): low values for quick tests; use much higher
# values (e.g. iter in the hundreds, pops 20-50) for final results
ga_iter <- 15   # number of generations
ga_pops <- 10   # population size per generation

# TRUE = range/nugget read from the legend (Phases 4, 5, 5b);
# FALSE = range = 0 in every domain, the legend is not read
spatial_component <- FALSE

use_lucas_grid <- TRUE   # TRUE to flag the cells of the LUCAS grid (POINT_ID > 0)

base_name <- "Test_6_soil_units_Italy"   # prefix of the output files



#----------------------------------------------------------------------------------------------------
# Phase 1 -- Set-up and domains
#----------------------------------------------------------------------------------------------------
# Sampling frame: one row per 1 km cell. Domain = soil unit.
grid_df <- frame_in
grid_df$x          <- grid_df$X_LAEA
grid_df$y          <- grid_df$Y_LAEA
grid_df$cell       <- seq_len(nrow(grid_df))
grid_df$dom_code   <- grid_df$soil_unit
grid_df$LUCAS_GRID <- if (use_lucas_grid) as.integer(grid_df$POINT_ID > 0) else 0L

#
# Frame creation
#
# Convert the grid to an sf spatial object
dat <- st_as_sf(grid_df, coords = c("X_LAEA", "Y_LAEA"), crs = 3035, remove = FALSE)

# ---------------------------------------------------------------------------
# Phase 2 -- Spatial parameters per soil unit
#
# spatial_component = TRUE: the legend gives, for each soil unit (soil_unit),
# the columns range_<var>_km (range parameter of the exponential variogram of
# the prediction errors, in km, >= 0) and nugget_<var> (relative nugget
# c0/(c0+c1), in [0, 1]) of every target variable in vars_cont_y. A missing
# nugget column is taken as 0. Other columns are ignored, except fitting_<var>
# (old format): values < 1 are rejected, because dividing by R2 together with
# the prediction-error variances would count the error variance twice.
# A row with range = 0 (or nugget = 1) makes that soil unit non-spatial.
# The tables range_by_su_km / nugget_by_su (one row per dom_code
# = soil_unit, one column per target variable) are used in Phases 3-5b.
#
# spatial_component = FALSE: the legend is not read; range = 0 is used for
# every domain (Phase 3).
# ---------------------------------------------------------------------------
if (spatial_component) {
  if (!file.exists(legend_csv_path))
    stop("spatial_component = TRUE but the legend was not found at: ", legend_csv_path,
         "\nProvide the legend in the input folder, or set spatial_component <- FALSE.")
  legend_par <- read.csv(legend_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  rng_cols <- paste0("range_", vars_cont_y, "_km")
  nug_cols <- paste0("nugget_", vars_cont_y)
  fit_cols <- paste0("fitting_", vars_cont_y)
  miss_cols <- setdiff(c("soil_unit", rng_cols), names(legend_par))
  if (length(miss_cols) > 0)
    stop("Legend: missing column(s): ", paste(miss_cols, collapse = ", "))
  for (nc in setdiff(nug_cols, names(legend_par))) {
    message("Legend: column ", nc, " not found, nugget = 0 assumed")
    legend_par[[nc]] <- 0
  }
  old_fit <- intersect(fit_cols, names(legend_par))
  if (length(old_fit) > 0 && any(as.matrix(legend_par[, old_fit]) < 1, na.rm = TRUE))
    stop("Legend: fitting values < 1 are no longer used. The prediction-error variances ",
         "of the frame already account for the levelling of the predictions, and dividing ",
         "by R2 as well would count the error variance twice: remove the fitting columns ",
         "and give range and nugget of the prediction errors.")
  legend_par$dom_code <- legend_par$soil_unit
  if (any(duplicated(legend_par$dom_code)))
    stop("Legend: duplicated soil unit(s): ",
         paste(unique(legend_par$dom_code[duplicated(legend_par$dom_code)]), collapse = ", "))
  miss_dom <- setdiff(unique(grid_df$dom_code), legend_par$dom_code)
  if (length(miss_dom) > 0)
    stop("Legend: no row for the soil unit(s) of the frame: ", paste(miss_dom, collapse = ", "),
         "\n(soil_unit must match the frame exactly)")
  rng_vals <- as.matrix(legend_par[, rng_cols])
  nug_vals <- as.matrix(legend_par[, nug_cols])
  if (any(!is.finite(rng_vals)) || any(rng_vals < 0))
    stop("Legend: every range value must be >= 0 (km)")
  if (any(!is.finite(nug_vals)) || any(nug_vals < 0 | nug_vals > 1))
    stop("Legend: every nugget value must be in [0, 1] (relative nugget c0/(c0+c1))")

  range_by_su_km <- as.data.frame(rng_vals, row.names = legend_par$dom_code)
  colnames(range_by_su_km) <- vars_cont_y
  nugget_by_su <- as.data.frame(nug_vals, row.names = legend_par$dom_code)
  colnames(nugget_by_su) <- vars_cont_y

  cat("\nSpatial parameters read from the legend (range in km, relative nugget):\n")
  print(cbind(setNames(range_by_su_km, paste0("range_", vars_cont_y, "_km")),
              setNames(nugget_by_su, paste0("nugget_", vars_cont_y))))
} else {
  message("spatial_component = FALSE: range = 0 in every domain; the legend is not read.")
}


# ---------------------------------------------------------------------------
# Phase 3 -- SamplingStrata frame construction. Y and variance come from the
# frame itself (predicted target variables and their variance); range/nugget come from
# the legend (Phase 2) if spatial_component = TRUE, else range = 0.
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
# df2$dom_code is already the soil_unit (built in grid_df above).
dom_code_orig <- as.character(df2$dom_code)   # preserve original soil_unit strings before integer conversion
df2$dom_code  <- as.integer(factor(df2$dom_code))

# Check the auxiliary (X) variables for negative values. SamplingStrata
# expects non-negative auxiliary variables, so any variable containing
# negative values is shifted upward (v - min(v)) to make all values >= 0.
neg_x <- vars_cont_x[sapply(vars_cont_x, function(v) any(df2[[v]] < 0, na.rm = TRUE))]
if (length(neg_x) > 0) {
  for (v in neg_x) {
    shift <- min(df2[[v]], na.rm = TRUE)
    message(sprintf("Variable '%s' has negative values (min = %.4f): shifted by %.4f.",
                    v, shift, -shift))
    df2[[v]] <- df2[[v]] - shift
  }
}

# Build the SamplingStrata spatial frame object: X = auxiliary variables,
# Y = predicted target variables, variance = their prediction variance,
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
  dom_code_str = dom_code_orig,       # original string (bridge for the legend lookup)
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

# Assign each re-indexed domain the range and nugget of its OWN soil unit:
# from the legend (Phase 2) if spatial_component = TRUE, otherwise range = 0
# and nugget = 0. The key is dom_code_str (soil_unit).
# range_by_domain_km / nugget_by_domain: one row per re-indexed domain
# (rownames "1", "2", ...), one column per target variable; the GA receives
# the range in metres (range_by_domain_m), like the grid coordinates.
new_to_str <- dom_code_to_num$dom_code_str[match(domain_map$old, dom_code_to_num$domainvalue)]

ny <- length(vars_cont_y)
if (spatial_component) {
  range_by_domain_km <- range_by_su_km[new_to_str, , drop = FALSE]
  nugget_by_domain   <- nugget_by_su[new_to_str, , drop = FALSE]
} else {
  range_by_domain_km <- as.data.frame(matrix(0, nrow = length(new_to_str), ncol = ny,
                                             dimnames = list(NULL, vars_cont_y)))
  nugget_by_domain   <- as.data.frame(matrix(0, nrow = length(new_to_str), ncol = ny,
                                             dimnames = list(NULL, vars_cont_y)))
}
rownames(range_by_domain_km) <- as.character(seq_len(nrow(range_by_domain_km)))
rownames(nugget_by_domain)   <- as.character(seq_len(nrow(nugget_by_domain)))
range_by_domain_m <- range_by_domain_km * 1000


# Set up precision targets: CV = cv_target for every target variable in every domain.
# The cv data frame is the input constraints table required by optimStrata.
ndom <- length(unique(frame$domainvalue))
cv_cols <- setNames(
  replicate(length(vars_cont_y), rep(cv_target, ndom), simplify = FALSE),
  paste0("CV", seq_along(vars_cont_y))
)
cv <- as.data.frame(c(
  list(DOM = rep("DOM1", ndom)),
  cv_cols,
  list(domainvalue = c(1:ndom))
))

# ---------------------------------------------------------------------------
# Phase 4 -- K-means clustering to provide strata initial solutions
# Run domain by domain (lapply over domain IDs):
#   spatial_component = TRUE : KmeansSolutionSpatial, with the domain's own
#                              soil-unit range/nugget from the legend;
#   spatial_component = FALSE: KmeansSolution2 (non-spatial; it clusters on
#                              the Y values and evaluates each clustering
#                              with the non-spatial Bethel allocation).
# Only the resulting number of strata per domain is passed to the GA.
# The maximum number of clusters per domain is capped at max_clusters and limited
# to (number of distinct Y-value combinations - 1) to avoid empty strata.
# A tryCatch wrapper (and an explicit row-count check, since the K-means
# functions can return a degenerate 0-strata result without raising an
# error) handles domains where clustering fails, assigning all units to a
# single stratum as a fallback.
# ---------------------------------------------------------------------------

# Determine the maximum feasible number of clusters per domain:
# at most max_clusters, at least 1, and strictly less than the number of distinct
# combinations of target variable values in that domain
y_cols <- paste0("Y", seq_along(vars_cont_y))
distinct_per_dom <- tapply(seq_len(nrow(frame)), frame$domainvalue, function(idx)
  nrow(unique(frame[idx, y_cols])))
max_k <- pmin(pmax(distinct_per_dom - 1L, 1L), max_clusters)

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
    range_dv   <- as.numeric(range_by_domain_m[as.character(dv), ])   # metres, like the coordinates
    nugget_dv  <- as.numeric(nugget_by_domain[as.character(dv), ])
    res <- tryCatch(
      KmeansSolutionSpatial(frame       = f_sub,
                            range       = range_dv,
                            kappa       = 1,     # range = range parameter of the exponential model
                            nugget      = nugget_dv,
                            errors      = cv_sub,
                            nstrata     = NA,
                            minnumstrat = min_sample_per_stratum,
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
                      minnumstrat = min_sample_per_stratum,
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
minnumstr_val <- min_sample_per_stratum
units_per_dom <- tapply(frame$domainvalue, frame$domainvalue, length)
max_feasible  <- floor(units_per_dom[names(nstrat)] / minnumstr_val)
nstrat        <- pmin(nstrat, max_feasible)
nstrat        <- pmax(nstrat, 2L)   # at least 2 strata per domain

# ---------------------------------------------------------------------------
# Phase 5 -- Genetic algorithm (GA) optimisation of the stratification
#
# spatial_component = TRUE : the GA runs with the soil-unit range/nugget read
#                            from the legend (range_by_domain_m /
#                            nugget_by_domain), i.e. it searches the optimal
#                            stratification under the spatial variance.
# spatial_component = FALSE: the GA runs with range = 0 in every domain
#                            (non-spatial).
# ---------------------------------------------------------------------------
set.seed(1234)    # the seed "kmeans1d" uses stats::kmeans()
sugg <- prepareSuggestion(frame, nStrata = nstrat)
attr(sugg, "seeds")                          # description of the seeds
evalSuggestion(frame, errors = cv, sugg, which = 1)  # 1-D k-means seed
evalSuggestion(frame, errors = cv, sugg, which = 2)  # univariate seed

set.seed(1234)
args <- list(
  method = "spatial",
  errors=cv,
  framesamp=frame,
  iter = ga_iter,
  pops = ga_pops,
  minnumstr = min_sample_per_stratum,
  nStrata = nstrat,    # initial number of strata (from K-means)
  range = range_by_domain_m,     # legend values in metres, or 0 everywhere (Phase 3)
  nugget = nugget_by_domain,     # legend values, or 0 everywhere (Phase 3)
  kappa = 1,                     # range = range parameter of the exponential model
  suggestions = sugg,
  writeFiles = FALSE,
  showPlot = FALSE,
  parallel = TRUE)

solution <- do.call(optimStrata, args)


# GA results
framenew  <- solution$framenew
outstrata <- solution$aggr_strata
ga_label  <- if (spatial_component) "spatial, range/nugget from the legend" else "non-spatial, range = 0 everywhere"
cat(sprintf("\nGA result (%s): total sample size = %.1f, %d strata\n",
            ga_label, sum(outstrata$SOLUZ), nrow(outstrata)))


# Keep the GA outstrata: it is the final allocation and the reference of the
# report.
outstrata_raw <- outstrata

# ---------------------------------------------------------------------------
# Phase 5b -- (only if spatial_component = TRUE) REPORT ONLY: reduction of the
# sample size due to the spatial component.
#
# For each domain, its own GA stratification (framenew$STRATO) is re-allocated
# (Bethel) without the spatial component (range = 0), via the same functions
# the GA uses internally (buildStrataDFSpatial + bethel_cpp, internal functions
# of SamplingStrata accessed via :::). The difference from the GA allocation is
# the reduction due to the spatial component on the same strata. The GA result
# is not changed: with the additive variance (no division by R2) the spatial
# component can only reduce the variance of a stratum, so no selection between
# spatial and non-spatial solutions is needed.
# ---------------------------------------------------------------------------
domain_log <- list()   # one row per domain, for the summary report
nospatial_total <- NA_real_
if (spatial_component) {
  buildStrataDFSpatial <- SamplingStrata:::buildStrataDFSpatial
  bethel_cpp <- SamplingStrata:::bethel_cpp
  nospatial_total <- 0
  for (dv in sort(unique(framenew$DOMAINVALUE))) {
    sub <- framenew[framenew$DOMAINVALUE == dv, ]
    strato_levels <- sort(unique(sub$STRATO))
    sub$X1 <- match(sub$STRATO, strato_levels)  # re-index to 1..k for buildStrataDFSpatial
    sub$DOMAINVALUE <- 1L
    colnames(sub) <- toupper(colnames(sub))
    st_ns <- buildStrataDFSpatial(dataset = sub[, c("ID", "X1", "LON", "LAT",
                                                    paste0("Y", seq_len(ny)), paste0("VAR", seq_len(ny)),
                                                    "DOMAINVALUE")],
                                  range = rep(0, ny), kappa = 1,
                                  progress = FALSE, verbose = FALSE)
    st_ns$DOM1 <- 1L
    cv_dv <- cv[cv$domainvalue == dv, , drop = FALSE]
    cv_dv$domainvalue <- 1L
    n_ns <- tryCatch(sum(bethel_cpp(st_ns, cv_dv, minnumstr = min_sample_per_stratum,
                                    realAllocation = TRUE)),
                     error = function(e) NA_real_)
    nospatial_total <- nospatial_total + n_ns
    ga_reg <- outstrata_raw[outstrata_raw$DOM1 == dv, ]
    domain_log[[length(domain_log) + 1]] <- data.frame(
      domain = dv, dom_code = new_to_str[dv],
      n_strata_ga = nrow(ga_reg), n_ga = sum(ga_reg$SOLUZ, na.rm = TRUE),
      n_nospatial = n_ns, stringsAsFactors = FALSE)
  }
  cat(sprintf("\nTotal sample size, GA result (spatial)          : %.1f\n", sum(outstrata_raw$SOLUZ)))
  cat(sprintf("Total sample size, same strata without spatial component: %.1f\n", nospatial_total))
} else {
  for (dv in sort(unique(framenew$DOMAINVALUE))) {
    ga_reg <- outstrata_raw[outstrata_raw$DOM1 == dv, ]
    domain_log[[length(domain_log) + 1]] <- data.frame(
      domain = dv, dom_code = new_to_str[dv],
      n_strata_ga = nrow(ga_reg), n_ga = sum(ga_reg$SOLUZ, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
}

# ---------------------------------------------------------------------------
# Phase 6 -- Sample selection and output assembly
# ---------------------------------------------------------------------------
strataStructure <- summaryStrata(framenew, outstrata, progress = FALSE)

cv_out <- as.data.frame(expected_CV(outstrata))
colnames(cv_out) <- paste0("cv_", vars_cont_y)

# framenew_samp: link framenew to grid_df through the cell id
framenew$cell <- df2$cell[match(as.integer(framenew$ID), df2$id)]
framenew_samp <- merge(
  subset(framenew, cell %in% grid_df$cell & !duplicated(cell)),
  grid_df[, c("cell", "dom_code", "x", "y",
              if (use_lucas_grid) "LUCAS_GRID" else NULL)],
  by = "cell"
)

key          <- unique(framenew_samp[, c("DOMAINVALUE", "STRATO")])
outstrata_samp <- merge(outstrata, key,
                        by.x = c("DOM1", "STRATO"),
                        by.y = c("DOMAINVALUE", "STRATO"))

# Spatial sampling
set.seed(1234)
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

# Selected points (dom_code renamed soil_unit for the outputs)
points_all <- points1 %>% rename(soil_unit = dom_code)

message("Nrow points_all: ",      nrow(points_all))
message("Sampled lucas grid: ",   sum(points_all$source == "sampled_lucas"))
message("Sampled: ",              sum(points_all$source == "sampled"))

# Final tab: one row per domain (soil unit); dom_val flags domains with at
# least min_point_per_domain cells (those actually stratified)
tab2 <- grid_df %>%
  count(soil_unit, name = "n_point") %>%
  mutate(dom_val = n_point >= min_point_per_domain)

tab3 <- tab2 %>%
  rename(n_cell = n_point) %>%
  left_join(points_all %>% count(soil_unit, name = "n_tot"),                                        by = "soil_unit") %>%
  left_join(points_all %>% filter(source == "sampled_lucas") %>% count(soil_unit, name = "n_lucas_grid"), by = "soil_unit") %>%
  mutate(
    n_tot          = coalesce(n_tot, 0L),
    n_lucas_grid   = coalesce(n_lucas_grid, 0L)
  ) %>%
  dplyr::select(
    soil_unit,
    n_cell, dom_val,
    n_tot, n_lucas_grid
  )



#--------------------------------------------------------------
#Save files and figures
#--------------------------------------------------------------
#--------------------------------------------------------------
# Summary report (text file, English): main results of the run
#--------------------------------------------------------------
dom_tab <- do.call(rbind, domain_log)
dom_tab$soil_unit <- dom_tab$dom_code

rep_lines <- character(0)
add  <- function(...) rep_lines <<- c(rep_lines, paste0(...))
addt <- function(df) {
  op <- options(width = 250)   # keep wide tables on one line
  on.exit(options(op))
  rep_lines <<- c(rep_lines, capture.output(print(df, row.names = FALSE)))
}
hdr  <- function(txt) add("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78))

add("SOIL MONITORING LAW (SML) - OPTIMAL SAMPLING: RESULTS SUMMARY")
add("Run date            : ", format(Sys.time(), "%Y-%m-%d %H:%M"))
add("Scenario            : ", base_name)
add("Domains             : ", nrow(dom_tab), " soil units")
add("Target variables    : ", paste(vars_cont_y, collapse = ", "), " (CV target = ", 100 * cv_target, "% in every domain)")
add("Stratification vars : ", paste(vars_cont_x, collapse = ", "))
add("spatial_component   : ", spatial_component)

hdr("1. SPATIAL PARAMETERS (Phase 2)")
if (spatial_component) {
  add("Read from ", legend_csv_path, ". Range parameter of the exponential variogram of the")
  add("prediction errors (km, kappa = 1) and relative nugget c0/(c0+c1).")
  add("Used in KmeansSolutionSpatial (Phase 4) and optimStrata (Phase 5).\n")
  par_out <- data.frame(domain = seq_along(new_to_str),
                        soil_unit = new_to_str,
                        round(setNames(range_by_domain_km, paste0("range_", vars_cont_y, "_km")), 1),
                        round(setNames(nugget_by_domain, paste0("nugget_", vars_cont_y)), 3),
                        check.names = FALSE)
  addt(par_out)
} else {
  add("Not used: range = 0 in every domain (non-spatial). Initial solution by")
  add("KmeansSolution2 (Phase 4); the legend is not read.")
}

hdr(paste0("2. OPTIMISATION (Phase 5: GA, ", ga_label, ")"))
add("Sample size per domain:\n")
ga_out <- dom_tab[, c("domain", "soil_unit", "n_strata_ga", "n_ga")]
ga_out$n_ga <- round(ga_out$n_ga, 1)
names(ga_out) <- c("domain", "soil unit", "strata", "sample size")
addt(ga_out)
add("\nTOTAL (GA result) : ", round(sum(outstrata_raw$SOLUZ), 1), " units in ", nrow(outstrata_raw), " strata")

hdr("3. REDUCTION DUE TO THE SPATIAL COMPONENT (Phase 5b)")
if (!spatial_component) {
  add("Not performed (spatial_component = FALSE).")
} else {
  add("The strata found by the GA are re-allocated (Bethel) without the spatial component")
  add("(range = 0). The GA result is the final allocation.\n")
  sp_out <- dom_tab[, c("domain", "soil_unit", "n_ga", "n_nospatial")]
  sp_out$difference <- sp_out$n_nospatial - sp_out$n_ga
  sp_out[, 3:5] <- round(sp_out[, 3:5], 1)
  names(sp_out) <- c("domain", "soil unit", "n spatial (GA)", "n without spatial component",
                     "reduction")
  addt(sp_out)
  add("\nTOTAL, GA result (spatial)                      : ", round(sum(outstrata_raw$SOLUZ), 1))
  add("TOTAL, same strata without spatial component    : ", round(nospatial_total, 1))
  add("Reduction due to the spatial component          : ", round(nospatial_total - sum(outstrata_raw$SOLUZ), 1),
      " units (", round(100 * (nospatial_total - sum(outstrata_raw$SOLUZ)) / nospatial_total, 1), "%)")
}

hdr("4. FINAL SAMPLE")
add("Total sites selected : ", nrow(points_all))
add("  sampled           : ", sum(points_all$source == "sampled"))
add("  sampled (LUCAS grid): ", sum(points_all$source == "sampled_lucas"))
add("\nExpected CV per domain (target ", 100 * cv_target, "%), maximum over domains:")
for (cn in colnames(cv_out)) add("  ", cn, " : ", sprintf("%.4f", max(cv_out[[cn]], na.rm = TRUE)))
add("\nSites per domain (domains with at least ", min_point_per_domain, " cells are stratified):\n")
addt(as.data.frame(tab3))

hdr("5. OTHER NOTES")
add("* The CV of the target variables is computed on the model-predicted values and variances of the frame.")
add("\nOutput files: ", base_name, "_units.csv, _allpoints.csv, _sample.csv, _figures.pdf, .gpkg, .RData, _strata_SML.gpkg")
writeLines(rep_lines, paste0(out_path, base_name, "_summary.txt"))

write.csv(tab3, paste0(out_path, paste0(base_name, "_units.csv")), row.names=FALSE)
write.csv(points_all, paste0(out_path, paste0(base_name, "_allpoints.csv")), row.names=FALSE)
write.table(sample1,file.path(out_path, paste0(base_name, "_sample.csv")), sep=",",row.names=F,quote=F)

# Save domains as geopackage
grid_sf <- st_as_sf(dplyr::select(grid_df, -dom_code), coords = c("x", "y"), crs = 3035) %>%
  left_join(tab3[, c("soil_unit", "dom_val", "n_tot", "n_lucas_grid")], by = "soil_unit")

st_write(grid_sf, paste0(out_path, paste0(base_name, ".gpkg")),
         delete_dsn = TRUE)


#################################################################################
# Figures
#################################################################################
# Correlation matrix of target and auxiliary variables (frame cells)
cor_vars <- c(vars_cont_y, vars_cont_x)
cor_mat  <- cor(grid_df[, cor_vars], use = "pairwise.complete.obs")
message("Correlation matrix (", paste(cor_vars, collapse = ", "), "):")
print(round(cor_mat, 3))

plot_cor_matrix <- function(m) {
  n   <- ncol(m)
  pal <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(201)
  op  <- par(mar = c(1, 6, 6, 1))
  on.exit(par(op))
  image(1:n, 1:n, t(m[n:1, ]), col = pal, zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "", main = "Correlation matrix - frame")
  axis(3, at = 1:n, labels = colnames(m), tick = FALSE, las = 2)
  axis(2, at = 1:n, labels = rev(rownames(m)), tick = FALSE, las = 1)
  txt_cex <- min(1, 8 / n)   # smaller labels when there are many variables
  for (i in 1:n) for (j in 1:n) {
    text(j, n - i + 1, sprintf("%.2f", m[i, j]), cex = txt_cex,
         col = ifelse(abs(m[i, j]) > 0.6, "white", "black"))
  }
}

pdf(paste0(out_path, paste0(base_name, "_figures.pdf")), width = 8.5, height = 8)

# Correlation matrix
plot_cor_matrix(cor_mat)

# CV boxplot
boxplot(cv_out,
        col  = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
        main = paste0("Coefficient of Variation Distribution (sample = ", nrow(points_all), ")"),
        ylab = "CV",
        xlab = "")
abline(h = cv_target, col = "red")

# Convert to sf object
s_sf1 <- st_as_sf(points_all, coords = c("x", "y"), crs = 3035) %>%
  st_transform(crs = 4326) %>%
  mutate(soil_unit = as.factor(soil_unit))

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
    n_sampled_only    = n_tot - n_lucas_grid,
    n_lucas_grid_only = n_lucas_grid
  ) %>%
  dplyr::select(soil_unit, n_sampled_only, n_lucas_grid_only) %>%
  tidyr::pivot_longer(cols      = c(n_sampled_only, n_lucas_grid_only),
                      names_to  = "type",
                      values_to = "count") %>%
  mutate(type = factor(type,
                       levels = c("n_sampled_only", "n_lucas_grid_only"),
                       labels = c("Sampled", "LUCAS Grid")))
print(
  ggplot(tab3_long, aes(x = soil_unit, y = count, fill = type)) +
    geom_col() +
    scale_fill_manual(values = c("Sampled"      = "steelblue",
                                 "LUCAS Grid"   = "darkorange")) +
    labs(title = paste0("Number of Points per Domain (total = ", sum(tab3$n_tot), ")"),
         x     = "Soil unit",
         y     = "Count",
         fill  = "Point Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

# Geographic distribution
pts_all_plot <- points_all %>%
  mutate(type = case_when(
    LUCAS_GRID > 0    ~ "LUCAS Grid",
    TRUE              ~ "Sampled"
  )) %>%
  dplyr::select(x, y, type)

pts_all_sf <- st_as_sf(pts_all_plot, coords = c("x", "y"), crs = 3035) %>%
  st_transform(crs = 4326)

print(
  ggplot(pts_all_sf) +
    geom_sf(aes(color = type, size = type, alpha = type)) +
    scale_color_manual(values = c("Sampled" = "steelblue", "LUCAS Grid" = "darkorange")) +
    scale_size_manual(values  = c("Sampled" = 0.5,         "LUCAS Grid" = 0.8)) +
    scale_alpha_manual(values = c("Sampled" = 0.4,         "LUCAS Grid" = 0.6)) +
    labs(title = paste0("Geographic Distribution of Points (total = ", nrow(pts_all_plot), ")"),
         color = "Point Type", size = "Point Type", alpha = "Point Type") +
    theme_minimal()
)

dev.off()


# ---------------------------------------------------------------------------
# Export optimized strata as polygons (to be copied to Files/Strata/ of the SML tool)
# ---------------------------------------------------------------------------
res_m <- 1000   # grid resolution (m)

strata_polygons <- framenew_samp %>%
  st_as_sf(coords = c("x", "y"), crs = 3035) %>%
  st_buffer(dist = res_m / 2, endCapStyle = "SQUARE", joinStyle = "MITRE") %>%
  group_by(soil_unit = dom_code, STRATO) %>%
  summarise(.groups = "drop") %>%
  mutate(stratum_id = paste(soil_unit, STRATO, sep = "*"))
st_write(strata_polygons,
         file.path(out_path, paste0(base_name, "_strata_SML.gpkg")),
         delete_dsn = TRUE)

# Full workspace
save.image(file = paste0(out_path, base_name, ".RData"))
