# =============================================================================
# Kriging-based spatial parameters for the SML sampling design
# EXAMPLE: 6 soil units in 2 zones (Toscana meridionale; Emilia orientale -
# delta del Po), on a 1 km grid
#
# Authors : Claudia Cagnarini, Giulio Barcaroli
# Adapted : 10/11/2026, split from SML_sampling_example_with_6_Italy_soil_units.R
#           (variogram fitting, formerly its Phase 3)
# =============================================================================
#
# PURPOSE
# -------
# Estimate, for every soil unit of the test area, the two spatial parameters
# used by SamplingStrata (method = "spatial") for each target variable
# (soil pH = PH, soil organic carbon = SOC):
#   - fitting : R2 of the variogram model fit
#   - range   : variogram range, in km
# and write them to legend.csv, which is read by the sampling script
# SML_sampling_example_with_6_Italy_soil_units.R when its flag
# spatial_component = TRUE. The script also assesses how reliable these
# parameters are and warns against using them when the fit is poor or the
# range is large compared with the size of the soil units.
#
# This script is OPTIONAL and independent of the sampling script: run it only
# when you want to (re)estimate the spatial parameters. The sampling script
# never runs it; it only reads legend.csv. legend.csv can also be edited by
# hand (e.g. fitting = 1 and range = 0 for a soil unit makes that domain
# non-spatial).
#
# INPUT (folder Input/)
# ---------------------
#   Test_6_soil_units_Italy.gpkg     legend: one polygon per soil unit, with
#                                    fields "zona" and "soil_unit"
#                                    (EPSG:3035). It gives the list of soil
#                                    units written to legend.csv, the zone of
#                                    each survey point (spatial join) and the
#                                    size of each soil unit.
#   SharewithMS20260825-IT.csv       point-level soil survey (X_LAEA, Y_LAEA,
#                                    pH_H2O, SOC), on which the variograms are
#                                    fitted
#   DEM_3035_1000m.tif               elevation raster, 1 km, EPSG:3035; its
#                                    value at each survey point is the DEM
#                                    covariate of the variogram trend
#
# PROCESSING PHASES
# -----------------------------------------------------------------------------
# Phase K1 - Set-up: parameters, legend and survey reading, DEM extraction at
#            the survey points, zone of each survey point (spatial join).
# Phase K2 - Variogram fitting. For PH and SOC a regression-kriging variogram
#            (trend: x + y + DEM, spherical model) is fitted on the raw survey
#            points: nationally, and separately for each zone. A zone with
#            fewer than min_obs_district valid points, or with a degenerate
#            fit (R2 <= 0), falls back to the national model. A fitted range
#            beyond range_frac_max of the largest empirical lag is capped
#            (and flagged). No kriging prediction is made: only the variogram
#            parameters are needed.
# Phase K3 - Parameters per soil unit. Each soil unit takes the fitting and
#            range of its zone. Its size is measured as the mean distance
#            between two random points of the soil unit, and the implied
#            correlation of the prediction errors at that distance,
#            exp(-kappa * d / range), is computed.
# Phase K4 - Reliability assessment and outputs: legend.csv, variogram plots,
#            text report with fit results, reliability class per soil unit
#            and variable, warnings and a recommendation on the flag
#            spatial_component of the sampling script.
#
# RELIABILITY RULES (thresholds set in the parameters block)
# ----------------------------------------------------------
#   NOT RELIABLE     fitting < fitting_low_threshold
#   USE WITH CAUTION fitting < fitting_warn_threshold, or the zone uses the
#                    national model (not confirmed by local data), or the
#                    range was capped at the data support, or the correlation
#                    at the mean within-unit distance exceeds rho_warn_threshold
#                    (large range compared with the soil unit)
#   ADEQUATE         none of the above
#
# OUTPUT
# ------
#   Input/legend.csv                  one row per soil unit: zona, soil_unit,
#     (path: legend_csv_path)         fitting_<var>, range_<var>_km for each
#                                     target variable (the only columns read
#                                     by the sampling script), plus
#                                     diagnostics (source and number of
#                                     points of the model, range capped,
#                                     mean distance, correlation, reliability)
#   Output/<base_nome>_variograms.pdf empirical variograms and fitted models
#                                     (national and zone fits)
#   Output/<base_nome>_kriging_report.txt
#                                     text report: settings, fit results,
#                                     parameters per soil unit, reliability,
#                                     warnings, recommendation
#
# =============================================================================

# -----------------------------------------------------------------------------
# Install the required packages, if not yet installed
# -----------------------------------------------------------------------------
cran_pkgs <- c("sf", "terra", "gstat")
to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(sf)
library(terra)
library(gstat)

# Modify the following paths:
inp_path <- "D:\\Google Drive\\ISPRA\\SML_sampling_procedure_example\\Input\\"
out_path <- "D:\\Google Drive\\ISPRA\\SML_sampling_procedure_example\\Output\\"
kriging_survey_path <- paste0(inp_path, "SharewithMS20260825-IT.csv")
# legend.csv is an INPUT of the sampling script: it is written in the input
# folder. The sampling script must read it from the same path.
legend_csv_path <- paste0(inp_path, "legend.csv")

#----------------------------------------------------------------------------------------------------
# Parameters & settings (to be modified by users)
#----------------------------------------------------------------------------------------------------
base_nome <- "Example_with_6_soil_units"   # prefix of the plot and report files

vars_cont_y <- c("PH", "SOC")   # target variables (names used in the frame and in legend.csv)
# vars_cont_y -> column of the point survey holding the same variable
kriging_source_col <- c(PH = "pH_H2O", SOC = "SOC")

# Trend covariates added to x + y in the regression-kriging variogram
# (columns of the survey after the DEM extraction).
trend_covars <- list(
  PH  = c("DEM"),
  SOC = c("DEM")
)

min_obs_district <- 30    # below this, a zone falls back to the national model
range_frac_max   <- 0.8   # cap on the fitted range, as a fraction of the largest empirical lag

# Reliability thresholds (see the header)
fitting_low_threshold  <- 0.5   # fitting below this: NOT RELIABLE
fitting_warn_threshold <- 0.9   # fitting below this: USE WITH CAUTION
rho_warn_threshold     <- 0.5   # correlation exp(-kappa*d/range) at the mean
                                # within-unit distance above this: USE WITH CAUTION
kappa <- 1                      # same value used in the sampling script
n_pts_extent <- 500             # random points per soil unit to measure its size

set.seed(1234)

# ---------------------------------------------------------------------------
# Phase K1 -- Set-up: legend, survey, DEM covariate, zone of each point
# ---------------------------------------------------------------------------
legend   <- st_read(paste0(inp_path, "Test_6_soil_units_Italy.gpkg"), quiet = TRUE)
dem_rast <- rast(paste0(inp_path, "DEM_3035_1000m.tif"))

soil_survey <- read.csv(kriging_survey_path)
survey_pts  <- vect(soil_survey, geom = c("X_LAEA", "Y_LAEA"), crs = "EPSG:3035")
survey_covs <- terra::extract(dem_rast, survey_pts, ID = FALSE)
names(survey_covs) <- c("DEM")
soil_survey <- cbind(soil_survey, survey_covs)
dat_survey  <- st_as_sf(soil_survey, coords = c("X_LAEA", "Y_LAEA"), crs = st_crs(legend), remove = FALSE)
# assign each survey point the "zona" of the legend polygon it falls in
dat_survey  <- st_join(dat_survey, legend[, "zona"], left = TRUE)

# ---------------------------------------------------------------------------
# Phase K2 -- Regression-kriging variogram fitting on POINT survey data,
# nationally and per zone.
#
# The variogram is fitted on the RAW point-level soil survey observations,
# not on the already-smoothed predicted values of the frame, which would show
# artificially strong spatial autocorrelation. For each target variable a
# trend v ~ x + y + covariates is removed and a spherical variogram is fitted
# on the residuals. The fitting returned is the R2 of the variogram model on
# the empirical semivariances (1 - SSE/SST); the R2 of the trend regression
# is also computed, for the report only.
# ---------------------------------------------------------------------------
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
  # the data support (the conservative "at least this far" reading) and flag it.
  range_cap_km <- range_frac_max * max(emp_v$dist) / 1000
  range_capped <- range_km > range_cap_km
  if (range_capped) {
    fit$range[2] <- range_cap_km * 1000
    range_km <- range_cap_km
  }

  # R2 of the trend regression (report only)
  trend_R2 <- summary(lm(form, data = st_drop_geometry(sub)))$r.squared

  list(sub = sub, formula = form, emp_v = emp_v, fit = fit, R2 = R2,
       range_km = range_km, range_capped = range_capped,
       trend_R2 = trend_R2, n = nrow(sub))
}

pdf(file.path(out_path, paste0(base_nome, "_variograms.pdf")), width = 8, height = 5)
national_vario <- list()
district_vario <- list()
zona_codes <- sort(unique(legend$zona[!is.na(legend$zona)]))
fit_log <- list()   # one row per scope x variable, for the report

for (v in vars_cont_y) {
  src_col <- kriging_source_col[[v]]
  covars  <- trend_covars[[v]]

  nat <- fit_variogram_model(src_col, dat_survey, covars)
  if (is.null(nat)) {
    stop(sprintf(paste(
      "fit_variogram_model() rejected the NATIONAL fit for %s (survey column '%s'):",
      "too few points, no convergence or R2 <= 0 on the full national sample. There is",
      "no further fallback below national scope -- inspect this variable's",
      "distribution/trend before retrying."
    ), v, src_col))
  }
  national_vario[[v]] <- nat
  fit_log[[length(fit_log) + 1]] <- data.frame(
    scope = "National", variable = v, R2 = nat$R2, trend_R2 = nat$trend_R2,
    range_km = nat$range_km, range_capped = nat$range_capped,
    n_used = nat$n, n_zone_valid = NA_integer_, source = "own fit", stringsAsFactors = FALSE)
  cat(sprintf("NATIONAL | %s (from '%s', covars: %s) -> R2 = %.3f, range = %.1f km, n = %d\n",
              v, src_col, paste(covars, collapse = ","), nat$R2, nat$range_km, nat$n))
  print(plot(nat$emp_v, nat$fit,
             main = sprintf("%s [survey column '%s'] NATIONAL\nrange=%.0f km, R2=%.3f",
                            v, src_col, nat$range_km, nat$R2)))

  district_vario[[v]] <- list()
  for (nc in zona_codes) {
    sub_nc  <- dat_survey[!is.na(dat_survey$zona) & dat_survey$zona == nc, ]
    ok_nc   <- !is.na(sub_nc[[src_col]])
    for (cc in covars) ok_nc <- ok_nc & !is.na(sub_nc[[cc]])
    n_valid <- sum(ok_nc)
    fit_d   <- if (n_valid >= min_obs_district) fit_variogram_model(src_col, sub_nc, covars) else NULL
    zone_fallback <- is.null(fit_d)
    if (zone_fallback) {
      reason <- if (n_valid < min_obs_district)
        sprintf("n = %d < %d valid observations", n_valid, min_obs_district)
      else
        sprintf("n = %d but the zone fit was degenerate (no convergence or R2 <= 0)", n_valid)
      message(sprintf("Zone %s: falling back to the national model for %s (%s)", nc, v, reason))
      fit_d <- nat
    } else {
      cat(sprintf("ZONE %s | %s -> R2 = %.3f, range = %.1f km, n = %d\n",
                  nc, v, fit_d$R2, fit_d$range_km, fit_d$n))
      print(plot(fit_d$emp_v, fit_d$fit,
                 main = sprintf("%s [survey column '%s'] ZONE %s\nrange=%.0f km, R2=%.3f",
                                v, src_col, nc, fit_d$range_km, fit_d$R2)))
    }
    district_vario[[v]][[as.character(nc)]] <- fit_d
    fit_log[[length(fit_log) + 1]] <- data.frame(
      scope = as.character(nc), variable = v, R2 = fit_d$R2, trend_R2 = fit_d$trend_R2,
      range_km = fit_d$range_km, range_capped = fit_d$range_capped,
      n_used = fit_d$n, n_zone_valid = n_valid,
      source = if (zone_fallback) "national model (fallback)" else "own fit",
      stringsAsFactors = FALSE)
  }
}
dev.off()
fit_tab <- do.call(rbind, fit_log)

# ---------------------------------------------------------------------------
# Phase K3 -- Parameters per soil unit, and size of each soil unit.
#
# Each soil unit takes the fit of its zone. Its size is measured as the mean
# distance d between two random points of the unit (n_pts_extent points drawn
# in the polygon). rho = exp(-kappa * d / range) is the correlation that the
# spatial variance of SamplingStrata assigns to the prediction errors of two
# cells at that distance. Strata are subsets of the soil unit, so distances
# inside a stratum are shorter and the correlation is even higher: rho is a
# lower bound for the strata.
# ---------------------------------------------------------------------------
su_key <- unique(st_drop_geometry(legend)[, c("zona", "soil_unit")])
su_key <- su_key[order(su_key$zona, su_key$soil_unit), ]
su_geom <- do.call(c, lapply(seq_len(nrow(su_key)), function(i)
  st_union(st_geometry(legend)[legend$zona == su_key$zona[i] &
                                 legend$soil_unit == su_key$soil_unit[i]])))
su <- st_sf(su_key, geometry = su_geom)

mean_dist_km <- sapply(seq_len(nrow(su)), function(i) {
  pts <- st_sample(st_geometry(su)[i], size = n_pts_extent, type = "random")
  d   <- st_distance(pts)
  mean(as.numeric(d[upper.tri(d)])) / 1000
})

leg_out <- data.frame(zona = su$zona, soil_unit = su$soil_unit, stringsAsFactors = FALSE)
for (v in vars_cont_y) {
  fz <- lapply(as.character(leg_out$zona), function(z) district_vario[[v]][[z]])
  ft <- fit_tab[fit_tab$variable == v, ]
  src <- ft$source[match(as.character(leg_out$zona), ft$scope)]
  leg_out[[paste0("fitting_", v)]]  <- round(sapply(fz, `[[`, "R2"), 4)
  leg_out[[paste0("range_", v, "_km")]] <- round(sapply(fz, `[[`, "range_km"), 2)
  leg_out[[paste0("source_", v)]]   <- ifelse(src == "own fit", "zone fit", "national fallback")
  leg_out[[paste0("n_", v)]]        <- sapply(fz, `[[`, "n")
  leg_out[[paste0("range_capped_", v)]] <- sapply(fz, `[[`, "range_capped")
}
leg_out$mean_dist_km <- round(mean_dist_km, 1)

# ---------------------------------------------------------------------------
# Phase K4 -- Reliability assessment, legend.csv and report
# ---------------------------------------------------------------------------
level_order <- c("ADEQUATE", "USE WITH CAUTION", "NOT RELIABLE")
issues_log  <- list()
for (v in vars_cont_y) {
  fit_v  <- leg_out[[paste0("fitting_", v)]]
  rng_v  <- leg_out[[paste0("range_", v, "_km")]]
  rho_v  <- exp(-kappa * leg_out$mean_dist_km / rng_v)
  leg_out[[paste0("rho_", v)]] <- round(rho_v, 3)
  rel <- character(nrow(leg_out))
  for (i in seq_len(nrow(leg_out))) {
    iss <- character(0)
    if (fit_v[i] < fitting_low_threshold)
      iss <- c(iss, sprintf("fitting %.2f < %.2f", fit_v[i], fitting_low_threshold))
    else if (fit_v[i] < fitting_warn_threshold)
      iss <- c(iss, sprintf("fitting %.2f < %.2f", fit_v[i], fitting_warn_threshold))
    if (leg_out[[paste0("source_", v)]][i] != "zone fit") iss <- c(iss, "national model")
    if (leg_out[[paste0("range_capped_", v)]][i]) iss <- c(iss, "range capped at data support")
    if (rho_v[i] > rho_warn_threshold)
      iss <- c(iss, sprintf("range %.0f km vs mean distance %.0f km (rho = %.2f)",
                            rng_v[i], leg_out$mean_dist_km[i], rho_v[i]))
    rel[i] <- if (fit_v[i] < fitting_low_threshold) "NOT RELIABLE"
              else if (length(iss) > 0) "USE WITH CAUTION" else "ADEQUATE"
    issues_log[[length(issues_log) + 1]] <- data.frame(
      zona = leg_out$zona[i], soil_unit = leg_out$soil_unit[i], variable = v,
      reliability = rel[i], issues = if (length(iss)) paste(iss, collapse = "; ") else "-",
      stringsAsFactors = FALSE)
  }
  leg_out[[paste0("reliability_", v)]] <- rel
}
rel_cols <- paste0("reliability_", vars_cont_y)
leg_out$reliability <- apply(leg_out[, rel_cols, drop = FALSE], 1, function(r)
  level_order[max(match(r, level_order))])
issues_tab <- do.call(rbind, issues_log)

write.csv(leg_out, legend_csv_path, row.names = FALSE)
message("legend.csv written to: ", legend_csv_path)

# --- Text report --------------------------------------------------------------
rep_lines <- character(0)
add  <- function(...) rep_lines <<- c(rep_lines, paste0(...))
addt <- function(df) {
  op <- options(width = 250)   # keep wide tables on one line
  on.exit(options(op))
  rep_lines <<- c(rep_lines, capture.output(print(df, row.names = FALSE)))
}
hdr  <- function(txt) add("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78))

add("SOIL MONITORING LEVEL (SML) - KRIGING-BASED SPATIAL PARAMETERS: REPORT")
add("Run date            : ", format(Sys.time(), "%Y-%m-%d %H:%M"))
add("Scenario            : ", base_nome)
add("Survey file         : ", basename(kriging_survey_path), " (", nrow(soil_survey), " points)")
add("Target variables    : ", paste(sprintf("%s (survey column '%s')", vars_cont_y,
                                            kriging_source_col[vars_cont_y]), collapse = ", "))
add("Variogram model     : spherical, on the residuals of the trend x + y + ",
    paste(unique(unlist(trend_covars[vars_cont_y])), collapse = " + "))
add("Soil units          : ", nrow(leg_out), " in ", length(zona_codes), " zones")
add("Output legend file  : ", legend_csv_path)

hdr("1. VARIOGRAM FITS (national and per zone)")
add("fitting = R2 of the spherical variogram on the empirical semivariances (the value")
add("written to legend.csv); trend R2 = R2 of the regression x + y + covariates (report")
add("only); range in km. A zone with fewer than ", min_obs_district, " valid survey points, or with a")
add("degenerate fit, uses the national model. A range beyond ", range_frac_max,
    " x the largest empirical lag is capped.\n")
fit_out <- fit_tab
fit_out$R2       <- round(fit_out$R2, 3)
fit_out$trend_R2 <- round(fit_out$trend_R2, 3)
fit_out$range_km <- round(fit_out$range_km, 1)
names(fit_out) <- c("scope", "variable", "fitting", "trend R2", "range (km)", "range capped",
                    "n of model", "n valid in zone", "source")
addt(fit_out)

hdr("2. PARAMETERS PER SOIL UNIT (written to legend.csv)")
add("mean dist = mean distance between two random points of the soil unit;")
add("rho = exp(-kappa * mean dist / range), kappa = ", kappa,
    ": correlation of the prediction errors at that distance.\n")
p_cols <- c("zona", "soil_unit", "mean_dist_km",
            as.vector(t(outer(vars_cont_y, c("fitting_", "range_", "rho_"), function(v, p)
              ifelse(p == "range_", paste0(p, v, "_km"), paste0(p, v))))))
addt(leg_out[, p_cols])

hdr("3. RELIABILITY")
add("NOT RELIABLE     : fitting < ", fitting_low_threshold)
add("USE WITH CAUTION : fitting < ", fitting_warn_threshold, ", national model used for the zone,")
add("                   range capped at the data support, or rho > ", rho_warn_threshold)
add("ADEQUATE         : none of the above\n")
addt(issues_tab)
add("\nOverall, per soil unit (worst over the variables):")
addt(leg_out[, c("zona", "soil_unit", "reliability")])

hdr("4. WARNINGS")
warn <- character(0)
n_not <- sum(leg_out$reliability == "NOT RELIABLE")
n_cau <- sum(leg_out$reliability == "USE WITH CAUTION")
n_su <- nrow(leg_out)
# Describe a subset of issues_tab (one row per soil unit x variable) as
# "PH: 6 of 6 soil units; SOC: 5 of 6 soil units (all but BE_Cropland)".
by_var <- function(sub) {
  paste(sapply(vars_cont_y, function(v) {
    su_v <- sub$soil_unit[sub$variable == v]
    txt  <- sprintf("%s: %d of %d soil units", v, length(su_v), n_su)
    if (length(su_v) > 0 && length(su_v) < n_su)
      txt <- paste0(txt, " (", paste(su_v, collapse = ", "), ")")
    txt
  }), collapse = "; ")
}
low_fit <- issues_tab[grepl("^fitting", issues_tab$issues), ]
if (nrow(low_fit) > 0) {
  fvals <- unlist(leg_out[, paste0("fitting_", vars_cont_y)])
  warn <- c(warn, sprintf(paste0(
    "LOW FITTING: fitting below %.2f for %s (values %.2f-%.2f). In the within-stratum ",
    "variance of SamplingStrata the fitting divides the term of the predicted values, ",
    "(z_i - z_j)^2 / fitting; a low fitting inflates that term, i.e. it says that the frame ",
    "values explain little of the true variability."),
    fitting_warn_threshold, by_var(low_fit), min(fvals), max(fvals)))
}
big_rng <- issues_tab[grepl("rho =", issues_tab$issues), ]
if (nrow(big_rng) > 0) {
  rvals <- unlist(leg_out[, paste0("range_", vars_cont_y, "_km")])
  warn <- c(warn, sprintf(paste0(
    "LARGE RANGE: rho > %.2f for %s. The range (%.0f-%.0f km) is large compared with the ",
    "mean distance inside the soil units (%.0f-%.0f km; see section 2 for each soil unit). ",
    "In the within-stratum variance the prediction-error term is ",
    "(s_i^2 + s_j^2 - 2 s_i s_j exp(-kappa d_ij / range)): when the range is much larger than ",
    "the distances inside a stratum, the errors are treated as almost perfectly correlated ",
    "and this term nearly cancels. The sample size then decreases because of the range, ",
    "independently of the fitting: a large range not supported by local data gives a reduction ",
    "that is not warranted."),
    rho_warn_threshold, by_var(big_rng), min(rvals), max(rvals),
    min(leg_out$mean_dist_km), max(leg_out$mean_dist_km)))
}
fb <- fit_tab[fit_tab$source != "own fit", ]
if (nrow(fb) > 0) {
  fb_issues <- issues_tab[grepl("national model", issues_tab$issues), ]
  n_zones   <- length(unique(fit_tab$scope[fit_tab$scope != "National"]))
  warn <- c(warn, sprintf(paste0(
    "NATIONAL MODEL USED AT ZONE LEVEL: the variogram is fitted per zone (%d zones x %d ",
    "variables = %d fits); in %d of these fits the zone model could not be estimated (fewer ",
    "than %d valid survey points, or degenerate fit) and the national model was used ",
    "(zones: %s). Soil units concerned: %s. The national spatial structure is not confirmed ",
    "by local data."),
    n_zones, length(vars_cont_y), sum(fit_tab$scope != "National"), nrow(fb),
    min_obs_district, paste(unique(fb$scope), collapse = "; "), by_var(fb_issues)))
}
cap <- fit_tab[fit_tab$range_capped, ]
if (nrow(cap) > 0) {
  warn <- c(warn, sprintf(paste0(
    "RANGE CAPPED: in %d fit(s) (%s) the fitted range exceeded %.1f x the largest empirical lag ",
    "and was capped: the range is not identified by the data."),
    nrow(cap), paste(unique(paste(cap$scope, cap$variable)), collapse = "; "), range_frac_max))
}
if (length(warn) == 0) add("No warnings.") else for (w in warn) add("* ", w, "\n")

hdr("5. RECOMMENDATION FOR THE SAMPLING SCRIPT")
if (n_not > 0) {
  add("Spatial parameters NOT RELIABLE for ", n_not, " of ", nrow(leg_out), " soil units.")
  add("Recommended: spatial_component <- FALSE in the sampling script (non-spatial design,")
  add("fitting = 1, range = 0). If spatial_component = TRUE is used anyway, the reduction of")
  add("the sample size should not be relied upon; alternatively set fitting = 1 and range = 0")
  add("in legend.csv for the soil units concerned, to make those domains non-spatial.")
} else if (n_cau > 0) {
  add("Spatial parameters to be USED WITH CAUTION for ", n_cau, " of ", nrow(leg_out), " soil units.")
  add("spatial_component = TRUE may be used, but compare the result with spatial_component = FALSE")
  add("and read the warnings above before relying on the reduction of the sample size.")
} else {
  add("Spatial parameters ADEQUATE for all soil units: spatial_component = TRUE can be used.")
}

hdr("6. NOTES")
add("* Each soil unit takes the variogram of its zone: soil units of the same zone share fitting and range.")
add("* The fitting written to legend.csv is the R2 of the variogram fit. In SamplingStrata the fitting")
add("  divides the variance of the predicted values, i.e. it plays the role of the R2 of the model that")
add("  produced the frame values; the trend R2 is reported in section 1 for comparison.")
add("* legend.csv can be edited by hand before running the sampling script; only the columns")
add("  zona, soil_unit, fitting_<var> and range_<var>_km are read (range in km).")
add("\nOutput files: ", basename(legend_csv_path), ", ", base_nome, "_variograms.pdf, ",
    base_nome, "_kriging_report.txt")
writeLines(rep_lines, paste0(out_path, base_nome, "_kriging_report.txt"))
message("Report written to: ", paste0(out_path, base_nome, "_kriging_report.txt"))
