# =============================================================================
# SML sampling site selection — Script 2: Spatial analysis
#
# Handles Phase 3 of the pipeline: soil
# district-level variogram fitting (and, optionally, domain-level
# variograms compared against the district estimates), producing the R²
# ("fitting") and range parameters later consumed by KmeansSolutionSpatial /
#                                          optimStrata(method = "spatial").
#
# Input (folder: intermediate/, subfolder of this script's own folder):
#   01_data_management_data.RData (from 01_data_management.R)
# Output (folder: intermediate/):
#   02_spatial_analysis_nuts_variogram.csv
#   02_spatial_analysis_domain_median_range.csv (if use_domain_range)
#   02_spatial_analysis_variograms.pdf
#   02_spatial_analysis_nuts_vs_domain_median.csv (if use_domain_range)
# =============================================================================

library(sf)
library(terra)
library(gstat)
library(dplyr)
library(tidyr)

#----------------------------------------------------------------------------------------------------
# Load script 1 output
#----------------------------------------------------------------------------------------------------
intermediate_path <- "intermediate/"
load(paste0(intermediate_path, "01_data_management_data.RData"))

# ---------------------------------------------------------------------------
# Soil district-level variogram — R² (fitting) and range
#
# Variograms are fitted separately for each soil district and each target
# variable. This respects spatial stationarity better than a single global
# model while providing stable estimates via large sample sizes.
# A random sample of max_pts_nuts_vario points is drawn when a district
# is too large for direct computation.
#
# For each district × variable combination, three trend models are evaluated:
#   - stationary  (v ~ 1)
#   - linear      (v ~ x + y)
#   - quadratic   (v ~ x + y + x² + y² + xy)
# The simplest model is retained unless a more complex one improves R² by
# more than r2_improvement_threshold (default: 0.05). This avoids
# over-fitting while capturing genuine large-scale non-stationarity.
#
# Two safeguards handle unsatisfactory fits:
#   1. Candidates whose fitted range exceeds range_frac_max × cutoff are
#      excluded from the selection: a spherical range at/beyond the data
#      support is extrapolated and usually reflects drift fitted as if it
#      were spatial structure (the detrended candidates take over).
#   2. If the best fit still has R² < r2_min_acceptable (or no candidate
#      has a valid range), the variogram is refitted at shorter cutoffs
#      (refit_cutoffs_m); the valid-range fit with the highest R² wins.
#
# Both R² and range are stored; domains are then assigned the values of
# their soil district in place of the global/within-domain estimates.
# Districts with fewer than min_pts_nuts_vario pixels are skipped and
# receive the cross-district column mean as fallback.
# ---------------------------------------------------------------------------
min_pts_nuts_vario <- 500
max_pts_nuts_vario <- 5000
nuts_codes <- sort(unique(stack5$NUTS))

nuts_vario_results <- list()


r2_improvement_threshold <- 0.05
vario_cutoff_m   <- 300000           # default empirical variogram cutoff (m)
vario_width_m    <- 10000            # lag width (m)
r2_min_acceptable <- 0.6             # fits with R² below this are refitted
range_frac_max    <- 0.8             # fitted range must stay below this fraction
# of the maximum empirical lag distance: beyond
# it the spherical model is extrapolated and
# signals unresolved drift
refit_cutoffs_m   <- c(150000, 100000)  # shorter cutoffs tried when refitting
# vario_cutoff_km defined once here (also duplicated, kept in sync, in
# 03_optimization_sampling.R where it is used to cap range_by_domain).
vario_cutoff_km   <- vario_cutoff_m / 1000
# Gates the domain-level variogram comparison, disabled by default.
use_domain_range  <- FALSE           # TRUE to compute domain-level variograms
# and compare with district-level estimates
# Toggles the SOURCE of the per-domain ranges fed to optimStrata
# (range_by_domain), consumed by 03_optimization_sampling.R:
#   FALSE -> soil district-level range for every domain (original behaviour).
#   TRUE  -> per-district MEDIAN of the domain-level ranges, falling back per
#            variable to the district range and then to the cross-district mean
#            where no median exists. Domain ranges capture within-domain
#            correlation only and are typically far shorter, hence more
#            conservative for the allocation.
# Requires use_domain_range = TRUE (domain-level variograms must be computed);
# the dependency is enforced just below.
# NOTE: this flag is duplicated (kept in sync) in 03_optimization_sampling.R,
# which needs it to decide whether to read the domain-median CSV.
use_domain_median_range <- TRUE
if (use_domain_median_range && !use_domain_range) {
  message("use_domain_median_range = TRUE requires use_domain_range = TRUE; enabling the latter.")
  use_domain_range <- TRUE
}

select_best_trend <- function(v, nuts_code, dat_nuts, sub_nuts,
                              cutoff = vario_cutoff_m, width = vario_width_m) {
  formulas <- list(
    stationary = as.formula(paste(v, "~ 1")),
    linear     = as.formula(paste(v, "~ x + y")),
    quadratic  = as.formula(paste(v, "~ x + y + I(x^2) + I(y^2) + I(x*y)"))
  )
  results <- list()
  for (trend_name in names(formulas)) {
    emp_v <- tryCatch(
      variogram(formulas[[trend_name]], data = dat_nuts, cutoff = cutoff, width = width),
      error = function(e) NULL
    )
    if (is.null(emp_v) || nrow(emp_v) < 3) next
    vvar <- var(sub_nuts[[v]], na.rm = TRUE)
    m0   <- vgm(psill = vvar * 0.9, model = "Sph", range = cutoff / 2, nugget = vvar * 0.1)
    fit  <- tryCatch(fit.variogram(emp_v, model = m0), error = function(e) NULL)
    if (is.null(fit)) next
    pred <- variogramLine(fit, dist_vector = emp_v$dist)
    SSE  <- sum((emp_v$gamma - pred$gamma)^2)
    SST  <- sum((emp_v$gamma - mean(emp_v$gamma))^2)
    R2   <- max(0, 1 - SSE / SST)
    # range validated against the actual data support (max empirical lag)
    # instead of the nominal cutoff. For districts the two coincide, but
    # domains can be far smaller than the cutoff, and a range extrapolated
    # beyond their extent would previously pass as valid.
    results[[trend_name]] <- list(formula = formulas[[trend_name]],
                                  emp_v = emp_v, fit = fit, R2 = R2,
                                  range_ok = fit$range[2] <= range_frac_max * max(emp_v$dist))
  }
  if (length(results) == 0) return(NULL)
  # Restrict the choice to candidates whose fitted range lies within the data
  # support: a high-R² "stationary" fit with range at/over the cutoff is drift
  # masquerading as spatial structure, and a detrended model must be preferred.
  valid <- Filter(function(x) isTRUE(x$range_ok), results)
  pool  <- if (length(valid) > 0) valid else results
  ord   <- intersect(c("stationary", "linear", "quadratic"), names(pool))
  best      <- pool[[ord[1]]]
  best_name <- ord[1]
  for (candidate in ord[-1]) {
    if (pool[[candidate]]$R2 > best$R2 + r2_improvement_threshold) {
      best <- pool[[candidate]]
      best_name <- candidate
    }
  }
  best$trend_name <- best_name
  best$cutoff <- cutoff
  best$R2_all <- sapply(results, function(x) round(x$R2, 4))
  best
}

# Refit unsatisfactory variograms. The fit at the default cutoff is kept only
# if R² >= r2_min_acceptable and the range lies within the data support.
# Otherwise the variogram is refitted at progressively shorter cutoffs: long
# lags carry few point pairs and hole effects that can destroy the fit of an
# otherwise sound short-range structure. The valid-range candidate with the
# highest R² across all cutoffs is retained.
fit_best_variogram <- function(v, code, dat_sub, sub_df) {
  best <- select_best_trend(v, code, dat_sub, sub_df)
  if (is.null(best)) return(NULL)
  best$refitted <- FALSE
  if (best$R2 >= r2_min_acceptable && isTRUE(best$range_ok)) return(best)
  for (co in refit_cutoffs_m) {
    cand <- select_best_trend(v, code, dat_sub, sub_df, cutoff = co)
    if (is.null(cand)) next
    if (isTRUE(cand$range_ok) && (!isTRUE(best$range_ok) || cand$R2 > best$R2)) {
      cand$refitted <- TRUE
      best <- cand
    }
  }
  best
}


pdf(paste0(intermediate_path, "02_spatial_analysis_variograms.pdf"), width = 8, height = 5)

for (nuts in nuts_codes) {
  sub_full <- stack5[stack5$NUTS == nuts, ]
  if (nrow(sub_full) < min_pts_nuts_vario) {
    cat("NUTS1", nuts, "- skipped (n =", nrow(sub_full), ")\n")
    next
  }
  set.seed(42)
  sub_nuts <- if (nrow(sub_full) > max_pts_nuts_vario)
    sub_full[sample(nrow(sub_full), max_pts_nuts_vario), ]
  else sub_full

  dat_nuts <- st_as_sf(sub_nuts, coords = c("x","y"),
                       crs = map_crs, remove = FALSE)

  for (v in vars_cont_y) {
    best <- fit_best_variogram(v, nuts, dat_nuts, sub_nuts)
    if (is.null(best)) next

    cat(sprintf("  DISTRICT %d | %s → trend: %s  (R²: stat=%.3f, lin=%.3f, quad=%.3f)%s\n",
                nuts, v, best$trend_name,
                best$R2_all["stationary"],
                best$R2_all["linear"],
                best$R2_all["quadratic"],
                if (best$refitted) sprintf("  [refitted, cutoff=%.0f km]", best$cutoff / 1000) else ""))

    print(plot(best$emp_v, best$fit,
               main = sprintf("DISTRICT %d — %s  (range=%.0f km, R²=%.3f) [%s%s]",
                              nuts, v, best$fit$range[2] / 1000,
                              best$R2, best$trend_name,
                              if (best$refitted) sprintf(", refit %.0f km", best$cutoff / 1000) else "")))

    nuts_vario_results[[paste(nuts, v, sep = "_")]] <- data.frame(
      NUTS      = nuts,
      variable  = v,
      n_pts     = nrow(sub_full),
      range_km  = round(best$fit$range[2] / 1000, 1),
      R2        = round(best$R2, 4),
      trend     = best$trend_name,
      cutoff_km = round(best$cutoff / 1000),
      refitted  = best$refitted
    )
  }
}

dev.off()



# ---------------------------------------------------------------------------
# Domain-level variogram range (supplementary — disabled by default)
# Set use_domain_range <- TRUE to compute domain-level variograms and
# compare R² and range with district-level estimates.
# For each domain × variable, the best trend model is selected automatically
# via select_best_trend (same logic as district level) for consistency.
# Only domains with R² >= 0.05 and range <= vario_cutoff_km are retained
# in the comparison. Results are saved to nuts_vs_domain_median.csv.
# # ---------------------------------------------------------------------------
nuts_results_df <- do.call(rbind, nuts_vario_results)
nuts_results_df$NUTS_CLASS <- nut_legend$NUTS_CLASS[
  match(nuts_results_df$NUTS, nut_legend$NUTS)
]
cat("\nVariogrammi distretto:\n")
print(nuts_results_df[, c("NUTS", "NUTS_CLASS", "variable", "n_pts", "range_km", "R2", "trend", "cutoff_km", "refitted")])

# Soil district lookup tables keyed by NUTS code (character rownames)
nuts_fitting_tbl <- nuts_results_df %>%
  dplyr::select(NUTS, variable, R2) %>%
  tidyr::pivot_wider(names_from = variable, values_from = R2) %>%
  as.data.frame()
rownames(nuts_fitting_tbl) <- as.character(nuts_fitting_tbl$NUTS)
nuts_fitting_tbl$NUTS <- NULL

nuts_range_tbl <- nuts_results_df %>%
  dplyr::select(NUTS, variable, range_km) %>%
  tidyr::pivot_wider(names_from = variable, values_from = range_km) %>%
  as.data.frame()
rownames(nuts_range_tbl) <- as.character(nuts_range_tbl$NUTS)
nuts_range_tbl$NUTS <- NULL

# Fallback = column means across districts (for any district without an estimate)
fitting_fallback  <- colMeans(nuts_fitting_tbl, na.rm = TRUE)
range_fallback_km <- colMeans(nuts_range_tbl,   na.rm = TRUE)
cat("\nFitting fallback (R²):\n");  print(fitting_fallback)
cat("Range fallback (km):\n");      print(range_fallback_km)

cat("\n--- nuts_range_tbl (district-level) ---\n");   print(nuts_range_tbl)
cat("\n--- nuts_fitting_tbl (district-level) ---\n"); print(nuts_fitting_tbl)

# Save the district-level fitting (R²) and range table — the "fitting" and
# "range" parameters used by KmeansSolutionSpatial / optimStrata — for
# 03_optimization_sampling.R to consume. Long format (one row per
# NUTS x variable) so it round-trips cleanly through a CSV; script 3
# pivots it back to the wide nuts_fitting_tbl / nuts_range_tbl layout and
# recomputes fitting_fallback / range_fallback_km via colMeans, exactly as
# done above.
write.csv(nuts_results_df, paste0(intermediate_path, "02_spatial_analysis_nuts_variogram.csv"), row.names = FALSE)
message("Saved: ", paste0(intermediate_path, "02_spatial_analysis_nuts_variogram.csv"))

# ---------------------------------------------------------------------------
# Domain-level variogram — comparison with district-level estimates
# For each domain with >= min_pts_dom_vario pixels, an empirical variogram
# is fitted using the same best-trend logic as the district level.
# Median R² and range across domains are then compared with the
# district-level estimates in a summary table saved to
# nuts_vs_domain_median.csv.
# ---------------------------------------------------------------------------
min_pts_dom_vario <- 500
dom_vario_results <- list()
# Per-district median-range lookup (rownames = NUTS, columns = variables);
# populated inside the block below when domain-level variograms are computed,
# and consumed by range_by_domain (in 03_optimization_sampling.R) when
# use_domain_median_range = TRUE. NULL means "no domain medians available".
nuts_range_dom_median_tbl <- NULL

# Block gated by use_domain_range (parameters block above): disabled by
# default. When enabled, fits variograms for every domain with
# >= min_pts_dom_vario pixels.
if (use_domain_range) {

  for (dc in sort(unique(stack5$dom_code))) {
    sub_dom <- stack5[stack5$dom_code == dc, ]
    if (nrow(sub_dom) < min_pts_dom_vario) next
    set.seed(42)
    sub_dom_s <- if (nrow(sub_dom) > max_pts_nuts_vario)
      sub_dom[sample(nrow(sub_dom), max_pts_nuts_vario), ] else sub_dom
    dat_dom <- st_as_sf(sub_dom_s, coords = c("x", "y"),
                        crs = map_crs, remove = FALSE)
    for (v in vars_cont_y) {
      vvar <- var(sub_dom_s[[v]], na.rm = TRUE)
      if (is.na(vvar) || vvar == 0) next
      best <- fit_best_variogram(v, dc, dat_dom, sub_dom_s)
      if (is.null(best)) next
      dom_vario_results[[paste(dc, v, sep = "_")]] <- data.frame(
        dom_code = dc,
        variable = v,
        n_pts    = nrow(sub_dom),
        range_km = round(best$fit$range[2] / 1000, 1),
        R2       = round(best$R2, 4),
        trend    = best$trend_name,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(dom_vario_results) > 0) {
    dom_results_df <- do.call(rbind, dom_vario_results)
    dom_results_df <- dom_results_df[dom_results_df$R2       >= 0.05, ]
    dom_results_df <- dom_results_df[dom_results_df$range_km <= vario_cutoff_km, ]

    domain_median_by_nuts <- dom_results_df %>%
      mutate(NUTS = as.integer(sapply(strsplit(dom_code, "\\*"), `[`, 1))) %>%
      group_by(NUTS, variable) %>%
      summarise(
        range_km_dom_median = median(range_km, na.rm = TRUE),
        R2_dom_median       = median(R2,       na.rm = TRUE),
        n_domains           = n(),
        .groups = "drop"
      )

    # Reshape the per-district median ranges into a wide lookup table (same
    # layout as nuts_range_tbl) so range_by_domain can use it.
    nuts_range_dom_median_tbl <- domain_median_by_nuts %>%
      dplyr::select(NUTS, variable, range_km_dom_median) %>%
      tidyr::pivot_wider(names_from = variable, values_from = range_km_dom_median) %>%
      as.data.frame()
    rownames(nuts_range_dom_median_tbl) <- as.character(nuts_range_dom_median_tbl$NUTS)
    nuts_range_dom_median_tbl$NUTS <- NULL

    # Save the long-format per-district domain-median table for
    # 03_optimization_sampling.R to consume (pivoted back to wide there).
    write.csv(domain_median_by_nuts, paste0(intermediate_path, "02_spatial_analysis_domain_median_range.csv"), row.names = FALSE)
    message("Saved: ", paste0(intermediate_path, "02_spatial_analysis_domain_median_range.csv"))

    comparison_by_nuts <- nuts_results_df %>%
      dplyr::select(NUTS, variable, range_km, R2) %>%
      rename(range_district_km = range_km, R2_district = R2) %>%
      left_join(domain_median_by_nuts, by = c("NUTS", "variable")) %>%
      dplyr::select(NUTS, variable, n_domains,
                    R2_district, R2_dom_median,
                    range_district_km, range_km_dom_median)

    cat("\n--- District vs domain median comparison ---\n")
    print(comparison_by_nuts)
    write.csv(comparison_by_nuts,
              paste0(intermediate_path, "02_spatial_analysis_nuts_vs_domain_median.csv"),
              row.names = FALSE)

  } else {
    cat("No domain had enough pixels (>=", min_pts_dom_vario,
        ") for domain-level variogram.\n")
  }

} # end of if (use_domain_range)
