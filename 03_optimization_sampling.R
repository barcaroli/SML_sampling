# =============================================================================
# SML sampling site selection — Script 3: Optimization & sampling
#
# Handles Phases 4-7 of the pipeline:
# SamplingStrata frame construction, spatial K-means initialisation, GA
# optimisation (optimStrata, method = "spatial"), spatial sample selection,
# and output assembly (tab3 summary, all-points table, figures, gpkg exports).
#
# Input (folder: intermediate/, subfolder of this script's own folder):
#   01_data_management_data.RData (from 01_data_management.R)
#   02_spatial_analysis_nuts_variogram.csv
#   02_spatial_analysis_domain_median_range.csv (if produced)
#   (both from 02_spatial_analysis.R)
# Output (folder: output/, subfolder of this script's own folder):
#   03_optimization_sampling_frame.RData
#   03_optimization_sampling_sample.csv
#   03_optimization_sampling_units.csv       (tab3)
#   03_optimization_sampling_allpoints.csv   (points_all)
#   03_optimization_sampling_figures.pdf
#   03_optimization_sampling.gpkg            (domains, with tab3 join)
#   03_optimization_sampling_strata.gpkg     (strata polygons)
# =============================================================================

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SamplingStrata)
library(bethel)
library(clhs)
library(parallel)
library(doParallel)

#----------------------------------------------------------------------------------------------------
# Parameters that must match 02_spatial_analysis.R
#----------------------------------------------------------------------------------------------------
# NOTE: these two values must match the ones used in 02_spatial_analysis.R —
# they gate how range_by_domain is capped/sourced and are not passed through
# the .RData/.csv hand-off because they belong to script 2's settings block.
vario_cutoff_km <- 300          # must match vario_cutoff_m / 1000 in script 2
use_domain_median_range <- TRUE # must match script 2's setting

#----------------------------------------------------------------------------------------------------
# Load script 1 and script 2 outputs; set up the output folder
#----------------------------------------------------------------------------------------------------
intermediate_path <- "intermediate/"
output_path <- "output/"
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

load(paste0(intermediate_path, "01_data_management_data.RData"))

nuts_results_df <- read.csv(paste0(intermediate_path, "02_spatial_analysis_nuts_variogram.csv"), stringsAsFactors = FALSE)

# Rebuild the wide district-level lookup tables and fallbacks exactly as
# 02_spatial_analysis.R computed them in-memory, just sourced from the CSV.
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

fitting_fallback  <- colMeans(nuts_fitting_tbl, na.rm = TRUE)
range_fallback_km <- colMeans(nuts_range_tbl,   na.rm = TRUE)

dom_median_file <- paste0(intermediate_path, "02_spatial_analysis_domain_median_range.csv")
nuts_range_dom_median_tbl <- NULL
if (use_domain_median_range && file.exists(dom_median_file)) {
  domain_median_by_nuts <- read.csv(dom_median_file, stringsAsFactors = FALSE)
  nuts_range_dom_median_tbl <- domain_median_by_nuts %>%
    dplyr::select(NUTS, variable, range_km_dom_median) %>%
    tidyr::pivot_wider(names_from = variable, values_from = range_km_dom_median) %>%
    as.data.frame()
  rownames(nuts_range_dom_median_tbl) <- as.character(nuts_range_dom_median_tbl$NUTS)
  nuts_range_dom_median_tbl$NUTS <- NULL
}

#----------------------------------------------------------------------------------------------------
# Building frame
#----------------------------------------------------------------------------------------------------
# Rebuild the sf wrapper (map_crs comes from script 1, since the raster
# stack1 itself is not serialized across scripts).
dat <- st_as_sf(stack5, coords = c("x", "y"), crs = map_crs, remove = FALSE)

df2 <- st_drop_geometry(dat)
df2$id <- c(1:nrow(df2))
dom_code_orig <- as.character(df2$dom_code)   # preserve original strings before integer conversion
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

# Build the SamplingStrata spatial frame object, specifying auxiliary
# variables (X), target variables (Y), their prediction variances, and
# spatial coordinates. domainvalue defines the stratification domain.
frame <- buildFrameSpatial(df=df2,
                           id="id",
                           X=vars_cont_x,
                           Y=vars_cont_y,
                           variance=variance_cont_y,
                           lon="x",
                           lat="y",
                           domainvalue = "dom_code")

# Filter out domains with fewer than 5 points: these are too small to
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

# Assign soil district-level fitting (R²) and range to each re-indexed domain.
# The soil district code is the first component of dom_code_str (format NUTS*IPCC*PEDO).
# Domains whose soil district has no variogram estimate receive the cross-district mean.
new_to_str <- dom_code_to_num$dom_code_str[
  match(domain_map$old, dom_code_to_num$domainvalue)
]
new_to_nuts <- as.integer(sapply(strsplit(new_to_str, "\\*"), `[`, 1))

# per-domain range source selectable via use_domain_median_range.
# Per variable the value is chosen as:
#   domain-median (if enabled & available) -> district-level -> cross-district mean.
get_range_row <- function(n) {
  key <- as.character(n)
  base <- if (!is.na(n) && key %in% rownames(nuts_range_tbl))
    nuts_range_tbl[key, , drop = FALSE]
  else
    as.data.frame(as.list(range_fallback_km))
  # element-wise fallback to the cross-district mean for any missing district value
  for (v in colnames(base)) if (is.na(base[[v]])) base[[v]] <- range_fallback_km[[v]]
  # optional override with the per-district domain-median where it exists
  if (use_domain_median_range && !is.null(nuts_range_dom_median_tbl) &&
      !is.na(n) && key %in% rownames(nuts_range_dom_median_tbl)) {
    med <- nuts_range_dom_median_tbl[key, , drop = FALSE]
    for (v in colnames(base))
      if (v %in% colnames(med) && !is.na(med[[v]])) base[[v]] <- med[[v]]
  }
  base
}
cat(sprintf("\nrange_by_domain source: %s\n",
            if (use_domain_median_range) "domain-median (district fallback)" else "district-level"))
range_by_domain <- do.call(rbind, lapply(new_to_nuts, get_range_row))
rownames(range_by_domain) <- as.character(seq_len(nrow(range_by_domain)))
colnames(range_by_domain) <- paste0("range_", colnames(range_by_domain))

# Cap range at the variogram cutoff (300 km): ranges fitted beyond the observation
# window are extrapolated and unreliable. This can occur especially when a
# detrended model (linear or quadratic) flattens the residual variogram,
# pushing the fitted range outside the data support.
# Very large ranges → near-zero effective within-stratum variance → Inf/NaN fitness in GA.
range_by_domain[] <- lapply(range_by_domain, function(col) pmin(col, vario_cutoff_km))
cat("range_by_domain after cap (max per variable):\n")
print(sapply(range_by_domain, max))

# R² per dominio: ogni dominio riceve l'R² del proprio soil district;
# i domini senza stima di district ricevono la media cross-regionale come fallback.
fitting_by_domain <- do.call(rbind, lapply(new_to_nuts, function(n) {
  key <- as.character(n)
  if (!is.na(n) && key %in% rownames(nuts_fitting_tbl))
    nuts_fitting_tbl[key, , drop = FALSE]
  else
    as.data.frame(as.list(fitting_fallback))
}))
rownames(fitting_by_domain) <- as.character(seq_len(nrow(fitting_by_domain)))

# Guard: R²=0 (failed variogram fit) causes optimStrata to compute zero or NaN
# within-stratum variance, making all GA fitness values Inf/NaN → plot crash.
# Replace any 0 with the cross-district column mean (fitting_fallback) as a safe floor.
for (col in colnames(fitting_by_domain)) {
  bad <- fitting_by_domain[[col]] == 0
  if (any(bad)) {
    message(sprintf("fitting_by_domain: R²=0 for %d domain(s) in '%s' — replacing with fallback %.4f",
                    sum(bad), col, fitting_fallback[[col]]))
    fitting_by_domain[[col]][bad] <- fitting_fallback[[col]]
  }
}

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
# Run KmeansSolutionSpatial domain by domain (lapply over domain IDs).
# The maximum number of clusters per domain is capped at 5 and limited
# to (number of distinct Y-value combinations - 1) to avoid empty strata.
# A tryCatch wrapper handles domains where clustering fails (e.g. too few
# distinct points), assigning all units to a single stratum as a fallback.
# ---------------------------------------------------------------------------

# Determine the maximum feasible number of clusters per domain:
# at most 10, at least 1, and strictly less than the number of distinct
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

# Named vector: new integer domain -> NUTS1 code (for per-domain lookup)
nuts_by_dv <- setNames(new_to_nuts, as.character(domain_map$new))

# Run spatial K-means clustering independently for each domain.
# NUTS1-specific R² (fitting) and range are passed so that domains within
# the same macro-region share the same spatial autocorrelation parameters.
kmean_list <- lapply(sort(unique(frame$domainvalue)), function(dv) {
  f_sub <- frame[frame$domainvalue == dv, ]
  f_sub$domainvalue <- 1L         # reset to 1 for the single-domain call
  cv_sub <- cv[cv$domainvalue == dv, , drop = FALSE]
  cv_sub$domainvalue <- 1L
  k_max <- min(max_k[[as.character(dv)]], nrow(f_sub) - 1L)
  k_max <- max(k_max, 1L)
  key <- as.character(nuts_by_dv[as.character(dv)])
  fitting_dv <- if (!is.na(key) && key %in% rownames(nuts_fitting_tbl))
    unlist(nuts_fitting_tbl[key, ])
  else
    fitting_fallback
  range_dv <- if (!is.na(key) && key %in% rownames(nuts_range_tbl))
    unlist(nuts_range_tbl[key, ]) * 1000   # km -> metres for KmeansSolutionSpatial
  else
    range_fallback_km * 1000
  res <- tryCatch(
    KmeansSolutionSpatial(frame       = f_sub,
                          fitting     = fitting_dv,
                          range       = range_dv,
                          kappa       = 1,
                          errors      = cv_sub,
                          nstrata     = NA,
                          minnumstrat = 2,
                          maxclusters = k_max,
                          showPlot    = FALSE),
    # Fallback: if clustering fails, assign all units to stratum "1"
    error = function(e) {
      cat("  [bypass] dv =", dv, " nrow =", nrow(f_sub),
          " k_max =", k_max, " error:", conditionMessage(e), "\n")
      id_col <- if ("id" %in% names(f_sub)) f_sub$id else seq_len(nrow(f_sub))
      data.frame(id          = as.character(id_col),
                 suggestions = factor(rep("1", nrow(f_sub))),
                 domainvalue = dv,
                 stringsAsFactors = FALSE)
    }
  )
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

nstrat
# ---------------------------------------------------------------------------
# Genetic algorithm (GA) optimisation of the stratification
# optimStrata searches for the allocation that minimises total sample size
# subject to the CV constraints defined in cv, using the spatial method
# which accounts for spatial autocorrelation via variogram R² and range.
# The K-means solution (nstrat) provides the initial number of strata.
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
  fitting = fitting_by_domain,  # data frame (ndom × n_vars): R² NUTS1-specifico per ogni dominio
  range = range_by_domain,
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

strataStructure <- summaryStrata(framenew, outstrata, progress = TRUE)
sum(strataStructure$Allocation)

cv_out <- as.data.frame(expected_CV(outstrata))
colnames(cv_out) <- paste0("cv_", vars_cont_y)
boxplot(cv_out)

stack_input <- if (use_lucas_points) stack5_samp else stack5

# framenew_samp: collega framenew a stack_input tramite cell
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

#--------------------------------------------------------------------------------------------------
# Save outputs: optimized frame (.RData) and selected sample (.csv)
#--------------------------------------------------------------------------------------------------
save(
  list = c("solution", "framenew", "outstrata", "framenew_samp", "outstrata_samp", "cv", "nstrat"),
  file = paste0(output_path, "03_optimization_sampling_frame.RData")
)
message("Saved: ", paste0(output_path, "03_optimization_sampling_frame.RData"))

write.table(sample1, file.path(output_path, "03_optimization_sampling_sample.csv"),
            sep = ",", row.names = FALSE, quote = FALSE)
message("Saved: ", file.path(output_path, "03_optimization_sampling_sample.csv"))

#---------------------------------------------------------------------------
# Output assembly: points_all, tab3, figures, gpkg (Phase 7 continued)
#---------------------------------------------------------------------------

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

# Final tab
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
    dom_code, NUTS, IPCC, PEDO,
    NUTS_CLASS, IPCC_CLASS, PEDO_CLASS,
    n_cell, dom_val,
    n_tot, n_lucas_grid, n_lucas_forced
  )

write.csv(tab3, paste0(output_path, "03_optimization_sampling_units.csv"), row.names = FALSE)
message("Saved: ", paste0(output_path, "03_optimization_sampling_units.csv"))
write.csv(points_all, paste0(output_path, "03_optimization_sampling_allpoints.csv"), row.names = FALSE)
message("Saved: ", paste0(output_path, "03_optimization_sampling_allpoints.csv"))

# Save domains as geopackage
stack5_sf <- st_as_sf(stack5, coords = c("x", "y"), crs = map_crs)
stack5_sf <- stack5_sf %>%
  left_join(tab3[, c("dom_code", "dom_val", "n_tot", "n_lucas_grid", "n_lucas_forced")],
            by = "dom_code")
st_write(stack5_sf, paste0(output_path, "03_optimization_sampling.gpkg"),
         delete_dsn = TRUE)
message("Saved: ", paste0(output_path, "03_optimization_sampling.gpkg"))

#################################################################################
# Figures
pdf(paste0(output_path, "03_optimization_sampling_figures.pdf"), width = 8.5, height = 8)

# CV boxplot
boxplot(cv_out,
        col  = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
        main = paste0("Coefficient of Variation Distribution (sample = ", nrow(points_all), ")"),
        ylab = "CV",
        xlab = "")
abline(a = 0.05, b = 0, col = "red")

# Recupera IPCC e PEDO da stack5 tramite cell
points_all$IPCC <- stack5$IPCC[match(points_all$cell, stack5$cell)]
points_all$PEDO <- stack5$PEDO[match(points_all$cell, stack5$cell)]

# Convert to sf object
s_sf1 <- st_as_sf(points_all, coords = c("x", "y"), crs = map_crs) %>%
  st_transform(crs = 4326) %>%
  mutate(IPCC = as.factor(IPCC),
         PEDO = as.factor(PEDO))

# Map by IPCC
print(
  ggplot(s_sf1) +
    geom_sf(aes(color = IPCC), size = 0.5, alpha = 0.5) +
    labs(title = paste0("IPCC Land Use Classes (sample = ", nrow(points_all), ")"),
         color = "IPCC") +
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
  # Converti coordinate WGS84 -> map CRS fuori da mutate
  lf_v <- project(
    vect(lucas_points_filt, geom = c("X_WGS84", "Y_WGS84"), crs = "EPSG:4326"),
    map_crs
  )
  lf_coords <- as.data.frame(crds(lf_v))
  colnames(lf_coords) <- c("x", "y")
  lf_coords$type <- "LUCAS Discarded"
  pts_all_plot <- bind_rows(pts_all_plot, lf_coords)
}

pts_all_sf <- st_as_sf(pts_all_plot, coords = c("x", "y"), crs = map_crs) %>%
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
message("Saved: ", paste0(output_path, "03_optimization_sampling_figures.pdf"))

# ---------------------------------------------------------------------------
# Export strati come poligoni per la cartella Files/Strata/ del tool SML
# ---------------------------------------------------------------------------
dv_to_str <- dom_code_to_num %>%
  left_join(domain_map, by = c("dom_code" = "old")) %>%
  dplyr::select(domainvalue = new, dom_code_str) %>%
  distinct()

strata_polygons <- framenew_samp %>%
  left_join(dv_to_str, by = c("DOMAINVALUE" = "domainvalue")) %>%
  st_as_sf(coords = c("x", "y"), crs = map_crs) %>%
  st_buffer(dist = res_m / 2, endCapStyle = "SQUARE", joinStyle = "MITRE") %>%
  group_by(dom_code_str, STRATO) %>%
  summarise(.groups = "drop") %>%
  mutate(stratum_id = paste(dom_code_str, STRATO, sep = "_"))

st_write(strata_polygons,
         file.path(output_path, "03_optimization_sampling_strata.gpkg"),
         delete_dsn = TRUE)
message("Saved: ", file.path(output_path, "03_optimization_sampling_strata.gpkg"))
