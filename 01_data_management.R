# =============================================================================
# SML sampling site selection — Script 1: Data management
#
# Handles Phase 1 of the pipeline (raster stack / domain
# definition) and Phase 2 (LUCAS point integration) of the pipeline: reads
# all raw inputs (rasters, legends, LUCAS points) and produces the cleaned
# analysis frame (stack5) plus the LUCAS-derived objects needed by the
# downstream scripts.
#
# Output (folder: intermediate/, subfolder of this script's own folder):
#   01_data_management_data.RData, consumed by 02_spatial_analysis.R and
#                                              03_optimization_sampling.R
#   01_data_management_correlation_matrix.pdf
# =============================================================================

library(sp)
library(raster)
library(dplyr)
library(tidyr)
library(SamplingStrata)
library(bethel)
library(clhs)
library(sf)
library(terra)
library(readxl)

#----------------------------------------------------------------------------------------------------
# Parameters & settings (to be modified by users)
#----------------------------------------------------------------------------------------------------
# Input data location (rasters, legends, LUCAS points) — unrelated to the
# intermediate/output hand-off folders below, left as the actual data path.
map_path <- "D:\\Google Drive\\ISPRA\\Spatial3\\SML_code_italy-main/Carte_finali_1000m_EPSG3035/"
# map_path <- "G:\\Il mio Drive\\ISPRA\\Spatial3\\SML_code_italy-main/Carte_finali_1000m_EPSG3035/"

min_point_per_domain <- 5 # if domains (soil units) are to be filtered according to a min number of points

vars_cont_y <- c("PH_FAO", "SOC_FAO")
variance_cont_y <- c("PH_VAR", "SOC_VAR")
vars_cont_x <- c("PH_JRC", "SOC_JRC")

use_lucas_points <- TRUE # TRUE if any lucas points is imposed
use_lucas_grid <- TRUE # TRUE if lucas grid is taken into account in the stack
base_names <- if (use_lucas_grid) c("LUCAS_GRID", "NUTS", "IPCC", "PEDO") else c("NUTS", "IPCC", "PEDO")
if (use_lucas_grid) {
  vars_to_stack_1 <- c("lucas_grid_3035_1000m_rif", "NUTS1_3035_1000m",
                       "CLC2018_6CL_3035_1000m", "SR_CE_3035_1000m")
} else {
  vars_to_stack_1 <- c("NUTS1_3035_1000m",
                       "CLC2018_6CL_3035_1000m", "SR_CE_3035_1000m")
}
vars_to_stack_2 <- c("GSASS_PH_3035_1000m", "GSOC_3035_1000m",
                     "GSASS_PH_VAR_3035_1000m", "GSOC_VAR_3035_1000m",
                     "JRC_PH_3035_1000m", "JRC_OC_3035_1000m")

vars_to_stack <- c(vars_to_stack_1, vars_to_stack_2)

# Hand-off folder for the .RData/.csv/.pdf outputs consumed by
# 02_spatial_analysis.R and 03_optimization_sampling.R — a subfolder of this
# script's own folder.
intermediate_path <- "intermediate/"
dir.create(intermediate_path, recursive = TRUE, showWarnings = FALSE)


pedo_legend   <- read.csv(paste0(map_path, "legend_SR_CE.csv"), stringsAsFactors = FALSE)
ipcc_legend <- read.csv(paste0(map_path, "legend_CLC2018_6CL.csv"), stringsAsFactors = FALSE)
nut_legend <- read.csv(paste0(map_path, "legend_NUTS1.csv"), stringsAsFactors = FALSE)

#--------------------------------------------------------------------------------------------------
#Load raster stack and legends
#--------------------------------------------------------------------------------------------------
all_rast<- list.files(path = map_path, pattern = ".tif$", full.names = T)
all_rast_t<-lapply(all_rast,rast)
names(all_rast_t) <- tools::file_path_sans_ext(basename(all_rast))
vars_to_stack_r <- (all_rast_t[vars_to_stack])
stack1 <- rast(vars_to_stack_r)
names(stack1) <- c(base_names, vars_cont_y, variance_cont_y, vars_cont_x)

# CRS of the raster stack, captured as a plain string because terra
# SpatRaster objects use external pointers and do not survive save()/load()
# across sessions — scripts 2 and 3 use map_crs wherever crs(stack1) was
# used in the original script.
map_crs <- crs(stack1)

# Pixel resolution (m), captured as a plain scalar for the same reason as
# map_crs — used by 03_optimization_sampling.R to buffer sample points into
# the strata polygon export.
res_m <- res(stack1)[1]


#----------------------------------------------------------------------------------------------------
# Define and clean soil units
#----------------------------------------------------------------------------------------------------
stack2_full <- as.data.frame(stack1, xy = TRUE, cells = TRUE)
stack2 <- stack2_full[!is.na(stack2_full$IPCC) & !is.na(stack2_full$NUTS), ]  #pedo maps is not complete
stack2$dom_code <- paste(stack2$NUTS,stack2$IPCC,stack2$PEDO,sep="*")
stack2$id <- c(1:nrow(stack2))
head(stack2 , 5)

stack3 <- stack2 %>%
  left_join(nut_legend,  by = "NUTS") %>%
  left_join(ipcc_legend, by = "IPCC") %>%
  left_join(pedo_legend, by = "PEDO") %>%
  dplyr::select(
    id, cell, dom_code, x, y,
    dplyr::all_of(base_names),
    NUTS_CLASS, IPCC_CLASS, PEDO_CLASS,
    dplyr::all_of(c(vars_cont_y, variance_cont_y, vars_cont_x))
  )
head(stack3 , 5)


tab1 <- stack3 %>%
  group_by(dom_code) %>%
  summarise(
    NUTS = first(NUTS),
    IPCC = first(IPCC),
    PEDO = first(PEDO),
    NUTS_CLASS = first(NUTS_CLASS),
    IPCC_CLASS = first(IPCC_CLASS),
    PEDO_CLASS = first(PEDO_CLASS),
    n_point = n(),
    .groups = "drop"
  )


# Clean for missing continuous variables
stack4 <- stack3[complete.cases(stack3[, c(vars_cont_y, variance_cont_y, vars_cont_x)]), ]


# Filter cells associated with continuous variables
tab1_valid <- stack4 %>%
  group_by(dom_code) %>%
  summarise(n_point_var_cont = n(), .groups = "drop")

tab2 <- tab1 %>%
  left_join(tab1_valid, by = "dom_code")
tab2$n_point_var_cont[is.na(tab2$n_point_var_cont)] <- 0


# Filter domains with less than the min points per domain
tab2$dom_val <- tab2$n_point_var_cont >= min_point_per_domain
valid_domains <- tab2$dom_code[tab2$dom_val]
stack5 <- stack4[stack4$dom_code %in% valid_domains, ]
head(stack5)

# Print domain statistics
message("Total number of domains = ", nrow(tab2))
message("Largest domain as n.cell = ", max(tab2$n_point))
message("Q1 domain as n.cell = ", quantile(tab2$n_point, 0.25))
message("Q2 domain as n.cell = ", quantile(tab2$n_point, 0.50))
message("Q3 domain as n.cell = ", quantile(tab2$n_point, 0.75))
message("Number of valid domains = ", sum(tab2$dom_val))
message("Number of total Settlements domains = ", sum(tab2$IPCC == 1))
message("Number of valid Settlements domains = ", sum(tab2$IPCC == 1 & tab2$dom_val))


#--------------------------------------------------------------------------------------------------
# Correlation matrix of continuous variables in stack5
#--------------------------------------------------------------------------------------------------
cor_vars <- c(vars_cont_x,vars_cont_x)
cor_mat  <- cor(stack5[, cor_vars], use = "pairwise.complete.obs")
message("Correlation matrix (PH_FAO, SOC_FAO, PH_JRC, OC_JRC):")
print(round(cor_mat, 3))

# Correlation matrix plot (colour-coded with coefficients)
plot_cor_matrix <- function(m) {
  n   <- ncol(m)
  pal <- colorRampPalette(c("#B2182B", "white", "#2166AC"))(201)
  op  <- par(mar = c(1, 5, 5, 1))
  on.exit(par(op))
  image(1:n, 1:n, t(m[n:1, ]), col = pal, zlim = c(-1, 1),
        axes = FALSE, xlab = "", ylab = "",
        main = "Correlation matrix - stack5")
  axis(3, at = 1:n, labels = colnames(m), tick = FALSE, las = 1)
  axis(2, at = 1:n, labels = rev(rownames(m)), tick = FALSE, las = 1)
  for (i in 1:n) for (j in 1:n) {
    text(j, n - i + 1, sprintf("%.2f", m[i, j]),
         col = ifelse(abs(m[i, j]) > 0.6, "white", "black"))
  }
}
plot_cor_matrix(cor_mat)

# Save the correlation matrix plot to PDF
pdf(paste0(intermediate_path, "01_data_management_correlation_matrix.pdf"), width = 6, height = 6)
plot_cor_matrix(cor_mat)
dev.off()


#--------------------------------------------------------------------------------------------------
# Read Lucas points to be retained
#--------------------------------------------------------------------------------------------------
if (use_lucas_points) {
  lucas_raw <- read_excel(
   "D:\\Google Drive\\ISPRA\\Spatial3\\SML_code_italy-main\\IT.xlsx",
    sheet = "Points"
  )
  plot(lucas_raw$X_WGS84, lucas_raw$Y_WGS84)

  lucas_soil <- subset(lucas_raw, SOIL == 1)

  lucas_v <- project(
    vect(lucas_soil, geom = c("X_WGS84", "Y_WGS84"), crs = "EPSG:4326"),
    crs(stack1)
  )

  lucas_pts <- data.frame(
    lucas_soil,
    crds(lucas_v) |> as.data.frame() |> setNames(c("x", "y")),
    as.data.frame(terra::extract(stack1, lucas_v, cells = TRUE, ID = FALSE))
  ) |>
    transform(dom_code = stack5$dom_code[match(cell, stack5$cell)])

  lucas_pts_fixed_val <- subset(lucas_pts,
#                                (GRASSLAND == 1 | LF == 1) &
!is.na(dom_code))

  ok <- complete.cases(lucas_pts_fixed_val[, c(vars_cont_y, variance_cont_y, vars_cont_x)])
  lucas_points      <- lucas_pts_fixed_val[ ok, ]
  lucas_points_filt <- lucas_pts_fixed_val[!ok, ]

  message(sprintf(
    "Lucas points summary:\n  - Outside valid domains (no dom_code): %d\n  - Excluded (missing continuous vars):  %d\n  - Total retained:                      %d\n  - Grand total (GRASSLAND==1 | LF==1):  %d",
    sum(is.na(lucas_pts$dom_code[lucas_pts$GRASSLAND == 1 | lucas_pts$LF == 1])),
    nrow(lucas_points_filt),
    nrow(lucas_points),
    nrow(lucas_pts_fixed_val) + sum(is.na(lucas_pts$dom_code[lucas_pts$GRASSLAND == 1 | lucas_pts$LF == 1]))
  ))

  lucas_points$id_lucas     <- seq_len(nrow(lucas_points))
  lucas_points$dom_code_num <- as.integer(factor(lucas_points$dom_code,
                                                 levels = levels(factor(stack5$dom_code))))

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



  stack5_samp <- subset(stack5, !cell %in% unique(lucas_points$cell))
}


#--------------------------------------------------------------------------------------------------
# Save outputs for script 2 (spatial analysis) and script 3 (optimization & sampling)
#--------------------------------------------------------------------------------------------------
save(
  list = c(
    "stack5", "tab2", "nut_legend", "map_crs", "res_m",
    "vars_cont_y", "variance_cont_y", "vars_cont_x", "base_names",
    "use_lucas_points", "use_lucas_grid", "min_point_per_domain",
    if (use_lucas_points) c("stack5_samp", "lucas_points", "lucas_points_filt", "frame_cens")
  ),
  file = paste0(intermediate_path, "01_data_management_data.RData")
)
message("Saved: ", paste0(intermediate_path, "01_data_management_data.RData"))
