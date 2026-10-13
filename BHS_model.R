# =============================================================================
# Bare-Bones Geographic Body Heat Storage (BHS) Model
# Version 0.1.0
# =============================================================================
#
# Required raster inputs:
#   Ta    ambient air temperature (deg C)
#   Trad  mean radiant / radiant temperature (deg C)
#   RH    relative humidity (%)
#   wind  air velocity (m/s)
#
# Optional terrain input:
#   dem   elevation (m), used to derive unsigned slope magnitude
#   grade signed grade as a fraction (+ uphill, - downhill); preferred when
#         direction of travel is known
#
# Outputs:
#   BHS_W            body heat-storage rate (W)
#   dTcore_C_per_h   estimated core-temperature change (deg C / h)
#   dTcore_C_per_m   estimated core-temperature change (deg C / m walked)
#
# This is a deliberately simplified raster implementation. It does not use
# LST, land cover, canopy, SVF, routing, or cumulative cost distance.
#
# Methodological lineage:
#
# Ravanelli N, Bongers CCWG, Jay O. 2019. The Biophysics of Human Heat
# Exchange. In: Heat Stress in Sport and Exercise, pp. 29-43.
# https://doi.org/10.1007/978-3-319-93515-7_2
#
# Chambers SN, McMahan B, Bongers CCWG. 2020. Developing a geospatial measure
# of change in core temperature for migrating persons in the Mexico-U.S.
# border region. Spatial and Spatio-temporal Epidemiology 35:100363.
# https://doi.org/10.1016/j.sste.2020.100363
#
# Chambers SN, Brown HE, Keith L, Austhof E. 2023. Development of a geographic
# human heat balance equation to support public health analyses: An Arizona
# urban sun corridor application. Remote Sensing Applications: Society and
# Environment 32:101009.
# https://doi.org/10.1016/j.rsase.2023.101009
#
# Chambers SN, Boyce GA, Martinez DE, Bongers CCWG, Keith L. 2023. The
# contribution of physical exertion to heat-related illness and death in the
# Arizona borderlands. Spatial and Spatio-temporal Epidemiology 46:100590.
# https://doi.org/10.1016/j.sste.2023.100590
#
# NOTE:
# This script is a simplified implementation based on the published model
# lineage above; it is not a verbatim reproduction of any one paper.
# =============================================================================

library(terra)

# =============================================================================
# 1. USER SETTINGS
#    Most users should only need to edit this section.
# =============================================================================

INPUT <- list(
  Ta    = "path/to/ambient_temperature_C.tif",
  Trad  = "path/to/MRT_or_radiant_temperature_C.tif",
  RH    = "path/to/relative_humidity_percent.tif",
  wind  = "path/to/wind_m_s.tif",
  dem   = NA_character_,   # optional
  grade = NA_character_    # optional signed fraction; preferred over DEM
)

OUTPUT <- list(
  file = "bhs_outputs.tif",
  write_components = FALSE,
  components_file = "bhs_components.tif"
)

BODY <- list(
  mass_kg = 70,
  BSA_m2 = 1.80,
  load_kg = 0
)

ACTIVITY <- list(
  mode = "fixed_M",             # "fixed_M", "fixed_vo2", "terrain_vo2"
  fixed_M_W = 120,
  fixed_VO2_ml_kg_min = 10,
  walking_speed_m_s = 0,        # > 0 required for deg C/m
  fixed_external_work_W = NA_real_
)

CLOTHING <- list(
  fcl = 1.31,                   # clothing area factor
  Re_cl_m2_kPa_W = 0,           # evaporative clothing resistance
  max_skin_wettedness = 0.90
)

# =============================================================================
# 2. MODEL FUNCTION
# =============================================================================

run_bhs_model <- function(
  INPUT,
  OUTPUT,
  BODY,
  ACTIVITY,
  CLOTHING
) {

  BODY_SPECIFIC_HEAT_J_KG_C <- 3470
  SIGMA <- 5.670374419e-8
  SKIN_EMISSIVITY <- 0.95
  RADIANT_AREA_RATIO <- 0.72

  read_one <- function(path, label) {
    if (length(path) != 1 || is.na(path) || !nzchar(path)) {
      stop(label, " path is missing.")
    }
    x <- rast(path)
    if (nlyr(x) != 1) stop(label, " must be a single-layer raster.")
    x
  }

  check_geom <- function(template, x, label) {
    ok <- compareGeom(
      template, x,
      crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
      stopOnError = FALSE
    )
    if (!isTRUE(ok)) {
      stop(
        label, " does not match the Ta raster geometry.\n",
        "Project/resample/align the input before running the model."
      )
    }
  }

  psat_kPa <- function(temp_C) {
    exp(18.956 - (4030.18 / (temp_C + 235.0))) / 10.0
  }

  if (!is.numeric(BODY$mass_kg) || BODY$mass_kg <= 0) {
    stop("BODY$mass_kg must be > 0.")
  }
  if (!is.numeric(BODY$BSA_m2) || BODY$BSA_m2 <= 0) {
    stop("BODY$BSA_m2 must be > 0.")
  }
  if (!is.numeric(BODY$load_kg) || BODY$load_kg < 0) {
    stop("BODY$load_kg must be >= 0.")
  }
  if (!is.numeric(ACTIVITY$walking_speed_m_s) ||
      ACTIVITY$walking_speed_m_s < 0) {
    stop("ACTIVITY$walking_speed_m_s must be >= 0.")
  }

  # ---------------------------------------------------------------------------
  # Climate inputs
  # ---------------------------------------------------------------------------

  Ta   <- read_one(INPUT$Ta,   "Ambient temperature (Ta)")
  Trad <- read_one(INPUT$Trad, "Radiant/MRT temperature (Trad)")
  RH   <- read_one(INPUT$RH,   "Relative humidity (RH)")
  wind <- read_one(INPUT$wind, "Wind/airflow")

  check_geom(Ta, Trad, "Trad")
  check_geom(Ta, RH,   "RH")
  check_geom(Ta, wind, "wind")

  RH   <- clamp(RH, lower = 0, upper = 100, values = TRUE)
  wind <- ifel(wind < 0, NA, wind)

  zero <- Ta * 0

  # ---------------------------------------------------------------------------
  # Terrain grade
  # ---------------------------------------------------------------------------

  grade <- zero

  if (!is.na(INPUT$grade) && nzchar(INPUT$grade)) {

    grade <- read_one(INPUT$grade, "Signed grade")
    check_geom(Ta, grade, "grade")

  } else if (!is.na(INPUT$dem) && nzchar(INPUT$dem)) {

    dem <- read_one(INPUT$dem, "DEM")
    check_geom(Ta, dem, "DEM")

    # A DEM gives slope magnitude but not direction of travel.
    slope_rad <- terrain(dem, v = "slope", unit = "radians", neighbors = 8)
    grade <- tan(slope_rad)
  }

  grade_percent <- grade * 100
  speed_m_s <- ACTIVITY$walking_speed_m_s
  speed_m_min <- speed_m_s * 60

  # ---------------------------------------------------------------------------
  # Metabolic heat production
  # ---------------------------------------------------------------------------

  VO2_ml_kg_min <- NULL

  if (ACTIVITY$mode == "fixed_M") {

    M_W <- zero + ACTIVITY$fixed_M_W

  } else if (ACTIVITY$mode == "fixed_vo2") {

    VO2_ml_kg_min <- zero + ACTIVITY$fixed_VO2_ml_kg_min
    VO2_L_min <- VO2_ml_kg_min * BODY$mass_kg / 1000
    M_W <- VO2_L_min * 341

  } else if (ACTIVITY$mode == "terrain_vo2") {

    if (speed_m_s <= 0) {
      stop("terrain_vo2 requires ACTIVITY$walking_speed_m_s > 0.")
    }

    # Uphill / level:
    # VO2 = (speed * 0.2) + (speed * grade * 0.9) + 3.5
    # speed in m/min; grade as a fraction.
    grade_up <- ifel(grade > 0, grade, 0)

    VO2_up <- 0.2 * speed_m_min +
              0.9 * speed_m_min * grade_up +
              3.5

    # Downhill:
    # VO2 = 6.8192 + (speed * 0.1313) + (1.2367 * grade)
    # Here grade is expressed in percent for the downhill formulation.
    VO2_down <- 6.8192 +
                0.1313 * speed_m_min +
                1.2367 * grade_percent

    VO2_ml_kg_min <- ifel(grade < 0, VO2_down, VO2_up)
    VO2_ml_kg_min <- ifel(VO2_ml_kg_min < 0, 0, VO2_ml_kg_min)

    VO2_L_min <- VO2_ml_kg_min * BODY$mass_kg / 1000
    M_W <- VO2_L_min * 341

  } else {
    stop(
      "ACTIVITY$mode must be 'fixed_M', 'fixed_vo2', or 'terrain_vo2'."
    )
  }

  # ---------------------------------------------------------------------------
  # External mechanical work
  # ---------------------------------------------------------------------------

  if (!is.na(ACTIVITY$fixed_external_work_W)) {

    W_W <- zero + ACTIVITY$fixed_external_work_W

  } else if (speed_m_s > 0) {

    moving_mass_kg <- BODY$mass_kg + BODY$load_kg
    W_W <- 9.81 * moving_mass_kg * speed_m_s * grade

  } else {

    W_W <- zero
  }

  # ---------------------------------------------------------------------------
  # Vapor pressure and estimated skin temperature
  # ---------------------------------------------------------------------------

  Pa_kPa <- psat_kPa(Ta) * (RH / 100)

  M_W_m2 <- M_W / BODY$BSA_m2

  Tskin_C <- 30.0 +
             0.138 * Ta +
             0.254 * Pa_kPa -
             0.57 * wind +
             (1.28 * 1.28 * 0.001) * M_W_m2 -
             0.553

  # ---------------------------------------------------------------------------
  # Convective and radiative exchange
  # ---------------------------------------------------------------------------

  hc_W_m2_K <- 8.3 * (wind ^ 0.6)

  # Radiative exchange is expressed explicitly against supplied Trad/MRT so
  # the model does not require LST.
  Tmean_K <- ((Tskin_C + Trad) / 2) + 273.15

  hr_W_m2_K <- 4 *
               SKIN_EMISSIVITY *
               SIGMA *
               RADIANT_AREA_RATIO *
               (Tmean_K ^ 3)

  # Signed total-body exchange:
  # positive = heat gain; negative = heat loss.
  C_W <- hc_W_m2_K * (Ta - Tskin_C) * BODY$BSA_m2
  R_W <- hr_W_m2_K * (Trad - Tskin_C) * BODY$BSA_m2

  # ---------------------------------------------------------------------------
  # Respiratory heat loss
  # ---------------------------------------------------------------------------

  Cres_W_m2 <- 0.0014 * M_W_m2 * (34.0 - Ta)
  Eres_W_m2 <- 0.0173 * M_W_m2 * (5.87 - Pa_kPa)

  Cres_W <- Cres_W_m2 * BODY$BSA_m2
  Eres_W <- Eres_W_m2 * BODY$BSA_m2

  # ---------------------------------------------------------------------------
  # Maximum evaporative heat loss
  # ---------------------------------------------------------------------------

  Psk_sat_kPa <- psat_kPa(Tskin_C)

  # Evaporative heat-transfer coefficient from convective transfer.
  he_W_m2_kPa <- 16.5 * hc_W_m2_K

  boundary_resistance <- ifel(
    he_W_m2_kPa > 0,
    1 / (he_W_m2_kPa * CLOTHING$fcl),
    Inf
  )

  evap_resistance <- CLOTHING$Re_cl_m2_kPa_W + boundary_resistance

  Emax_W_m2 <- CLOTHING$max_skin_wettedness *
               (Psk_sat_kPa - Pa_kPa) /
               evap_resistance

  Emax_W_m2 <- ifel(Emax_W_m2 < 0, 0, Emax_W_m2)
  Emax_W <- Emax_W_m2 * BODY$BSA_m2

  # ---------------------------------------------------------------------------
  # Body heat storage
  # ---------------------------------------------------------------------------

  BHS_W <- M_W -
           W_W +
           C_W +
           R_W -
           Cres_W -
           Eres_W -
           Emax_W

  names(BHS_W) <- "BHS_W"

  # ---------------------------------------------------------------------------
  # Unit conversions
  # ---------------------------------------------------------------------------

  body_heat_capacity_J_C <- BODY$mass_kg * BODY_SPECIFIC_HEAT_J_KG_C

  dTcore_C_per_h <- BHS_W * 3600 / body_heat_capacity_J_C
  names(dTcore_C_per_h) <- "dTcore_C_per_h"

  if (speed_m_s > 0) {
    dTcore_C_per_m <- BHS_W /
      (body_heat_capacity_J_C * speed_m_s)
  } else {
    dTcore_C_per_m <- BHS_W * NA_real_
  }
  names(dTcore_C_per_m) <- "dTcore_C_per_m"

  # ---------------------------------------------------------------------------
  # Write output
  # ---------------------------------------------------------------------------

  output <- c(BHS_W, dTcore_C_per_h, dTcore_C_per_m)

  out_dir <- dirname(OUTPUT$file)
  if (!identical(out_dir, ".") && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  writeRaster(
    output,
    OUTPUT$file,
    overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE")
  )

  if (isTRUE(OUTPUT$write_components)) {

    components <- c(
      grade,
      M_W,
      W_W,
      Tskin_C,
      Pa_kPa,
      C_W,
      R_W,
      Cres_W,
      Eres_W,
      Emax_W
    )

    names(components) <- c(
      "grade_fraction",
      "M_W",
      "W_W",
      "Tskin_C",
      "Pa_kPa",
      "C_W",
      "R_W",
      "Cres_W",
      "Eres_W",
      "Emax_W"
    )

    if (!is.null(VO2_ml_kg_min)) {
      names(VO2_ml_kg_min) <- "VO2_ml_kg_min"
      components <- c(components, VO2_ml_kg_min)
    }

    comp_dir <- dirname(OUTPUT$components_file)
    if (!identical(comp_dir, ".") && !dir.exists(comp_dir)) {
      dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
    }

    writeRaster(
      components,
      OUTPUT$components_file,
      overwrite = TRUE,
      gdal = c("COMPRESS=DEFLATE")
    )
  }

  invisible(output)
}

# =============================================================================
# 3. RUN
# =============================================================================
#
# When this file is executed directly with Rscript, the USER SETTINGS above
# are used. When this file is sourced from another R script, it only defines
# run_bhs_model(), allowing examples or other workflows to provide their own
# settings.

if (sys.nframe() == 0) {
  run_bhs_model(INPUT, OUTPUT, BODY, ACTIVITY, CLOTHING)
}
