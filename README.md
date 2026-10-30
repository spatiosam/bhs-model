# Geographic Body Heat Storage (BHS) Model

A small R implementation of a spatially explicit human body heat-storage model. The repository is intentionally minimal: users provide aligned raster layers for the environmental variables, specify body and activity assumptions, and receive raster outputs in watts, degrees Celsius per hour, and degrees Celsius per meter walked.

This version is designed as a transparent research script rather than an R package. It omits land-surface temperature (LST), land cover, canopy, sky-view factor, route reconstruction, and other project-specific inputs.

A versioned release of the model is permanently archived in Zenodo: https://doi.org/10.5281/zenodo.23816752

## What the model uses

Required raster inputs:

|Input|Symbol|Units|Required|
|-|-|-:|-|
|Ambient air temperature|`Ta`|°C|Yes|
|Mean radiant / radiant temperature|`Trad`|°C|Yes|
|Relative humidity|`RH`|%|Yes|
|Wind / airflow|`wind`|m/s|Yes|
|Digital elevation model|`dem`|m|No|
|Signed terrain grade|`grade`|fraction|No|

All required rasters must have the same CRS, extent, resolution, and row/column geometry.

`Trad` is supplied directly. The script does **not** use LST as a proxy for radiant temperature.

For terrain-aware activity, a signed grade raster is preferred when the direction of movement is known. Positive values represent uphill movement and negative values downhill movement. A DEM alone can provide slope magnitude but cannot determine whether a traveler is moving uphill or downhill.

## Body and activity settings

The user-editable block near the top of `BHS\_model.R` contains all normal settings.

```r
BODY <- list(
  mass\_kg = 70,
  BSA\_m2 = 1.80,
  load\_kg = 0
)

ACTIVITY <- list(
  mode = "fixed\_M",
  fixed\_M\_W = 120,
  fixed\_VO2\_ml\_kg\_min = 10,
  walking\_speed\_m\_s = 0,
  fixed\_external\_work\_W = NA\_real\_
)
```

Three activity modes are available:

* `fixed\_M`: directly specify metabolic heat production in watts.
* `fixed\_vo2`: specify VO2 in mL kg^-1 min^-1 and allow the script to convert it to metabolic heat.
* `terrain\_vo2`: estimate VO2 from the supplied walking speed and terrain grade using the uphill/level and downhill formulations used in the exertion implementation.

Carried load is included in the external mechanical work term when work is calculated from movement and terrain. The model does not impose an additional uncited metabolic load multiplier. If metabolic cost while carrying a load is known or assumed, use `fixed\_vo2` or `fixed\_M`.

## Clothing settings

```r
CLOTHING <- list(
  fcl = 1.31,
  Re\_cl\_m2\_kPa\_W = 0,
  max\_skin\_wettedness = 0.90
)
```

These values are intentionally exposed. Users should edit them when a study has defensible clothing or skin-wettedness assumptions.

## Outputs

The main GeoTIFF contains three layers:

|Layer|Units|Interpretation|
|-|-:|-|
|`BHS\_W`|W|Total-body heat-storage rate|
|`dTcore\_C\_per\_h`|°C/h|Estimated temperature change per hour|
|`dTcore\_C\_per\_m`|°C/m|Estimated temperature change per meter walked|

Positive BHS indicates net heat storage; negative BHS indicates net heat loss.

`dTcore\_C\_per\_m` requires a walking speed greater than zero. When `walking\_speed\_m\_s = 0`, that output layer is `NA`.

The script can optionally write intermediate components such as metabolic heat, external work, skin temperature, vapor pressure, convective exchange, radiative exchange, respiratory heat loss, and evaporative capacity by setting:

```r
OUTPUT$write\_components <- TRUE
```

## Installation

The only required R package is [`terra`](https://rspatial.org/terra/).

```r
install.packages("terra")
```

R 4.x is recommended.

## Quick start

1. Clone or download this repository.
2. Open `BHS\_model.R`.
3. Edit the `INPUT`, `OUTPUT`, `BODY`, `ACTIVITY`, and `CLOTHING` objects at the top of the file.
4. Run:

```bash
Rscript BHS\_model.R
```

Or source the model from another R script:

```r
source("BHS\_model.R")

run\_bhs\_model(
  INPUT = INPUT,
  OUTPUT = OUTPUT,
  BODY = BODY,
  ACTIVITY = ACTIVITY,
  CLOTHING = CLOTHING
)
```

See `examples/example\_run.R` for a fully reproducible synthetic example.

## Important implementation choices

This repository is a **simplified implementation**, not a verbatim reproduction of any single published analysis.

The core heat-balance structure follows the published BHS/ΔTcore lineage:

\[
BHS = M - W + C + R - (C\_{res} + E\_{res}) - E\_{max}
\]

where `M` is metabolic heat production, `W` is external work, `C` and `R` are convective and radiative exchange, `Cres` and `Eres` are respiratory heat exchange, and `Emax` is evaporative heat-loss capacity.

The temperature-rate output is calculated directly from total heat storage and whole-body heat capacity:

\[
\\Delta T / \\Delta t =
\\frac{BHS}{m c}
\]

and converted to °C/h by multiplying by 3600 seconds per hour.

For a person walking at speed (v),

\[
\\Delta T / \\Delta x =
\\frac{BHS}{m c v}
\]

which produces °C per meter.

Radiant temperature is supplied directly and radiative exchange is calculated against the supplied `Trad`/MRT rather than deriving air temperature or radiant temperature from LST.

## Methodological lineage

The repository intentionally limits its methodological citations to five sources.

**Physiological heat-exchange basis**

Ravanelli, N., Bongers, C. C. W. G., \& Jay, O. (2019). The Biophysics of Human Heat Exchange. In *Heat Stress in Sport and Exercise* (pp. 29–43). Springer. https://doi.org/10.1007/978-3-319-93515-7\_2

**Original geographic ΔTcore model and cost-distance implementation in this model lineage**

Chambers, S. N., McMahan, B., \& Bongers, C. C. W. G. (2020). Developing a geospatial measure of change in core temperature for migrating persons in the Mexico-U.S. border region. *Spatial and Spatio-temporal Epidemiology, 35*, 100363. https://doi.org/10.1016/j.sste.2020.100363

**Geographic body heat-storage raster model**

Chambers, S. N., Brown, H. E., Keith, L., \& Austhof, E. (2023). Development of a geographic human heat balance equation to support public health analyses: An Arizona urban sun corridor application. *Remote Sensing Applications: Society and Environment, 32*, 101009. https://doi.org/10.1016/j.rsase.2023.101009

**Terrain/exertion extension**

Chambers, S. N., Boyce, G. A., Martínez, D. E., Bongers, C. C. W. G., \& Keith, L. (2023). The contribution of physical exertion to heat-related illness and death in the Arizona borderlands. *Spatial and Spatio-temporal Epidemiology, 46*, 100590. https://doi.org/10.1016/j.sste.2023.100590

**Spatiotemporal radiation and heat-exposure extension**

Chambers, S. N. (2027). Reconstructing the environmental and temporal contexts of heat-related deaths on the Arizona Borderscape through a multivariate typology of heat exposure. *Applied Geography*, 104219. https://doi.org/10.1016/j.apgeog.2026.104219

## Scope and limitations

This model estimates heat storage under the supplied environmental, physiological, clothing, and activity assumptions. It is not a clinical prediction model and should not be interpreted as an individual diagnosis or a precise forecast of a person's core temperature.

Important limitations include:

* input rasters represent environmental conditions at their own spatial and temporal resolutions;
* a DEM alone cannot distinguish uphill from downhill travel direction;
* a single walking speed is applied unless the user constructs a spatial speed surface outside the model;
* hydration, acclimatization, fatigue, changing clothing, behavioral thermoregulation, and individual clinical factors are not dynamically modeled;
* °C/m represents the local rate of estimated temperature change per distance traveled, not cumulative route exposure unless integrated along a path;
* parameter choices should be justified for the population and activity being studied.

## Repository structure

```text
bhs-model/
├── BHS\_model.R
├── README.md
├── LICENSE
├── CITATION.cff
├── CHANGELOG.md
├── .gitignore
├── examples/
│   └── example\_run.R
└── data/
    └── README.md
```

## License

Code in this repository is released under the MIT License. The cited publications retain their respective copyrights and licenses.

## Citation

GitHub will display a **Cite this repository** option from `CITATION.cff`.

A versioned release of the software is permanently archived in Zenodo:

Chambers, S. N. (2026). *Geographic Body Heat Storage (BHS) Model* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23816752

When using the model scientifically, cite the archived software release together with the methodological publication or publications relevant to the implementation being used.

The GitHub repository contains the current development version of the model, while the Zenodo record provides a permanent, citable archive of the released software version.
