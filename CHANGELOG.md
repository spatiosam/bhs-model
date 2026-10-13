# Changelog

All notable changes to this repository will be documented here.

The project uses semantic versioning where practical.

## 0.1.0 — 2026-09-30

Initial public-ready bare-bones implementation.

### Added

- raster inputs for ambient temperature, radiant/MRT temperature, relative humidity, and wind;
- optional DEM or signed-grade input;
- editable body mass, body surface area, and carried-load settings;
- fixed metabolic heat, fixed VO2, and terrain-dependent VO2 modes;
- optional fixed external work;
- editable clothing area factor, evaporative clothing resistance, and maximum skin wettedness;
- body heat storage output in watts;
- estimated temperature-change output in °C/h;
- estimated temperature-change output in °C/m when walking speed is greater than zero;
- optional intermediate-component rasters;
- runnable synthetic example;
- README, citation metadata, MIT license, and Git ignore rules.

### Scope

This release intentionally excludes LST-derived temperature, land cover, canopy,
SVF, routing, cost-distance accumulation, and other project-specific extensions.
