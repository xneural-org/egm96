# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0

### Added

- `egm96.undulation(lat, lon)` returning the EGM96 geoid undulation in meters
  (height of the geoid above the WGS 84 ellipsoid). Accepts scalars or NumPy
  arrays; longitude wraps at ±180° and latitude is clamped to ±90°.
- The full EGM96 5' grid, bundled inside the wheel and interpolated bilinearly.
  Measured worst-case error against the full model is < 0.06 m — well below
  EGM96's own ~0.5–1 m uncertainty relative to the true geoid.
- Compact, lossless grid storage: undulations are quantized to int16 centimeters
  and second-difference encoded, shrinking the grid from ~13 MB to ~3 MB and
  decoded with two cumulative sums at load time.
- `generate-egm96-grid`, which builds the bundled grid by decimating the
  official GeographicLib 5' PGM. The source PGM is downloaded and cached on
  demand, and the generated grid is kept out of the repository.

## [v0.1.0] - Initial Release

- No features

[v0.1.0]: https://github.com/xneural-org/egm96/releases/tag/v0.1.0
[v0.2.0]: https://github.com/xneural-org/egm96/releases/tag/v0.2.0
