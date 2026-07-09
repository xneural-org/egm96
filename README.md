# egm96

GPS reports altitude as height above the WGS 84 reference ellipsoid, a smooth
mathematical surface that can sit tens of meters away from actual sea level.
[EGM96](https://en.wikipedia.org/wiki/Earth_Gravitational_Model), the 1996 Earth
Gravitational Model from the U.S. National Geospatial-Intelligence Agency,
describes where mean sea level lies relative to that ellipsoid. The gap between
the two is the *geoid undulation*; subtracting it from a GPS height gives an
approximate height above sea level.

This package looks up the EGM96 undulation for any coordinate.

## Usage

```python
import egm96

egm96.undulation(48.14, 11.58)   # 45.48  (meters, near Munich)
```

The undulation is the height of the geoid above the WGS 84 ellipsoid. Subtract
it from a GPS ellipsoidal height `h` to get a height above mean sea level:

```python
H = h - egm96.undulation(lat, lon)
```

Coordinates are in degrees. Longitude wraps at ±180° and latitude is clamped to
±90°. Both arguments also accept NumPy arrays, which broadcast together and
return an array.

## Accuracy and size

The package bundles the complete EGM96 5' grid (2161 × 4321 points) and
interpolates it bilinearly. Against GeographicLib's cubic evaluation of the same
model, over 6000 random points, the largest difference is < 0.06 m — negligible
next to EGM96's own ~0.5–1 m uncertainty relative to the real geoid.

Bilinear error grows quickly if the grid is coarsened, so it is worth keeping
the full resolution:

| Grid spacing      | Max error vs EGM96 |
| ----------------- | ------------------ |
| 0.25° (15')       | 0.5 m              |
| 0.17° (10')       | 0.2 m              |
| 0.083° (5', used) | 0.05 m             |

Stored plainly the grid is about 13 MB. The geoid is smooth, so encoding it as
second differences (in centimeters) and deflating the result shrinks it to
3.1 MB with no loss; the loader reconstructs the values with two cumulative
sums.

## Grid Data

`src/egm96/data/egm96-5.npz` is built from the source geoid at packaging time and
is not stored in the repository. `generate-egm96-grid` downloads the 5' PGM from
GeographicLib — only if it is not already cached locally — and writes the
`.npz`. CI runs the generator before `uv build`, so the grid ends up inside the
wheel, and installing from PyPI needs no network access at runtime.

To build locally:

```bash
uv sync
uv run ./generate-egm96-grid
uv build
```

## License

MIT (see [LICENSE](LICENSE)). EGM96 is a work of the U.S. Government (NGA, NASA,
and Ohio State University) in the public domain; the 5' dataset is redistributed
by the MIT-licensed [GeographicLib](https://geographiclib.sourceforge.io/).
