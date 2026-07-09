# egm96

GPS reports altitude as height above the WGS 84 reference ellipsoid, a smooth
mathematical surface that can sit tens of meters away from actual sea level.
[EGM96](https://en.wikipedia.org/wiki/Earth_Gravitational_Model), the 1996 Earth
Gravitational Model from the U.S. National Geospatial-Intelligence Agency,
describes where mean sea level lies relative to that ellipsoid. The gap between
the two is the _geoid undulation_; subtracting it from a GPS height gives an
approximate height above sea level.

This package looks up the EGM96 undulation for any coordinate.

## License

MIT (see [LICENSE](LICENSE)). EGM96 is a work of the U.S. Government (NGA, NASA,
and Ohio State University) in the public domain; the 5′ dataset is redistributed
by the MIT-licensed [GeographicLib](https://geographiclib.sourceforge.io/).
