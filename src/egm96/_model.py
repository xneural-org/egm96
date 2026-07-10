"""Loading and interpolation of the pre-sampled EGM96 grid."""

from __future__ import annotations

import functools
import importlib.resources

import numpy as np
import numpy.typing as npt

#: Name of the grid file bundled inside :mod:`egm96.data`.
_GRID_FILE = "egm96-5.npz"

#: Grid spacing in degrees (5 arc-minutes). Kept in sync with the generated data
#: file and used only for documentation/introspection, the actual spacing is
#: derived from the stored axes at runtime.
GRID_RESOLUTION_DEG = 5.0 / 60.0

#: Empirically measured worst-case bilinear-interpolation error against the full
#: 5' EGM96 model, over 6000 random global points.
MAX_INTERPOLATION_ERROR_M = 0.06


class _Grid:
    """A regularly spaced, global latitude/longitude grid of undulations.

    :param latitudes: Ascending latitude axis in degrees.
    :param longitudes: Ascending longitude axis in degrees.
    :param values: ``(latitudes.size, longitudes.size)`` array of undulations in meters,
        where ``values[i, j]`` corresponds to ``latitudes[i]``/``longitudes[j]``.
    """

    __slots__ = (
        "_dlatitude",
        "_dlongitude",
        "_latitude_0",
        "_longitude_0",
        "latitudes",
        "longitudes",
        "values",
    )

    def __init__(
        self,
        latitudes: npt.NDArray[np.float64],
        longitudes: npt.NDArray[np.float64],
        values: npt.NDArray[np.float32],
    ) -> None:
        self.latitudes = latitudes
        self.longitudes = longitudes
        self.values = values
        self._latitude_0 = float(latitudes[0])
        self._longitude_0 = float(longitudes[0])
        self._dlatitude = float(latitudes[1] - latitudes[0])
        self._dlongitude = float(longitudes[1] - longitudes[0])

    def interpolate(
        self, latitude: npt.NDArray[np.float64], longitude: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Bilinearly interpolate the grid at the given coordinates.

        Longitude is wrapped modulo 360 degrees onto the grid's span and
        latitude is clamped to the grid's range before interpolation.

        Why bilinear is enough
        ----------------------
        The grid is the full EGM96 model at its native 5' (~9 km) spacing, and
        the geoid is a very smooth field at that scale -- it has no structure
        finer than the model itself, so a straight-line blend between two
        adjacent samples tracks it closely. Measured against GeographicLib's
        cubic evaluation of the same model, bilinear here is off by at most
        <0.06 m anywhere on Earth (see the tests). That is an order of magnitude
        below EGM96's own ~0.5-1 m uncertainty relative to the true geoid, so a
        higher-order scheme (bicubic, or pygeodesy's cubic) would only refine a
        number that is already far more precise than the model it comes from.
        Bilinear also stays a handful of vectorized array ops, so array inputs
        are fast and it needs nothing beyond numpy.

        :param latitude: Latitude(s) in degrees.
        :param longitude: Longitude(s) in degrees, broadcastable against ``latitude``.
        :returns: Interpolated undulation(s) in meters.
        :rtype: numpy.ndarray
        """
        # Longitude is periodic: fold it into the grid's [lon0, lon0 + 360) span
        # so e.g. +181 deg maps onto -179 deg. The grid carries a duplicate +180
        # column (identical to -180), so points in the last cell interpolate
        # across the seam without special-casing.
        longitude = self._longitude_0 + np.mod(longitude - self._longitude_0, 360.0)
        # Latitude does not wrap; clamp to [-90, 90] so the poles stay in range.
        latitude = np.clip(latitude, self._latitude_0, self.latitudes[-1])

        # Fractional grid coordinates: fi/fj are the (real-valued) row/column
        # positions of the query point, e.g. fi = 3.7 lies 70% of the way from
        # row 3 to row 4.
        fi = (latitude - self._latitude_0) / self._dlatitude
        fj = (longitude - self._longitude_0) / self._dlongitude

        # Index of the lower-left node of the enclosing cell. Clamp to size - 2
        # so the +1 neighbors below always exist (also guards the clamped poles /
        # the exact +180 edge, where fi/fj can land on the last node).
        i0 = np.clip(np.floor(fi).astype(np.intp), 0, self.latitudes.size - 2)
        j0 = np.clip(np.floor(fj).astype(np.intp), 0, self.longitudes.size - 2)
        # Position within the cell, in [0, 1]: di along latitude, dj along longitude.
        di = fi - i0
        dj = fj - j0

        # Weighted average of the four surrounding nodes; the four weights are
        # the areas of the opposite sub-rectangles and sum to 1. Reduces to the
        # exact node value when the point sits on a node (di or dj == 0).
        return (
            self.values[i0, j0] * (1.0 - di) * (1.0 - dj)  # lower-left
            + self.values[i0, j0 + 1] * (1.0 - di) * dj  # lower-right
            + self.values[i0 + 1, j0] * di * (1.0 - dj)  # upper-left
            + self.values[i0 + 1, j0 + 1] * di * dj  # upper-right
        )


def _decode(second_diff: npt.NDArray[np.int16]) -> npt.NDArray[np.int32]:
    """Invert the second-difference encoding produced by the grid generator.

    The generator differences the centimeter grid along longitude and then
    along latitude; two cumulative sums (in int32, to avoid overflow) recover
    the original values exactly.

    :param second_diff: The stored second-difference array in int16 centimeters.
    :returns: The undulation grid in int32 centimeters.
    :rtype: numpy.ndarray
    """
    cm = np.cumsum(second_diff.astype(np.int32), axis=0)  # undo latitude diff
    np.cumsum(cm, axis=1, out=cm)  # undo longitude diff
    return cm


@functools.lru_cache(maxsize=1)
def _load_grid() -> _Grid:
    """Load and cache the bundled grid.

    Called at most once per process thanks to :func:`functools.lru_cache`.

    :returns: The bundled EGM96 grid.
    :rtype: _Grid
    """
    resource = importlib.resources.files("egm96.data").joinpath(_GRID_FILE)
    with resource.open("rb") as fh, np.load(fh) as data:
        latitudes = data["latitudes"].astype(np.float64)
        longitudes = data["longitudes"].astype(np.float64)
        # Stored as second differences of int16 centimeters to keep the file
        # small; decode and recover meters.
        values = _decode(data["undulation_d2"]).astype(np.float32) / 100.0
    return _Grid(latitudes, longitudes, values)


def undulation(
    latitude: npt.ArrayLike, longitude: npt.ArrayLike
) -> float | npt.NDArray[np.float64]:
    """Return the EGM96 geoid undulation in meters at ``latitude``/``longitude``.

    The undulation is the height of the geoid (mean sea level) above the WGS84
    ellipsoid. To convert a GPS ellipsoidal height ``h`` to an orthometric
    (above-mean-sea-level) height ``H``::

        H = h - undulation(lat, lon)

    Longitude is wrapped modulo 360 degrees and latitude is clamped to
    ``[-90, 90]``.

    :param latitude: Latitude in degrees. A scalar or an array-like broadcastable
        against ``longitude``.
    :param longitude: Longitude in degrees. A scalar or an array-like broadcastable
        against ``latitude``.
    :returns: The undulation in meters - a ``float`` for scalar inputs, or an
        ``ndarray`` when either input is array-like.
    :rtype: float | numpy.ndarray
    """
    latitudes = np.asarray(latitude, dtype=np.float64)
    longitudes = np.asarray(longitude, dtype=np.float64)
    result = _load_grid().interpolate(latitudes, longitudes)
    if np.isscalar(latitude) and np.isscalar(longitude):
        return float(result)
    return result
