"""Access to the EGM96 geoid model.

The public entry point is :func:`undulation`, which returns the geoid
undulation (height of the geoid above the WGS84 ellipsoid, in meters) for a
given latitude/longitude.

The model is shipped as the full 5' grid (see :mod:`egm96.data`). Values are
recovered with bilinear interpolation, which reproduces EGM96 to a maximum
error of about 0.05 m globally.
"""

from __future__ import annotations

from ._model import GRID_RESOLUTION_DEG, MAX_INTERPOLATION_ERROR_M, undulation

__all__ = ["GRID_RESOLUTION_DEG", "MAX_INTERPOLATION_ERROR_M", "undulation"]
