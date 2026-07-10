"""Tests for the public :mod:`egm96` API.

The accuracy tests compare against the full 5' EGM96 model via ``pygeodesy``
(the ``dev`` dependency group). They are skipped automatically if the reference
PGM is not available locally.
"""

from __future__ import annotations

import pathlib

import numpy as np
import numpy.typing as npt
import pytest

import egm96

#: Canonical location of the reference PGM, written by ``./generate-egm96-grid``.
REFERENCE_PGM = pathlib.Path(__file__).parent.parent / "geoids" / "egm96-5.pgm"

needs_reference = pytest.mark.skipif(
    not REFERENCE_PGM.exists(),
    reason="reference egm96-5.pgm not available (run ./generate-egm96-grid first)",
)


def test_scalar_returns_float() -> None:
    n = egm96.undulation(48.14, 11.58)
    assert isinstance(n, float)
    assert 45.4 < n < 55.6  # southern Germany is well-known ~+45.5 m


def test_array_broadcasts() -> None:
    lats = np.array([48.14, 0.0, -33.87])
    lons = np.array([11.58, 0.0, 151.21])
    out = egm96.undulation(lats, lons)
    assert isinstance(out, np.ndarray)
    assert out.shape == (3,)


def test_float_broadcasts() -> None:
    lats = np.array([48.14, 0.0, -33.87])
    lons = 11.58
    out = egm96.undulation(lats, lons)
    assert isinstance(out, np.ndarray)
    assert out.shape == (3,)
    assert np.array_equal(out, egm96.undulation(lats, np.full_like(lats, lons)))


def test_longitude_wraps() -> None:
    a = egm96.undulation(10.0, 179.9)
    b = egm96.undulation(10.0, 179.9 - 360.0)  # same meridian, wrapped
    assert a == pytest.approx(b, abs=1e-9)


def test_poles_do_not_raise() -> None:
    assert np.isfinite(egm96.undulation(90.0, 0.0))
    assert np.isfinite(egm96.undulation(-90.0, 123.0))


def _reference_undulations(lats: npt.ArrayLike, lons: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Evaluate the full 5' EGM96 model at the given coordinates."""
    from pygeodesy.ellipsoidalKarney import LatLon
    from pygeodesy.geoids import GeoidKarney

    geoid = GeoidKarney(str(REFERENCE_PGM))
    lats = np.asarray(lats, dtype=float)
    lons = np.asarray(lons, dtype=float)
    return np.fromiter(
        (geoid(LatLon(float(a), float(o))) for a, o in zip(lats, lons, strict=True)),
        dtype=float,
        count=lats.size,
    )


@needs_reference
def test_accuracy_matches_full_model() -> None:
    """Random global points must track the full model to within the bound."""
    pytest.importorskip("pygeodesy", reason="install the 'dev' dependency group")
    rng = np.random.default_rng(42)
    lats = rng.uniform(-89.5, 89.5, 3000)
    lons = rng.uniform(-179.99, 179.99, 3000)

    err = np.abs(egm96.undulation(lats, lons) - _reference_undulations(lats, lons))
    assert err.max() < egm96.MAX_INTERPOLATION_ERROR_M, f"max error {err.max():.3f} m"


@needs_reference
def test_accuracy_across_global_grid() -> None:
    """A dense sweep of the whole globe -- including the poles and the
    antimeridian -- must stay far below the 1 m requirement.
    """
    pytest.importorskip("pygeodesy", reason="install the 'dev' dependency group")
    # Sweep on a 2.6deg mesh shifted by half a 5' cell so points fall *between*
    # grid nodes, where bilinear error peaks (2.6deg is not a multiple of the 5'
    # spacing, so the samples do not silently land on nodes). Add the exact
    # edges to cover the poles and both ends of the antimeridian.
    half_cell = egm96.GRID_RESOLUTION_DEG / 2
    lat_axis = np.union1d(np.arange(-90.0, 89.0, 2.6) + half_cell, [-90.0, 0.0, 90.0])
    lon_axis = np.union1d(np.arange(-180.0, 179.0, 2.6) + half_cell, [-180.0, 0.0, 180.0])
    lat_grid, lon_grid = np.meshgrid(lat_axis, lon_axis, indexing="ij")
    lats = lat_grid.ravel()
    lons = lon_grid.ravel()

    err = np.abs(egm96.undulation(lats, lons) - _reference_undulations(lats, lons))
    worst = int(np.argmax(err))
    assert err.max() < 0.25, (
        f"max error {err.max():.3f} m at lat={lats[worst]:.2f}, lon={lons[worst]:.2f} "
        f"(mean {err.mean():.4f} m over {err.size} points)"
    )
