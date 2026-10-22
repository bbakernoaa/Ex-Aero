"""Program: modal table validation; Author: EX-aero contributors.

History: 2026-10. Usage: pytest tests/test_modal_optics_generator.py.
Input: committed GEOSmie snapshot; output: test assertions only.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "modal_generator", ROOT / "tools/generate_modal_optics.py"
)
GENERATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GENERATOR)


class Sphere:
    """Scalar Mie test double for exercising eager/lazy dataflow."""

    def __init__(self, x: float, m: complex) -> None:
        self.x = x

    def qext(self) -> float:
        return self.x * 2

    def qsca(self) -> float:
        return self.x

    def qb(self) -> float:
        return self.x / 10

    def asy(self) -> float:
        return 0.4


@pytest.mark.parametrize("lazy", [False, True])
def test_eager_and_lazy_optical_generation(lazy: bool) -> None:
    """Verify identical coordinates/results without forcing lazy inputs eager."""
    radius = xr.DataArray([1e-8, 1e-7], dims="radius")
    growth = xr.DataArray([1.0, 1.5], dims="rh")
    wave = xr.DataArray([5.5e-7, 8.7e-7], dims="wavelength")
    index = xr.DataArray([1.5 + 0.01j, 1.6 + 0.02j], dims="wavelength")
    eager = GENERATOR.optical_coefficients(radius, growth, wave, index, Sphere)
    if lazy:
        radius, growth, wave, index = (v.chunk() for v in (radius, growth, wave, index))
    result = GENERATOR.optical_coefficients(radius, growth, wave, index, Sphere)
    if lazy:
        assert result.qext.chunks is not None
    xr.testing.assert_allclose(result.compute(), eager)


def test_generated_snapshot_rayleigh_limit_and_mass_normalization() -> None:
    """Check the Mie small-sphere limit independently of the GEOSmie engine."""
    with xr.open_dataset(ROOT / "tools/modal_optics/modal_mie_v1.nc") as ds:
        GENERATOR.validate(ds)
        recipe = json.loads(ds.attrs["recipe"])
        assert ds.attrs["geosmie_commit"] == recipe["geosmie_commit"]
        for material in recipe["materials"]:
            row = ds.sel(material=material["label"])
            dry = row.isel(radius=slice(0, 20), rh=0)
            index = dry.refreal + 1j * dry.refimag
            x = 2 * np.pi * dry.radius / dry.wavelength
            rayleigh = 8 / 3 * x**4 * np.abs((index**2 - 1) / (index**2 + 2)) ** 2
            np.testing.assert_allclose(dry.qsca, rayleigh, rtol=0.003, atol=1e-30)
            gf = (1 + material["kappa"] * row.rh / (1 - row.rh)) ** (1 / 3)
            expected = (
                row.qext
                * np.pi
                * (row.radius * gf) ** 2
                / (4 * np.pi / 3 * material["density"] * row.radius**3)
            )
            xr.testing.assert_allclose(row.bext, expected.transpose(*row.bext.dims))


def test_committed_header_is_deterministic() -> None:
    """Regenerate text only and verify exact bytes without modifying the repo."""
    GENERATOR.bake(
        ROOT / "tools/modal_optics/modal_mie_v1.nc",
        ROOT / "src/src_impl/data/generated/ModalOpticsData.hpp",
        True,
    )
