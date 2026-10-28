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


def test_native_rRTMG_band_average_and_extreme_lw_support() -> None:
    """Band mean uses 100 uniform-wavenumber samples and never extrapolates."""
    radius = xr.DataArray([1e-7], dims="radius")
    growth = xr.DataArray([1.0], dims="rh")
    refractive = np.array([[0.1, 1.5, 0.01], [200.0, 1.5, 0.01]])
    result = GENERATOR.band_optical_coefficients(
        radius, growth, 1800.0, refractive, refractive, Sphere
    )
    assert result.band_available.sel(band=3).item()
    assert not result.band_available.sel(band=15).item()
    assert np.isnan(result.band_bext.sel(band=15)).all()
    wavenumber = np.linspace(4000, 4650, 100)
    wavelength = 1e-2 / wavenumber
    x = 2 * np.pi * 1e-7 / wavelength
    expected = np.mean(2 * x * np.pi * (1e-7) ** 2 /
                       (4 * np.pi / 3 * 1800 * (1e-7) ** 3))
    assert result.band_bext.sel(band=3).item() == pytest.approx(expected)
    assert result.band_g.sel(band=3).item() == pytest.approx(0.4)


def test_modal_band_metadata_is_pinned_and_wavelength_interpolation_is_bounded() -> None:
    recipe = json.loads((ROOT / "tools/modal_optics/manifest.json").read_text())[
        "band_average"
    ]
    assert recipe["set"] == "RRTMG_30"
    assert recipe["definition_commit"] == "fbc507e6edc09cf6dcaa5a1d88d63ffda4f30652"
    assert len(recipe["bands"]) == 30
    assert recipe["bands"][12]["lower_wavenumber_cm1"] == 38000
    source = np.array([[0.25, 1.5, .01], [40., 1.6, .02]])
    assert not GENERATOR.band_coverage(np.array([0.2e-6, 0.25e-6]), source,
                                       np.array([[.2, 1.3, .01], [30000, 1.2, .01]]))
    with pytest.raises(ValueError, match="exceeds source domain"):
        GENERATOR.interpolate_index(np.array([0.2e-6]), source, "dust")


@pytest.mark.parametrize("snapshot", ["modal_mie_v1.nc", "mam4xx_mie_v1.nc"])
def test_generated_snapshot_rayleigh_limit_and_mass_normalization(
    snapshot: str,
) -> None:
    """Check the Mie small-sphere limit independently of the GEOSmie engine."""
    with xr.open_dataset(ROOT / "tools/modal_optics" / snapshot) as ds:
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


def test_mam4xx_recipe_and_header() -> None:
    """Verify native material properties, surrogate labels and deterministic bake."""
    recipe = json.loads((ROOT / "tools/modal_optics/mam4xx_manifest.json").read_text())
    assert [m["density"] for m in recipe["materials"]] == [
        1000,
        1770,
        1000,
        1700,
        1900,
        2600,
        1601,
    ]
    assert [m["kappa"] for m in recipe["materials"]] == [
        0.1,
        0.507,
        1e-10,
        1e-10,
        1.16,
        0.14,
        0.1,
    ]
    for m in recipe["materials"]:
        if m["label"] in ("MAM4_SOA", "MAM4_POM", "MAM4_MOM"):
            assert m["surrogate"] and "surrogate" in m["description"]
    with xr.open_dataset(ROOT / "tools/modal_optics/mam4xx_mie_v1.nc") as ds:
        snapshot_recipe = json.loads(ds.attrs["recipe"])
        # The checked-in MAM4xx snapshot predates native bands. Keep asserting
        # its exact legacy recipe while allowing a regenerated snapshot to
        # carry the newly pinned band recipe.
        assert {k: v for k, v in recipe.items() if k != "band_average"} == snapshot_recipe
        if "band_average" in snapshot_recipe:
            assert snapshot_recipe == recipe
            assert ds.band.size == 30
        else:
            assert "band" not in ds.coords
    GENERATOR.bake(
        ROOT / "tools/modal_optics/mam4xx_mie_v1.nc",
        ROOT / "src/src_impl/data/generated/Mam4xxOpticsData.hpp",
        True,
        "mie_mam4xx",
    )
