#!/usr/bin/env python3
"""Generate and bake traceable particle-resolved GEOSmie modal optical curves.

Program: generate_modal_optics; Author: EX-aero contributors.
History: 2026-10 initial implementation. Input: pinned GEOSmie + manifest or
committed NetCDF snapshot. Output: snapshot and/or deterministic C++ header.
Usage: --source CHECKOUT --snapshot FILE; --snapshot FILE --header FILE [--check].
Normal builds only bake the committed snapshot; they never run or fetch GEOSmie.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore", category=FutureWarning)

FIELDS = ("qext", "qsca", "bext", "bsca", "bbck", "g", "refreal", "refimag")
RRTMG_BAND_EDGES = (
    (2600,3250),(3250,4000),(4000,4650),(4650,5150),(5150,6150),
    (6150,7700),(7700,8050),(8050,12850),(12850,16000),(16000,22650),
    (22650,29000),(29000,38000),(38000,50000),(820,2600),
    (10,350),(350,500),(500,630),(630,700),(700,820),(820,980),
    (980,1080),(1080,1180),(1180,1390),(1390,1480),(1480,1800),
    (1800,2080),(2080,2250),(2250,2380),(2380,2600),(2600,3250),
)


def interpolate_index(wavelength: np.ndarray, table: np.ndarray, label: str,
                      require_coverage: bool = True) -> np.ndarray:
    """Interpolate a complex index, optionally requiring source coverage."""
    source_wave = table[:, 0] * 1e-6
    order = np.argsort(source_wave)
    source_wave = source_wave[order]
    if require_coverage and (wavelength.min() < source_wave[0] or
                             wavelength.max() > source_wave[-1]):
        raise ValueError(f"FATAL ERROR: {label} wavelength query exceeds source domain")
    real = np.interp(wavelength, source_wave, table[order, 1])
    imag = np.abs(np.interp(wavelength, source_wave, table[order, 2]))
    return real + 1j * imag


def band_coverage(wavelength: np.ndarray, material: np.ndarray,
                  water: np.ndarray) -> bool:
    """A full band's sample grid must fit both refractive-index domains."""
    for source in (material, water):
        axis = source[:, 0] * 1e-6
        if wavelength.min() < axis.min() or wavelength.max() > axis.max():
            return False
    return True


def band_optical_coefficients(radius: xr.DataArray, growth: xr.DataArray,
                              density: float, material: np.ndarray,
                              water: np.ndarray, mie_function: object,
                              samples_per_band: int = 100) -> xr.Dataset:
    """Calculate native band means only where both index sources cover the band."""
    if samples_per_band != 100:
        raise ValueError("FATAL ERROR: pinned RRTMG recipe requires 100 samples")
    shape = (radius.size, growth.size, len(RRTMG_BAND_EDGES))
    values = {key: np.full(shape, np.nan) for key in
              ("band_bext", "band_bsca", "band_bbck", "band_g")}
    available = np.zeros(len(RRTMG_BAND_EDGES), dtype=bool)
    dry_fraction = growth ** -3
    for bi, (lower, upper) in enumerate(RRTMG_BAND_EDGES):
        wavenumber = np.linspace(lower, upper, samples_per_band)
        band_wave = 1e-2 / wavenumber
        if not band_coverage(band_wave, material, water):
            continue
        available[bi] = True
        dry_index = xr.DataArray(interpolate_index(band_wave, material, "material"),
                                 dims="band_sample")
        water_index = xr.DataArray(interpolate_index(band_wave, water, "water"),
                                   dims="band_sample")
        index = dry_fraction * dry_index + (1 - dry_fraction) * water_index
        sample = optical_coefficients(
            radius, growth, xr.DataArray(band_wave, dims="band_sample"),
            index, mie_function)
        cross = np.pi * (radius * growth) ** 2
        dry_mass = 4 * np.pi / 3 * density * radius**3
        bext = sample.qext * cross / dry_mass
        bsca = sample.qsca * cross / dry_mass
        bbck = sample.qb * cross / (4 * np.pi * dry_mass)
        values["band_bext"][:, :, bi] = bext.mean("band_sample").transpose("radius", "rh").to_numpy()
        values["band_bsca"][:, :, bi] = bsca.mean("band_sample").transpose("radius", "rh").to_numpy()
        values["band_bbck"][:, :, bi] = bbck.mean("band_sample").transpose("radius", "rh").to_numpy()
        weighted = (bsca * sample.g).sum("band_sample")
        total = bsca.sum("band_sample")
        values["band_g"][:, :, bi] = xr.where(
            total > 0, weighted / total, 0
        ).transpose("radius", "rh").to_numpy()
    result = xr.Dataset(
        {name: (("radius", "rh", "band"), data)
         for name, data in values.items()}
    )
    result["band_available"] = (("band",), available)
    result = result.assign_coords(band=np.arange(1, 31, dtype=np.int32))
    return result


def optical_coefficients(
    radius: xr.DataArray,
    growth: xr.DataArray,
    wavelength: xr.DataArray,
    refractive: xr.DataArray,
    mie_function: object,
) -> xr.Dataset:
    """Evaluate GEOSmie on eager or lazy coordinate arrays.

    Parameters
    ----------
    radius, growth, wavelength, refractive : xarray.DataArray
        Dry radii [m], wet/dry radius ratio, wavelengths [m], complex indices.
    mie_function : callable
        GEOSmie scalar Mie constructor (or test double).

    Returns
    -------
    xarray.Dataset
        Lazy-preserving efficiencies, asymmetry, and backscatter efficiency.

    Examples
    --------
    See tests/test_modal_optics_generator.py for eager and Dask examples.
    """

    def evaluate(r: float, gf: float, wave: float, index: complex) -> tuple:
        """Return four single-sphere quantities in the GEOSmie convention."""
        mie = mie_function(x=float(2 * np.pi * r * gf / wave), m=complex(index))
        return mie.qext(), mie.qsca(), mie.qb(), mie.asy()

    outputs = xr.apply_ufunc(
        evaluate,
        radius,
        growth,
        wavelength,
        refractive,
        input_core_dims=[[], [], [], []],
        output_core_dims=[[], [], [], []],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float, float, float, float],
    )
    result = xr.Dataset(dict(zip(("qext", "qsca", "qb", "g"), outputs)))
    # Nonabsorbing small-x cases in the pinned solver have bounded relative
    # roundoff (observed <1.6e-6). Correct only that numerical excess; larger
    # violations survive to the artifact validation gate and fail generation.
    result["qsca"] = xr.where(
        (result.qsca > result.qext) & (result.qsca <= result.qext * (1 + 2e-6)),
        result.qext,
        result.qsca,
    )
    result.attrs["history"] = "GEOSmie single-sphere calculation; no PSD integration"
    return result


def generate(source: Path, manifest_path: Path, snapshot: Path) -> None:
    """Generate a snapshot from a clean pinned GEOSmie checkout.

    Parameters
    ----------
    source, manifest_path, snapshot : pathlib.Path
        Checkout, scientific recipe, and output NetCDF path.

    Returns
    -------
    None

    Examples
    --------
    generate(Path("GEOSmie"), Path("manifest.json"), Path("modal.nc"))
    """
    manifest_text = manifest_path.read_text()
    manifest = json.loads(manifest_text)
    band_recipe = manifest.get("band_average", {})
    if band_recipe:
        edges = band_recipe.get("bands", [])
        expected = [
            {"source_id": i + 1, "lower_wavenumber_cm1": lo,
             "upper_wavenumber_cm1": hi}
            for i, (lo, hi) in enumerate(RRTMG_BAND_EDGES)
        ]
        if edges != expected or band_recipe.get("set") != "RRTMG_30":
            raise ValueError("FATAL ERROR: manifest RRTMG metadata differs from canonical bands")
        if band_recipe.get("samples_per_band") != 100 or band_recipe.get("weighting") != "uniform_wavenumber_arithmetic_mean":
            raise ValueError("FATAL ERROR: unsupported RRTMG band averaging recipe")
    revision = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != manifest["geosmie_commit"]:
        raise ValueError("FATAL ERROR: GEOSmie revision does not match the manifest")
    if subprocess.check_output(
        ["git", "-C", str(source), "status", "--porcelain", "--untracked-files=no"],
        text=True,
    ).strip():
        raise ValueError("FATAL ERROR: GEOSmie tracked sources must be clean")
    sys.path.insert(0, str(source / "src" / "pymiecoated"))
    from pymiecoated import Mie

    spec = manifest["radius"]
    radius = xr.DataArray(
        np.geomspace(spec["minimum"], spec["maximum"], spec["count"]), dims="radius"
    )
    rh = xr.DataArray(manifest["rh"], dims="rh")
    wave = xr.DataArray(manifest["wavelength"], dims="wavelength")
    directory = source / "src" / "geosmie" / "data"
    water_path = directory / "refrac.water.txt"
    water = np.loadtxt(water_path, skiprows=13, usecols=[0, 1, 2])
    # Preserve the legacy wavelength-table endpoint behavior. Native-band
    # queries below are stricter and reject any material band outside the data.
    water_index = xr.DataArray(interpolate_index(wave.data, water, "water", False),
                               dims="wavelength")
    records = []
    hashes = {"refrac.water.txt": hashlib.sha256(water_path.read_bytes()).hexdigest()}
    for material in manifest["materials"]:
        path = directory / material["refractive_index"]
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        ri = (
            np.loadtxt(path, usecols=[0, 1, 2])
            if material.get("format") == "wsv"
            else np.loadtxt(path, skiprows=17, usecols=[1, 8, 9])
        )
        dry_index = xr.DataArray(interpolate_index(wave.data, ri, path.name, False),
                                 dims="wavelength")
        growth = (1 + material["kappa"] * rh / (1 - rh)) ** (1 / 3)
        dry_fraction = growth**-3
        index = dry_fraction * dry_index + (1 - dry_fraction) * water_index
        ds = optical_coefficients(radius, growth, wave, index, Mie)
        cross_section = np.pi * (radius * growth) ** 2
        dry_mass = 4 * np.pi / 3 * material["density"] * radius**3
        ds["bext"] = ds.qext * cross_section / dry_mass
        ds["bsca"] = ds.qsca * cross_section / dry_mass
        # GEOSmie qb is the radar backscatter efficiency. Divide by 4*pi
        # to obtain differential 180-degree backscatter per steradian.
        ds["bbck"] = ds.qb * cross_section / (4 * np.pi * dry_mass)
        ds["refreal"], ds["refimag"] = (
            index.real.broadcast_like(ds.qext),
            np.abs(index.imag).broadcast_like(ds.qext),
        )
        ds["rEff"] = radius * growth
        ds["particle_mass"] = dry_mass + 1000 * 4 * np.pi / 3 * radius**3 * (
            growth**3 - 1
        )
        if band_recipe:
            band_result=band_optical_coefficients(
                radius,growth,material["density"],ri,water,Mie,
                band_recipe["samples_per_band"])
            ds=ds.assign_coords(band=band_result.band)
            for name in ("band_bext","band_bsca","band_bbck","band_g",
                         "band_available"):
                ds[name]=band_result[name]
        ds = ds.drop_vars("qb").expand_dims(material=[material["label"]])
        records.append(ds)
    result = xr.concat(records, dim="material").assign_coords(
        radius=radius, rh=rh, wavelength=wave
    )
    result.attrs.update(
        {
            "format": manifest["format"],
            "geosmie_commit": revision,
            "recipe": manifest_text,
            "recipe_sha256": hashlib.sha256(manifest_text.encode()).hexdigest(),
            "input_sha256": json.dumps(hashes, sort_keys=True),
            "citation": "GEOSmie (Kemppinen et al. 2022); OPAC refractive indices; EX-aero modal recipe",
            "history": "Particle-resolved Mie spheres; SI dry-mass normalization; kappa growth once; qsca excess <=2e-6 relative clipped to qext",
            "numpy_version": np.__version__,
            "band_average_recipe": json.dumps(band_recipe, sort_keys=True),
        }
    )
    result.radius.attrs["units"] = "m"
    result.wavelength.attrs["units"] = "m"
    result.rh.attrs["units"] = "fraction"
    # Materialize only at this explicit, bounded artifact-writing boundary.
    validate(result)
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    result.to_netcdf(
        snapshot,
        engine="netcdf4",
        encoding={key: {"zlib": True, "complevel": 6} for key in result.data_vars},
    )


def validate(ds: xr.Dataset) -> None:
    """Validate a bounded generated artifact before publication or baking.

    Parameters
    ----------
    ds : xarray.Dataset
        Generated curves.

    Returns
    -------
    None

    Examples
    --------
    validate(xr.open_dataset("modal.nc"))
    """
    for axis in ("radius", "rh", "wavelength"):
        values = ds[axis].to_numpy()
        if not np.isfinite(values).all() or not (np.diff(values) > 0).all():
            raise ValueError(f"FATAL ERROR: invalid {axis} axis")
    for key in FIELDS:
        values = ds[key].to_numpy()
        if not np.isfinite(values).all() or (key != "g" and (values < 0).any()):
            raise ValueError(f"FATAL ERROR: invalid {key} values")
    if (ds.qsca > ds.qext * (1 + 1e-10)).any().item() or (
        np.abs(ds.g) > 1
    ).any().item():
        ratio = (ds.qsca / ds.qext).max().item()
        asym = np.abs(ds.g).max().item()
        raise ValueError(
            f"FATAL ERROR: Mie optical bounds violated: max qsca/qext={ratio}, max abs(g)={asym}"
        )
    if "band_available" in ds:
        ids = ds.band.to_numpy()
        if not np.array_equal(ids, np.arange(1, 31)):
            raise ValueError("FATAL ERROR: invalid RRTMG band source IDs")
        available = ds.band_available.to_numpy().astype(bool)
        if available.shape != (ds.sizes["material"], 30):
            raise ValueError("FATAL ERROR: invalid per-material band support map")
        for name in ("band_bext", "band_bsca", "band_bbck", "band_g"):
            if name not in ds or ds[name].dims != ("material", "radius", "rh", "band"):
                raise ValueError(f"FATAL ERROR: invalid {name} shape")
            values = ds[name].to_numpy()
            for material in range(ds.sizes["material"]):
                for band in range(30):
                    values_for_band = values[material, :, :, band]
                    if available[material, band]:
                        if not np.isfinite(values_for_band).all():
                            raise ValueError(f"FATAL ERROR: nonfinite supported {name}")
                        if name != "band_g" and (values_for_band < 0).any():
                            raise ValueError(f"FATAL ERROR: negative supported {name}")
                        if name == "band_g" and (np.abs(values_for_band) > 1).any():
                            raise ValueError("FATAL ERROR: band asymmetry outside [-1,1]")
                    elif np.isfinite(values_for_band).any():
                        raise ValueError(f"FATAL ERROR: unsupported {name} must be NaN")
        if (ds.band_bsca.where(ds.band_available) >
                ds.band_bext.where(ds.band_available) * (1 + 1e-10)).any().item():
            raise ValueError("FATAL ERROR: band scattering exceeds extinction")


def bake(
    snapshot: Path, header: Path, check: bool, namespace: str = "mie_modal"
) -> None:
    """Serialize a committed snapshot into a deterministic C++ data header.

    Parameters
    ----------
    snapshot, header : pathlib.Path
        Committed input NetCDF and output header.
    check : bool
        Compare exact bytes without changing files.
    namespace : str
        C++ namespace for this independently generated dataset.

    Returns
    -------
    None

    Examples
    --------
    bake(Path("modal.nc"), Path("ModalOpticsData.hpp"), True)
    """
    if namespace not in ("mie_modal", "mie_mam4xx"):
        raise ValueError("FATAL ERROR: unsupported optical dataset namespace")
    with xr.open_dataset(snapshot) as ds:
        validate(ds)
        manifest = json.loads(ds.attrs["recipe"])
        lines = [
            "// Generated by tools/generate_modal_optics.py; do not edit.",
            "#pragma once",
            f"namespace exaero::{namespace} {{",
            "struct Curve { const char *label; int nr, nh, nl; double density, kappa;",
            "const double *radius, *rh, *wavelength, *rEff, *particle_mass,",
            "*qext, *qsca, *bext, *bsca, *bbck, *g, *refreal, *refimag;",
            "int nb = 0; const double *band_source_ids = nullptr, *band_bext = nullptr,",
            "*band_bsca = nullptr, *band_bbck = nullptr, *band_g = nullptr; };",
            f'inline constexpr char source_version[] = "GEOSmie/{ds.attrs["geosmie_commit"]}";',
            f'inline constexpr char recipe_sha256[] = "{ds.attrs["recipe_sha256"]}";',
            f'inline constexpr char citation[] = "{ds.attrs["citation"]}";',
        ]
        descriptors = []
        for material in manifest["materials"]:
            label = material["label"]
            row = ds.sel(material=label)
            symbols = []
            for key in ("radius", "rh", "wavelength", "rEff", "particle_mass") + FIELDS:
                field = row[key]
                dims = tuple(
                    d for d in ("radius", "rh", "wavelength") if d in field.dims
                )
                values = field.transpose(*dims).to_numpy().ravel()
                symbol = f"{label}_{key}"
                lines.append(f"inline constexpr double {symbol}[] = {{")
                literals = [format(float(v), ".17g") for v in values]
                lines.extend(
                    "  " + ", ".join(literals[i : i + 4]) + ","
                    for i in range(0, len(literals), 4)
                )
                lines.append("};")
                symbols.append(symbol)
            band_symbols = []
            n_band = 0
            if "band_available" in ds:
                supported = row.band_available.to_numpy().astype(bool)
                band_ids = ds.band.to_numpy()[supported]
                n_band = len(band_ids)
                if n_band:
                    for key in ("band_source_ids", "band_bext", "band_bsca",
                                "band_bbck", "band_g"):
                        symbol = f"{label}_{key}"
                        if key == "band_source_ids":
                            values = band_ids
                        else:
                            values = row[key.replace("_source_ids", "")].sel(
                                band=band_ids
                            ).transpose("radius", "rh", "band").to_numpy().ravel()
                        lines.append(f"inline constexpr double {symbol}[] = {{")
                        literals = [format(float(v), ".17g") for v in values]
                        lines.extend(
                            "  " + ", ".join(literals[i : i + 4]) + ","
                            for i in range(0, len(literals), 4)
                        )
                        lines.append("};")
                        band_symbols.append(symbol)
            descriptors.append(
                f'{{"{label}", {ds.sizes["radius"]}, {ds.sizes["rh"]}, {ds.sizes["wavelength"]}, '
                f"{material['density']}, {material['kappa']}, "
                + ", ".join(symbols)
                + (f", {n_band}, " + ", ".join(band_symbols)
                   if n_band else ", 0, nullptr, nullptr, nullptr, nullptr, nullptr")
                + "}"
            )
        lines.append(
            "inline constexpr Curve curves[] = {" + ",\n".join(descriptors) + "};"
        )
        lines.append(f"}} // namespace exaero::{namespace}\n")
        text = "\n".join(lines)
    if check:
        if not header.exists() or header.read_text() != text:
            raise ValueError("FATAL ERROR: generated modal optical header differs")
        print("Modal optical tables deterministic (check mode)")
    elif not header.exists() or header.read_text() != text:
        header.parent.mkdir(parents=True, exist_ok=True)
        header.write_text(text)


def main() -> None:
    """Run the offline generation or deterministic bake CLI.

    Returns
    -------
    None

    Examples
    --------
    python tools/generate_modal_optics.py --help
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).parent / "modal_optics" / "manifest.json",
    )
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--header", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--namespace", choices=["mie_modal", "mie_mam4xx"], default="mie_modal"
    )
    args = parser.parse_args()
    if args.source:
        if args.check:
            parser.error("--check is a read-only snapshot/header comparison")
        generate(args.source, args.manifest, args.snapshot)
    if args.header:
        bake(args.snapshot, args.header, args.check, args.namespace)
    elif not args.source:
        parser.error("provide --source or --header")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
