#!/usr/bin/env python3
"""Build-time generator: pinned GEOSmie RRTMG optics tables -> baked-in C++ float64 header.

Reads the version-tagged snapshot under ``tools/geosmie_snapshot/`` and emits a
deterministic ``src/src_impl/data/generated/MieTableData.hpp`` holding, per species, the
RAW lookup arrays that actually exist in the source (axes + microphysical 2D + spectral 3D)
plus a PROVENANCE block. Derived attributes (growth factor, wet density, per-mass volume /
area, qabs, ssa, lidar ratio) are computed EXACTLY at grid points by MieTableStore from
these raw arrays -- see the file's "Derivations" note -- so this generator stays the single
auditable reader of the pinned source.

Determinism (FR-016, SC-007): the output is a pure function of the committed inputs. Values
are printed with ``%.17g`` (round-trip exact for float64) and every species/variable is
emitted in a fixed order. CI re-runs this and fails on any diff (T039).

Usage:
    python3 tools/generate_mie_tables.py <output_header> [snapshot_dir]
    python3 tools/generate_mie_tables.py --check <output_header> [snapshot_dir]
        ^ determinism gate (T039): regenerate and compare, exit non-zero on any diff.
"""
from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# --- Pinned source mapping (single authoritative version, FR-016) ------------------------
SOURCE_VERSION = "ufs-regtests-input-data-20260617/GOCART/p8c_5d"
CITATION = ("NASA GOCART2G GEOSmie optics tables (Kemppinen et al. 2022; "
            "Colarco/Tsigaridis Chin GOCART); UFS regression-test input-data-20260617")

# Canonical GOCART species -> snapshot file. BR/BRC ships no RRTMG band file in this pinned
# release (PROVENANCE.md), so BR is monochromatic/file-only (FR-017): it is intentionally
# absent from the baked-in default set and reported NotInSource until a file supplies it.
SPECIES_FILES = {
    "DU": "opticsBands_DU.v15_3.RRTMG.nc",
    "SS": "opticsBands_SS.v3_3.RRTMG.nc",
    "SU": "opticsBands_SU.v1_3.RRTMG.nc",
    "BC": "opticsBands_BC.v1_3.RRTMG.nc",
    "OC": "opticsBands_OC.v1_3.RRTMG.nc",
    "NI": "opticsBands_NI.v2_5.RRTMG.nc",
}

# Raw microphysical variables present in the band tables: (name, dims-rank).
MICRO_1D = ["rLow", "rUp"]          # (radius,)
MICRO_2D = ["rEff", "rMass"]        # (radius, rh)
# Raw spectral variables: (radius, rh, lambda).
SPECTRAL_3D = ["qext", "qsca", "bext", "bsca", "bbck", "g", "refreal", "refimag"]

# Canonical attribute units (data-model.md) recorded verbatim for audit.
UNITS = {
    "radius": "m", "rh": "fraction", "lambda": "dimensionless(band_index)",
    "rLow": "m", "rUp": "m", "rEff": "m", "rMass": "kg",
    "qext": "1", "qsca": "1", "bext": "m2 (kg dry mass)-1",
    "bsca": "m2 (kg dry mass)-1", "bbck": "m2 (kg dry mass)-1 sr-1",
    "g": "1", "refreal": "1", "refimag": "1",
}


def _fmt(v: float) -> str:
    """Round-trip-exact float64 literal (deterministic across platforms)."""
    return f"{float(v):.17g}"


def _emit_array(out: list[str], cname: str, length: int, flat: np.ndarray) -> None:
    """Emit one flat float64 C array in fixed row-major (radius, rh, lambda) order.

    The consumer treats these as layout_left Fortran-order fields via the store's explicit
    strides; here we emit in the NetCDF read order so the header is a faithful copy.
    """
    vals = np.asarray(flat, dtype=np.float64).ravel(order="C")
    assert vals.size == length, f"{cname}: expected {length} values, got {vals.size}"
    out.append(f"inline constexpr std::size_t {cname}_N = {length};")
    out.append(f"inline constexpr double {cname}[{length}] = {{")
    line: list[str] = []
    for x in vals:
        line.append(_fmt(x))
        if len(line) == 4:
            out.append("    " + ", ".join(line) + ",")
            line = []
    if line:
        out.append("    " + ", ".join(line) + ",")
    out.append("};")
    out.append("")


def _load_species(snapshot: Path, label: str, fname: str) -> dict:
    """Read one species band table into raw float64 arrays + axes."""
    path = snapshot / fname
    if not path.is_file():
        raise FileNotFoundError(f"FATAL ERROR: pinned snapshot missing: {path}")
    ds = Dataset(path)
    try:
        n_rad = len(ds.dimensions["radius"])
        n_rh = len(ds.dimensions["rh"])
        n_lam = len(ds.dimensions["lambda"])

        def read(var: str) -> np.ndarray:
            # np.asarray(var[:]) avoids the masked-array dtype kwarg (netCDF4 >=1.7
            # __array__ takes no positional args); filled gives exact float64 bits.
            return np.asarray(ds.variables[var][:], dtype=np.float64)

        data = {
            "n_radius": n_rad, "n_rh": n_rh, "n_lambda": n_lam,
            "radius": read("radius"),
            "rh": read("rh"),
            "lambda": read("lambda"),
        }
        for v in MICRO_1D:
            data[v] = read(v)
        for v in MICRO_2D:
            data[v] = read(v)
        for v in SPECTRAL_3D:
            data[v] = read(v)
    finally:
        ds.close()
    _validate(label, data)
    return data


def _validate(label: str, d: dict) -> None:
    """Load-time physical-bounds gate (FR-009, data-model Validation Rules). Aborts loudly.

    The pinned source is float32; widened to float64 a nominal 0.99 is 0.9900000095..., so
    bound checks use a float32-appropriate slack (FR-005 guarantees exact reproduction of the
    stored float32 bits, not of the decimal literal).
    """
    f32 = 1e-6  # slack for float32-sourced grid values near 0 / 1 / 0.99

    def fail(msg: str) -> None:
        raise ValueError(f"FATAL ERROR: {label}: {msg}")

    if np.any(~np.isfinite(np.concatenate([d["radius"], d["rh"], d["lambda"]]))):
        fail("non-finite axis value")
    if np.any(d["radius"] <= 0):
        fail("radius entries must be > 0")
    # Radius may be non-STRICTLY monotone: multi-mode species (BC, OC) carry duplicate
    # representative radii (hydrophobic + hydrophilic share a size), so radius is an index
    # axis there, not a curve. The store disables radius resampling for non-strictly-increasing
    # axes and surfaces an explicit status (R9); the generator only rejects decreasing order.
    if np.any(np.diff(d["radius"]) < 0):
        fail("radius axis not monotonically non-decreasing")
    if np.any(np.diff(d["rh"]) <= 0):
        fail("rh axis not strictly increasing")
    if d["rh"][0] < -f32 or d["rh"][-1] > 0.99 + f32:
        fail(f"rh out of [0,0.99]: [{d['rh'][0]},{d['rh'][-1]}]")
    if np.any(d["rEff"] <= 0):
        fail("rEff must be > 0")
    if np.any(d["rMass"] <= 0):
        fail("rMass must be > 0")
    qext, qsca, g = d["qext"], d["qsca"], d["g"]
    if np.any(qext < -f32) or np.any(qsca < -f32):
        fail("negative efficiency")
    if np.any(qsca > qext + 1e-5):
        fail("qsca > qext")
    ssa = np.where(qext > 1e-30, qsca / np.where(qext <= 1e-30, 1.0, qext), 0.0)
    if np.any(ssa < -f32) or np.any(ssa > 1.0 + f32):
        fail("derived ssa outside [0,1]")
    if np.any(np.abs(g) > 1.0 + f32):
        fail("asymmetry g outside [-1,1]")


def generate(snapshot: Path) -> str:
    """Return the full generated header text (deterministic)."""
    species = {lbl: _load_species(snapshot, lbl, fn) for lbl, fn in SPECIES_FILES.items()}

    out: list[str] = []
    out.append("// AUTO-GENERATED by tools/generate_mie_tables.py -- DO NOT EDIT.")
    out.append("// Source: pinned GEOSmie RRTMG band tables (see PROVENANCE block below).")
    out.append("//")
    out.append("// Baked-in float64 lookup data for the EX-aero GEOSmie MIE attribute surface")
    out.append("// (specs/001-geosmie-lut-attributes). Header-only, no Kokkos, no runtime I/O.")
    out.append("//")
    out.append("// Derivations (computed by MieTableStore, exact at grid points):")
    out.append("//   growth_factor(bin,rh)   = rEff(bin,rh) / rEff(bin,0)")
    out.append("//   wet_particle_density    = rMass / (4/3 pi rEff^3)")
    out.append("//   volume_per_mass         = (4/3 pi rEff^3) / rMass")
    out.append("//   area_per_mass           = pi rEff^2 / rMass")
    out.append("//   qabs = qext-qsca ; ssa = qsca/qext ; lidar_ratio = bext/bbck (guarded)")
    out.append("#pragma once")
    out.append("#include <cstddef>")
    out.append("")
    out.append("namespace exaero {")
    out.append("namespace mie_baked {")
    out.append("")
    out.append("// --- PROVENANCE (FR-004, FR-016) ---")
    out.append(f'inline constexpr const char* SOURCE_VERSION = "{SOURCE_VERSION}";')
    out.append(f'inline constexpr const char* CITATION = "{CITATION}";')
    out.append("inline constexpr int NUM_SPECIES = %d;" % len(species))
    out.append("inline constexpr const char* SPECIES_LABELS[NUM_SPECIES] = {%s};" %
               ", ".join('"%s"' % s for s in species))
    out.append("")

    # Per-species axis lengths as named constants (data, not hardcoded in consumers).
    out.append("struct SpeciesDims { int n_radius; int n_rh; int n_lambda; };")
    out.append("inline constexpr SpeciesDims SPECIES_DIMS[NUM_SPECIES] = {")
    for lbl, d in species.items():
        out.append("    {%d, %d, %d}, // %s" %
                   (d["n_radius"], d["n_rh"], d["n_lambda"], lbl))
    out.append("};")
    out.append("")

    for lbl, d in species.items():
        nr, nh, nl = d["n_radius"], d["n_rh"], d["n_lambda"]
        out.append(f"// ===== species {lbl}: radius={nr} rh={nh} lambda={nl} =====")
        _emit_array(out, f"{lbl}_radius", nr, d["radius"])
        _emit_array(out, f"{lbl}_rh", nh, d["rh"])
        _emit_array(out, f"{lbl}_lambda", nl, d["lambda"])
        for v in MICRO_1D:
            _emit_array(out, f"{lbl}_{v}", nr, d[v])
        for v in MICRO_2D:
            _emit_array(out, f"{lbl}_{v}", nr * nh, d[v])
        for v in SPECTRAL_3D:
            _emit_array(out, f"{lbl}_{v}", nr * nh * nl, d[v])

    out.append("} // namespace mie_baked")
    out.append("} // namespace exaero")
    out.append("")
    return "\n".join(out)


def main(argv: list[str]) -> int:
    args = [a for a in argv[1:] if a != "--check"]
    check = "--check" in argv[1:]
    if len(args) < 1:
        sys.stderr.write("usage: generate_mie_tables.py [--check] <output_header> "
                         "[snapshot_dir]\n")
        return 2
    out_path = Path(args[0])
    default_snapshot = Path(__file__).resolve().parent / "geosmie_snapshot"
    snapshot = Path(args[1]) if len(args) > 1 else default_snapshot

    text = generate(snapshot)
    if check:
        # CI determinism gate (T039, FR-016/SC-007): regenerate in memory and compare
        # against the committed header; any diff fails the gate. Never writes.
        if not out_path.is_file():
            sys.stderr.write(f"FATAL ERROR: {out_path} missing; regenerate with "
                             f"tools/generate_mie_tables.py\n")
            return 1
        if out_path.read_text() != text:
            digest_new = hashlib.sha256(text.encode()).hexdigest()[:12]
            digest_old = hashlib.sha256(out_path.read_bytes()).hexdigest()[:12]
            sys.stderr.write(f"FATAL ERROR: {out_path} is stale (sha256:{digest_old}, "
                             f"expected sha256:{digest_new}); run "
                             f"python3 tools/generate_mie_tables.py {out_path}\n")
            return 1
        sys.stderr.write("MieTableData.hpp is deterministic (check mode, no diff).\n")
        return 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Only rewrite when content changed -> stable timestamps, CI-diffable (T039).
    if out_path.is_file() and out_path.read_text() == text:
        sys.stderr.write("MieTableData.hpp up to date (no diff).\n")
        return 0
    out_path.write_text(text)
    digest = hashlib.sha256(text.encode()).hexdigest()[:12]
    sys.stderr.write(f"Wrote {out_path} (sha256:{digest}, {len(text)} bytes)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
