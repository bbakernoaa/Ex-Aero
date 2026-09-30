#!/usr/bin/env python3
"""Deterministic generator for the GEOSmie MIE runtime-file fixtures (tests/data/).

Part of specs/001-geosmie-lut-attributes (US2 T029, US3 T031/T034, polish T036/T037).
All valid fixtures are derived from the pinned snapshot under tools/geosmie_snapshot/
at full float64 precision so file-parsed values equal the source bits exactly
(FR-005 grid tolerance 1e-7). The corrupt fixtures are hand-authored (see README)
and intentionally invalid - this script writes only the valid ones.

Schema: the native portable text format parsed by MieFileLoader
(`# @format exaero_mie_v1`, `# @source_version`, `# @citation`, `# @species`,
`# @axis <name> <units> <n> [vals...]`, `# @field <name> <dims-csv> <units> [count]`
followed by whitespace value rows). The loader derives each field's element count from
the declared axis lengths; the trailing count is readability only. Field names are the
store's canonical query names so a merged file field is actually surfaced by
queryAttribute/computeAttributes.

Usage:  python3 tools/make_test_fixtures.py
"""
import glob
import os
import sys
import numpy as np
import netCDF4 as nc

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAP = os.path.join(REPO, "tools", "geosmie_snapshot")
OUT = os.path.join(REPO, "tests", "data")

SOURCE_VERSION = "ufs-regtests-input-data-20260617/GOCART/p8c_5d"
CITATION = ("NASA GOCART2G GEOSmie optics tables (Kemppinen et al. 2022; "
            "Colarco/Tsigaridis Chin GOCART); UFS regression-test input-data-20260617")

UNITS = {
    "radius": "m", "rh": "fraction", "lambda": "m", "pol": "1", "moment": "1",
    "qext": "1", "qsca": "1", "bext": "m^2 (kg dry mass)^-1",
    "bsca": "m^2 (kg dry mass)^-1", "bbck": "m^2 (kg dry mass)^-1 sr^-1",
    "g": "1", "refreal": "1", "refimag": "1",
    "wet_particle_density": "kg m^-3", "growth_factor": "1", "pmom": "1",
}
DIMS = {
    "qext": "radius,rh,lambda", "qsca": "radius,rh,lambda",
    "bext": "radius,rh,lambda", "bsca": "radius,rh,lambda",
    "bbck": "radius,rh,lambda", "g": "radius,rh,lambda",
    "refreal": "radius,rh,lambda", "refimag": "radius,rh,lambda",
    "wet_particle_density": "radius,rh", "growth_factor": "radius,rh",
    "pmom": "radius,rh,lambda,pol,moment",
}


def fmt(v):
    return repr(float(v))  # shortest round-trip float64


def emit_axis(lines, name, vals):
    vals = np.asarray(vals, dtype=np.float64)
    lines.append(f"# @axis {name} {UNITS[name]} {vals.size} " +
                 " ".join(fmt(x) for x in vals))


def emit_field(lines, name, arr):
    arr = np.asarray(arr, dtype=np.float64).ravel(order="C")
    lines.append(f"# @field {name} {DIMS[name]} {UNITS[name]} {arr.size}")
    for i in range(0, arr.size, 8):
        lines.append(" ".join(fmt(x) for x in arr[i:i + 8]))


def load_species(lbl):
    hits = glob.glob(os.path.join(SNAP, f"opticsBands_{lbl}.*.RRTMG.nc"))
    assert len(hits) == 1, hits
    ds = nc.Dataset(hits[0])
    d = {k: np.asarray(ds.variables[k][:], dtype=np.float64)
         for k in ds.variables if k != "qname"}
    ds.close()
    return d


def header(species):
    return ["# @format exaero_mie_v1",
            f"# @source_version {SOURCE_VERSION}",
            f"# @citation {CITATION}",
            f"# @species {species}"]


def fixture_dust_monochromatic():
    """DU: replace the 30-band RRTMG lambda axis with a 3-wavelength monochromatic set
    (0.355, 0.55, 1.33 um) and add polarized moments (a file-only category). radius/rh
    are inherited from baked (not emitted). Proves file values EXTEND the baked set on
    the new axis (C7) and that file-only pmom surfaces (US3)."""
    d = load_species("DU")
    nR, nH = d["radius"].size, d["rh"].size
    mono_lam = [0.355e-6, 0.55e-6, 1.33e-6]
    nL = len(mono_lam)
    band2 = 2  # ~visible RRTMG band reused as the 0.55 um value

    def mono(field):
        base = d[field][:, :, band2]        # (radius, rh)
        out = np.empty((nR, nH, nL))
        out[:, :, 1] = base
        out[:, :, 0] = base * 1.30          # 0.355 um: higher (shorter wavelength)
        out[:, :, 2] = base * 0.60          # 1.33 um: lower
        if field == "g":                     # asymmetry is bounded [-1,1] (SC-005):
            out = np.clip(out, -0.999, 0.999)  # keep the synthetic set physically valid
        return out

    lines = header("DU")
    # Files are self-describing: the loader derives each field's element count from the
    # axes declared IN THE FILE, so radius/rh (identical to baked) are emitted too. The
    # store merge keeps baked radius/rh and replaces only the spectral axis + fields.
    emit_axis(lines, "radius", d["radius"])
    emit_axis(lines, "rh", d["rh"])
    emit_axis(lines, "lambda", mono_lam)
    for f in ("qext", "qsca", "bext", "bsca", "bbck", "g", "refreal", "refimag"):
        emit_field(lines, f, mono(f))

    # Polarized moments: rank-5 (radius, rh, lambda, pol, moment), 6 elements
    # (P11,P12,P33,P34,P22,P44) x 3 moments. value = qext(bin,rh) * 0.5^moment *
    # (1 + 0.1*element): deterministic, ordering- and count-testable.
    n_pol, n_mom = 6, 3
    qe = d["qext"][:, :, band2]             # (radius, rh)
    pm = np.empty((nR, nH, nL, n_pol, n_mom))
    for p in range(n_pol):
        for m in range(n_mom):
            pm[:, :, :, p, m] = (qe[:, :, None] * (0.5 ** m) * (1.0 + 0.1 * p)) \
                .repeat(nL, axis=2)
    # pol element labels are the documented ordering P11,P12,P33,P34,P22,P44 -> 1..6
    # (data-model); moment labels 0..M-1. Axis VALUES are labels; the store addresses
    # them positionally (0-based) and treats the count as data.
    emit_axis(lines, "pol", [1, 2, 3, 4, 5, 6])
    emit_axis(lines, "moment", [0, 1, 2])
    emit_field(lines, "pmom", pm)
    return "\n".join(lines) + "\n", dict(nR=nR, nH=nH, nL=nL, n_pol=n_pol,
                                         n_mom=n_mom, mono_lam=mono_lam,
                                         qe00=float(qe[0, 10]))


def fixture_dust_override():
    """DU: override bext on the SAME 30-band axis with a distinct constant field so a
    query reports delivery_source=runtime-file (C7) and differs from baked."""
    d = load_species("DU")
    nR, nH, nL = d["radius"].size, d["rh"].size, d["lambda"].size
    lines = header("DU")
    emit_axis(lines, "radius", d["radius"])
    emit_axis(lines, "rh", d["rh"])
    emit_axis(lines, "lambda", d["lambda"])
    emit_field(lines, "bext", np.full((nR, nH, nL), 1234.5))
    return "\n".join(lines) + "\n"


def fixture_nitrate_density():
    """NI: supply wet_particle_density + growth_factor as explicit file fields (rank-2
    radius,rh) so the store reports them with delivery_source=runtime-file. The pinned
    RRTMG band tables carry no source density/growth var - these prove config/file-
    supplied microphysics are surfaced (FR-001)."""
    d = load_species("NI")
    rEff, rMass = d["rEff"], d["rMass"]
    vol = (4.0 / 3.0) * np.pi * rEff ** 3
    wd = np.where(vol > 0, rMass / vol, 0.0)
    gf = np.where(rEff[:, :1] > 0, rEff / rEff[:, :1], 1.0)
    lines = header("NI")
    emit_axis(lines, "radius", d["radius"])
    emit_axis(lines, "rh", d["rh"])
    emit_field(lines, "wet_particle_density", wd)
    emit_field(lines, "growth_factor", gf)
    return "\n".join(lines) + "\n"


def main():
    os.makedirs(OUT, exist_ok=True)
    dust_mono, meta = fixture_dust_monochromatic()
    with open(os.path.join(OUT, "mie_dust_monochromatic.txt"), "w") as f:
        f.write(dust_mono)
    with open(os.path.join(OUT, "mie_dust_override.txt"), "w") as f:
        f.write(fixture_dust_override())
    with open(os.path.join(OUT, "mie_nitrate_density.txt"), "w") as f:
        f.write(fixture_nitrate_density())

    # Deliberately-malformed fixtures for the fail-fast gate. Each is
    # valid up to exactly one fault so the loader aborts on a specific check.
    with open(os.path.join(OUT, "mie_corrupt_schema.txt"), "w") as fh:
        fh.write(corrupt_schema())
    with open(os.path.join(OUT, "mie_corrupt_unit.txt"), "w") as fh:
        fh.write(corrupt_unit())
    with open(os.path.join(OUT, "mie_corrupt_version.txt"), "w") as fh:
        fh.write(corrupt_version())
    with open(os.path.join(OUT, "mie_truncated.txt"), "w") as fh:
        fh.write(corrupt_truncated())
    with open(os.path.join(OUT, "mie_dust7.txt"), "w") as fh:
        fh.write(fixture_dust7())
    print("wrote fixtures to", OUT)
    print("dust mono meta:", meta)


def fixture_dust7():
    """A NEW species label 'DUST7' (not in the baked set) declaring a PRIME bin count
    (7 radius nodes) and a NON-SOURCE RH length (7 levels) over 3 bands. Because the
    label is absent from baked, the file fully defines the curve -> proves no bin/RH/band
    count is hardcoded on any path (C12/R10). Deterministic synthetic values."""
    nR, nH, nL = 7, 7, 3
    radius = [1.0e-7 * (i + 1) for i in range(nR)]     # 1e-7 .. 7e-7 m (strictly inc.)
    rh = [i / (nH - 1) * 0.99 for i in range(nH)]      # 0 .. 0.99 (7 levels)
    lam = [1.0, 2.0, 3.0]
    lines = header("DUST7")
    emit_axis(lines, "radius", np.asarray(radius))
    emit_axis(lines, "rh", np.asarray(rh))
    emit_axis(lines, "lambda", np.asarray(lam))
    # qext[b,h,l] = (b+1) + 0.1*h + 0.01*l -> distinct, monotone, finite everywhere.
    qe = np.array([[[ (b + 1) + 0.1 * h + 0.01 * l for l in range(nL)]
                    for h in range(nH)] for b in range(nR)])
    emit_field(lines, "qext", qe)
    emit_field(lines, "qsca", qe * 0.5)
    # Mass cross-sections + asymmetry so the derived-field pass (qabs/ssa/lidar) and the
    # device spectral-block wiring have complete raw inputs on this prime-count axis.
    emit_field(lines, "bext", qe * 2.0e-11)
    emit_field(lines, "bsca", qe * 1.0e-11)
    emit_field(lines, "bbck", qe * 1.0e-13)
    emit_field(lines, "g", np.full(qe.shape, 0.6))
    return "\n".join(lines) + "\n"


# --- corrupt-file fixtures each carries exactly one schema/unit/version/
# truncation fault; the loader must abort with "FATAL ERROR:" and no silent fallback. ---

def corrupt_schema():
    """Unknown directive -> loader_fatal('unknown directive ...')."""
    return "\n".join(header("DU") + [
        "# @axis radius m 1 1.0e-06",
        "# @axis rh fraction 1 0.0",
        "# @bogus_this_is_not_a_directive oops",
    ]) + "\n"


def corrupt_unit():
    """@field declares a unit disagreeing with the canonical one -> unit-mismatch abort."""
    lines = header("DU")
    emit_axis(lines, "radius", [1.0e-06])
    emit_axis(lines, "rh", [0.0])
    # qext canonical unit is "1"; declare it in metres to trip the check.
    lines.append("# @field qext radius,rh m 1")
    lines.append("1.0")
    return "\n".join(lines) + "\n"


def corrupt_version():
    """Missing @source_version -> loader_fatal('missing @source_version (FR-016)')."""
    return "\n".join([
        "# @format exaero_mie_v1",
        "# @citation no version line here",
        "# @species DU",
        "# @axis radius m 1 1.0e-06",
        "# @axis rh fraction 1 0.0",
        "# @field qext radius,rh 1 1",
        "1.0",
    ]) + "\n"


def corrupt_truncated():
    """@field supplies fewer values than the declared axes require -> truncation abort."""
    lines = header("DU")
    emit_axis(lines, "radius", [1.0e-06, 2.0e-06])
    emit_axis(lines, "rh", [0.0, 0.5])
    # qext over (radius=2, rh=2) needs 4 values; supply only 3.
    lines.append("# @field qext radius,rh 1 4")
    lines.append("1.0 2.0 3.0")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
