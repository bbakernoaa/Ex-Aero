# Generated modal optical snapshot

`modal_mie_v1.nc` is produced from the pinned `bbakernoaa/GEOSmie` Mie engine
and its refractive-index input files using `manifest.json`. It contains spherical
single-particle optical curves for sulfate, soot, water-soluble organic material,
sea salt, mineral dust, and ammonium nitrate. Kappa growth and volume-weighted
mixing with water are explicit EX-aero recipe choices. They do not reproduce the
legacy GOCART size distributions or humidity-growth prescriptions.

The NetCDF attributes record the GEOSmie revision, exact recipe and SHA-256,
input-file hashes, numerical convention, normalization, and NumPy version.
The generated C++ header preserves float64 values exactly. Normal builds bake
this committed snapshot and require no GEOSmie checkout or network.

Offline generation (use an isolated Python environment with NumPy, SciPy,
Numba, netCDF4, and xarray):

```sh
git clone https://github.com/bbakernoaa/GEOSmie.git GEOSmie
git -C GEOSmie checkout c535388356c7e8982ead98deeebff33985254f9c
python tools/generate_modal_optics.py --source GEOSmie \
  --snapshot tools/modal_optics/modal_mie_v1.nc \
  --header src/src_impl/data/generated/ModalOpticsData.hpp
```

Deterministic bake check:

```sh
python tools/generate_modal_optics.py \
  --snapshot tools/modal_optics/modal_mie_v1.nc \
  --header src/src_impl/data/generated/ModalOpticsData.hpp --check
python -m pytest tests/test_modal_optics_generator.py
```

The pinned Mie solver has small roundoff excess of scattering over extinction
in effectively nonabsorbing small-particle cases. Only excess up to `1e-6`
relative is clipped to extinction; larger violations fail generation. No
placeholder optical coefficients are used.

These are visible/NIR starter tables, not full RRTMG bands. The dust calculation
assumes spheres. Nitrate is normalized to complete dry ammonium-nitrate-equivalent
mass, rather than nitrate-ion mass. Full scope, physical assumptions, domains,
and coupling limitations are in [the backend guide](../../docs/guides/generic-modal.md).

GEOSmie is Apache-2.0, with its pymiecoated component carrying the MIT license.
Refractive indices are distributed with the pinned GEOSmie source. See
[`GEOSmie/LICENSE`](https://github.com/bbakernoaa/GEOSmie/blob/c535388356c7e8982ead98deeebff33985254f9c/LICENSE)
and the notices in its pymiecoated source. The snapshot contains derived
numerical results; it does not vendor the solver source.
