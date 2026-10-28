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
in effectively nonabsorbing small-particle cases. The generic snapshot clips
excess up to `1e-6` relative. The current generator and MAM4xx snapshot allow up
to `2e-6` relative (observed maximum below `1.6e-6`) to be clipped to extinction.
Larger violations fail generation. No placeholder optical coefficients are used.

These are visible/NIR starter tables, not full RRTMG bands. The dust calculation
assumes spheres. Nitrate is normalized to complete dry ammonium-nitrate-equivalent
mass, rather than nitrate-ion mass. Full scope, physical assumptions, domains,
and coupling limitations are in [the backend guide](../../docs/guides/generic-modal.md).

The manifests now pin native `RRTMG_30` band edges and GEOSmie's averaging
recipe: 100 equally spaced wavenumber samples per band, arithmetic means for
extinction, scattering, and backscatter, and scattering-weighted asymmetry.
Generation checks every sample against both the material and water refractive
index domains; an uncovered material/band pair is omitted from that material's
band axis and is reported unavailable. Unsupported ranges are never endpoint
clamped or filled at band centers. The current checked-in snapshots have not
yet been regenerated with those band fields, so the native modal interface
reports their bands unavailable until the pinned source and generator
environment are used to regenerate them.

GEOSmie is Apache-2.0, with its pymiecoated component carrying the MIT license.
Refractive indices are distributed with the pinned GEOSmie source. See
[`GEOSmie/LICENSE`](https://github.com/bbakernoaa/GEOSmie/blob/c535388356c7e8982ead98deeebff33985254f9c/LICENSE)
and the notices in its pymiecoated source. The snapshot contains derived
numerical results; it does not vendor the solver source.


## MAM4xx material variant

`mam4xx_manifest.json` defines seven native MAM4xx materials with densities and
kappa matching the pinned upstream revision. SOA, POM and MOM explicitly use
configurable OPAC `waso00` surrogates. Set each material's `refractive_index`
(and `format: wsv` for a three-column input) to change the surrogate, then
regenerate. The output records the exact recipe and input hashes. This does
not change native MAM4xx microphysics or its water-uptake hysteresis.

```sh
python tools/generate_modal_optics.py --source GEOSmie \
  --manifest tools/modal_optics/mam4xx_manifest.json \
  --snapshot tools/modal_optics/mam4xx_mie_v1.nc \
  --header src/src_impl/data/generated/Mam4xxOpticsData.hpp \
  --namespace mie_mam4xx
python tools/generate_modal_optics.py \
  --snapshot tools/modal_optics/mam4xx_mie_v1.nc \
  --header src/src_impl/data/generated/Mam4xxOpticsData.hpp \
  --namespace mie_mam4xx --check
```

The radius range extends to 200 micrometres to cover the coarse-mode mass tail.
The same physical bounds and dry-mass normalization checks apply to both sets.
Normal builds use committed snapshots and deterministic headers; GEOSmie is
required only for offline regeneration. See the
[MAM4xx guide](../../docs/guides/aerosol-backends.md#mam4xx-component-optics)
for configuration, phase selection, validity domain and optical approximations.
