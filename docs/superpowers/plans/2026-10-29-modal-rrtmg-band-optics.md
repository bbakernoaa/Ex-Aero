# Modal RRTMG Band Optics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Ex-Aero's generic modal and MAM4xx component optics usable on native RRTMG shortwave and longwave bands wherever the pinned refractive-index sources cover the full band.

**Architecture:** Keep `computeOptics` on physical wavelengths and use `IOpticsBandCoupling` for the separate RRTMG-30 axis. Generate auditable band-averaged coefficients from the pinned GEOSmie source, retain per-material availability, and never substitute a wavelength or extrapolate unsupported bands.

**Tech Stack:** Python, NumPy, SciPy, xarray, Dask, netCDF4, pinned GEOSmie/pymiecoated, C++23, Catch2 v2.

**Spec:** `docs/superpowers/specs/2026-10-28-cam-sima-aerosol-parity-design.md`

## Global Constraints

- Reference recipe: GEOSmie commit `c535388356c7e8982ead98deeebff33985254f9c` and CAM-SIMA commit `52723a73a2a49e37fdd981f795d08084bb5bfda8` for the target optics scope.
- RRTMG set: exactly 30 native identifiers, bands 1–14 shortwave and 15–30 longwave; API selectors remain zero-based and source IDs remain one-based.
- Retain 100 equally spaced wavenumber samples per band, arithmetic means for extinction/scattering/backscatter, and scattering-weighted asymmetry.
- A material/band is supported only if all samples are within both material and water refractive-index domains; no interpolation beyond source endpoints.
- Preserve deterministic NetCDF-to-header serialization, dry-mass normalization, and atomic failure for positive-mass unsupported components.
- Do not change the existing physical-wavelength `computeOptics` behavior or `IAerosolPackage` vtable.

## Review Focus

- **Partial source coverage:** unsupported bands must be represented as unavailable, never filled with NaNs in baked runtime arrays or clamped values. Test in Tasks 1–3.
- **Coefficient shape:** every generated table is ordered `(radius,rh,band)` with increasing source IDs. Test in Task 2.
- **Asymmetry aggregation:** band asymmetry must use scattering-weighted `g`, not an arithmetic average. Test in Task 1.
- **Band/wavelength separation:** physical wavelength APIs must reject band IDs as metres and vice versa. Test in Task 4.
- **Deterministic optimized evaluation:** batching/parallelization must reproduce scalar pinned GEOSmie results within declared floating tolerance. Test in Task 2.

---

### Task 1: Lock the native-band recipe and reference cases

**Files:**
- Modify: `tests/test_modal_optics_generator.py`
- Reference: `tools/modal_optics/manifest.json`, `tools/modal_optics/mam4xx_manifest.json`, pinned GEOSmie sources

**Interfaces:**
- Consumes: `band_optical_coefficients(radius, growth, density, material, water, mie_function, samples_per_band=100)`.
- Produces: tests pinning RRTMG IDs/edges, source coverage, normalization, and aggregation before optimizing generation.

- [ ] **Step 1: Add `test_band_asymmetry_uses_scattering_weighting`** with a synthetic Mie function whose `qsca` and `asy` vary over the 100 wavelength samples; assert output `band_g` equals `sum(bsca_j*g_j)/sum(bsca_j)` and differs from the unweighted mean. Extend `test_native_rRTMG_band_average_and_extreme_lw_support` to assert 30 ordered bands, 14/16 SW/LW partition, 100 uniform-wavenumber samples, strict material-and-water coverage, NaN intermediate values for unsupported bands, and arithmetic mean `bext/bsca/bbck`.
- [ ] **Step 2: Run** `python -m pytest tests/test_modal_optics_generator.py -q`; expected result before fixes is a failure for any unpinned or incorrectly aggregated case.
- [ ] **Step 3: Correct the manifests/generator only where tests show a mismatch**; maintain the exact pinned band recipe and separate legacy wavelength interpolation behavior.
- [ ] **Step 4: Rerun generator tests**; expected result is all recipe and coverage tests pass on synthetic scalar Mie functions.
- [ ] **Step 5: Commit** as `test: pin modal RRTMG band averaging recipe`.

### Task 2: Make pinned generation reproducible at practical cost

**Files:**
- Modify: `tools/generate_modal_optics.py`
- Create: `tools/benchmark_modal_optics_generation.py`
- Test: `tests/test_modal_optics_generator.py`
- Reference: `/tmp` or another clean checkout at the pinned GEOSmie revision; no solver source is vendored into Ex-Aero.

**Interfaces:**
- Consumes: locked scalar GEOSmie equations and Task 1 expectations.
- Produces: a deterministic batch/chunk evaluation path for `(radius,rh,band_sample)` grids that agrees with the scalar Mie implementation.

- [ ] **Step 1: Create `tools/benchmark_modal_optics_generation.py`** with CLI arguments `--source`, `--manifest`, and `--material`; it must evaluate the selected material on a deterministic 4-radius × 2-RH × 3-band subset using the pinned 100 samples per band and report elapsed time, peak resident memory, GEOSmie revision, Python/package versions, and scalar evaluation count.
- [ ] **Step 2: Add `test_batched_mie_matches_scalar_geosmie`** for small-radius, coarse-radius, low-RH, and high-RH points; assert `qext`, `qsca`, `qb`, and `g` agree with pinned scalar `Mie` outputs using `rtol=1e-12, atol=1e-30`.
- [ ] **Step 3: Profile the pinned implementation** and select one deterministic optimization that preserves equations (chunked evaluation, reusable intermediates, or deterministic independent work scheduling); do not reduce sample count or radius/RH resolution.
- [ ] **Step 4: Implement the optimized evaluator** behind the existing `band_optical_coefficients` contract and preserve band order when chunks complete.
- [ ] **Step 5: Run** `python -m pytest tests/test_modal_optics_generator.py -q` and `python tools/benchmark_modal_optics_generation.py --source GEOSmie --manifest tools/modal_optics/manifest.json --material MODAL_SU`; expected result is scalar-equivalent output, stable memory use, and a recorded runtime before full regeneration.
- [ ] **Step 6: Commit** as `perf: make modal band generation deterministic and bounded`.

### Task 3: Regenerate generic modal and MAM4xx band datasets

**Files:**
- Modify: `tools/modal_optics/modal_mie_v1.nc`
- Modify: `tools/modal_optics/mam4xx_mie_v1.nc`
- Modify: `src/src_impl/data/generated/ModalOpticsData.hpp`
- Modify: `src/src_impl/data/generated/Mam4xxOpticsData.hpp`
- Modify: `tests/test_modal_optics_generator.py`

**Interfaces:**
- Consumes: Task 2 deterministic evaluator and the existing `band_average` manifest recipes.
- Produces: generated per-material `band_source_ids`, `band_bext`, `band_bsca`, `band_bbck`, `band_g`, and exact availability metadata; old wavelength tables remain intact.

- [ ] **Step 1: Run full generic generation** with `python tools/generate_modal_optics.py --source GEOSmie --manifest tools/modal_optics/manifest.json --snapshot tools/modal_optics/modal_mie_v1.nc --header src/src_impl/data/generated/ModalOpticsData.hpp` using the exact pinned source commit.
- [ ] **Step 2: Run full MAM4xx generation** with `python tools/generate_modal_optics.py --source GEOSmie --manifest tools/modal_optics/mam4xx_manifest.json --snapshot tools/modal_optics/mam4xx_mie_v1.nc --header src/src_impl/data/generated/Mam4xxOpticsData.hpp --namespace mie_mam4xx`.
- [ ] **Step 3: Assert generated coverage** per material/band from both manifest domains; verify every unavailable coefficient is excluded from the baked material axis and no band source exceeds the original input domain.
- [ ] **Step 4: Verify deterministic serialization** with both generator `--check` commands and `python -m pytest tests/test_modal_optics_generator.py -q`; expected result is exact header bytes and passing Rayleigh/mass-normalization tests.
- [ ] **Step 5: Record artifact checksums, data-generation runtime, peak memory, and source commit** in `tools/modal_optics/README.md`.
- [ ] **Step 6: Commit** the snapshots, generated headers, tests, and provenance as `data: add native modal RRTMG optics`.

### Task 4: Verify runtime band coverage and parity outputs

**Files:**
- Modify if needed: `src/src_impl/modal/ModalPackage.cpp`, `src/src_impl/loader/MieTableStore.cpp`, `src/src_impl/mam4xx/Mam4xxPackage.cpp`
- Test: `tests/test_modal_optics_bands.cpp`, `tests/test_mam4xx.cpp`
- Modify: `docs/guides/generic-modal.md`, `docs/guides/aerosol-backends.md`, `docs/guides/cam-sima-parity.md`

**Interfaces:**
- Consumes: generated table data and existing `IOpticsBandCoupling` interface.
- Produces: active `OpticsBands` capability for data-backed modal packages, per-component coverage, and SW/LW computation through native IDs.

- [ ] **Step 1: Add failing runtime tests** for full band metadata, representative supported SW/LW material values, unsupported source-domain statuses, zero-mass unsupported components, positive-mass atomic failure, provider-backed internally mixed modes, and wavelength/band-axis separation.
- [ ] **Step 2: Run** `ctest --test-dir build -R 'ModalOpticsBands|ExaeroMam4xx' --output-on-failure`; expected result before any fix is a targeted availability or coefficient mismatch.
- [ ] **Step 3: Correct runtime loading/integration only where data-backed test results require it**; never substitute wavelength tables for native bands.
- [ ] **Step 4: Run the focused tests and both deterministic bake checks**; expected result is enabled capability only when usable band data/provider exists, complete output on supported requests, and unchanged outputs on unsupported positive-mass requests.
- [ ] **Step 5: Update the parity matrix and guides** with actual supported material/band coverage and remaining documented source limitations.
- [ ] **Step 6: Commit** as `feat: enable data-backed modal RRTMG optics`.

### Task 5: Full optics validation

**Files:**
- Test: `tests/test_modal_optics_bands.cpp`, `tests/test_modal_optics_generator.py`, `tests/test_mam4xx.cpp`, existing `tests/test_gocart_optics_bands.cpp`

- [ ] **Step 1: Run** `python -m pytest tests/test_modal_optics_generator.py -q`; expected result is all generator cases pass.
- [ ] **Step 2: Run** `ctest --test-dir build -R '(GocartOpticsBands|ModalOpticsBands|ExaeroMam4xx|ExaeroModalMieDeterminism|ExaeroMam4xxMieDeterminism)' --output-on-failure`; expected result is all available band, wrapper, and data checks pass.
- [ ] **Step 3: Run a full supported CMake build with and without MAM4xx** and record any environment limitation separately from code failures.
- [ ] **Step 4: Commit** final verification/provenance updates as `docs: record modal RRTMG parity coverage`.
