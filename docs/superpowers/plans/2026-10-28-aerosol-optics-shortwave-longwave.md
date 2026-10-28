# Shortwave and Longwave Aerosol Optics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an explicit RRTMG-30 band optics API alongside physical-wavelength optics for bulk and generic modal aerosol, covering 14 shortwave and 16 longwave bands where validated data exist.

**Architecture:** Keep `computeOptics` and the `IModalOpticsProvider` wavelength contract intact. Add a separate public radiation-band metadata and coupling interface implemented by GOCART and generic modal packages, backed by native GEOSmie band averages for bulk and reproducible RRTMG-band data for modal optics. Availability is queried per package component; a band calculation fails atomically if any positive-mass component is unsupported.

**Tech Stack:** C++23, Kokkos, mdspan public views, C and Fortran ISO_C_BINDING, Python/NumPy/xarray NetCDF optics generator, CMake/CTest.

**Spec:** `docs/superpowers/specs/2026-10-28-aerosol-optics-shortwave-longwave-design.md`

## Global Constraints

- “`computeOptics` remains the physical-wavelength operation.”
- “The first supported identifier is `RRTMG_30`; its entries are zero-based at the API boundary and map explicitly to source IDs 1–30.”
- “Band calculations aggregate species, so they fail atomically if any species with positive mass in the requested cell/level lacks the requested band data. On failure the caller's output remains untouched.”
- “No implicit conversion between wavelength and band index is allowed.”
- “The legacy RH-only lookup and ADT path remain available to the physical-wavelength operation; they are not an automatic substitute for missing RRTMG data.”
- “Do not extrapolate beyond the material source domain.”
- “A band may be marked unsupported if its required source/weighting range is missing; it must not be populated from its center wavelength.”
- Preserve `IAerosolPackage` vtable and current optics output slots/layout.

## Review Focus

- A positive-mass species with no RRTMG band source must fail the whole computation without changing any output values. Test in Task 2 (bulk) and Task 4 (modal).
- A zero-mass unsupported component must not block a result. Test in Task 2 and Task 4.
- A physical wavelength must never be interpreted as a band ID, and a band ID must never enter the wavelength interpolation/fallback path. Test in Tasks 2 and 4.
- A request at an unsupported extreme longwave band must report unavailable rather than extrapolate spectral optical constants. Test in Task 3 and Task 4.
- Invalid band-set codes, duplicate/out-of-range requested band indices, buffer shapes, and band metadata counts must fail before any caller output is written. Test in Tasks 1, 2, 4, and 6.

---

## File Map

- Create `src/exaero/RadiationBands.hpp` for Kokkos-free band-set identifiers, region/status enums, and immutable RRTMG-30 metadata.
- Create `src/exaero/OpticsBandCoupling.hpp` for the additive `IOpticsBandCoupling` interface, component availability records, and checked accessor.
- Modify `src/exaero/StateField.hpp` to add the `OpticsBands = 32` capability bit; do not add a virtual or alter the `IAerosolPackage` vtable.
- Modify `src/src_impl/gocart/GocartPackage.hpp/.cpp`, `GocartSpeciesParams.hpp`, and `GocartSolver.cpp` to track spectral-axis semantics, report per-species support, and evaluate native RRTMG bulk fields.
- Modify `tools/modal_optics/manifest.json`, `tools/modal_optics/mam4xx_manifest.json`, `tools/generate_modal_optics.py`, and `tests/test_modal_optics_generator.py` to validate the physical spectral source domain and generate auditable generic-modal and MAM4xx RRTMG band data without midpoint substitution.
- Modify `src/src_impl/modal/ModalSolver.hpp/.cpp` and `src/src_impl/modal/ModalPackage.hpp/.cpp` to upload, evaluate, and atomically commit generic-modal band products while retaining current wavelength operations.
- Modify `src/src_impl/mam4xx/Mam4xxPackage.hpp/.cpp` to delegate the band interface to the configured modal optics package.
- Modify `src/exaero/IAerosolPackage_C.h/.cpp`, `src/exaero/exaero_interface.f90`, and `tests/test_cloud_ice_fortran.f90` or a dedicated optics test to expose band metadata, availability, and computations across C/Fortran.
- Add focused C++ tests under `tests/` and register targets in `tests/CMakeLists.txt`; add public-header checks in `tests/standalone_public_check.cpp`.
- Update `docs/api-contracts.md`, `docs/guides/generic-modal.md`, `docs/guides/aerosol-backends.md`, and `docs/guides/using-exaero-in-ccpp.md` to distinguish wavelength and band coordinates and explain SW/LW usage.

## Tasks

### Task 1: Public RRTMG band metadata and coupling contract

**Files:**
- Create: `src/exaero/RadiationBands.hpp`
- Create: `src/exaero/OpticsBandCoupling.hpp`
- Modify: `src/exaero/StateField.hpp`
- Test: `tests/test_radiation_bands.cpp`
- Test: `tests/standalone_public_check.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Produces `enum class RadiationBandSet : int { RRTMG_30 = 0 }`.
- Produces `enum class RadiationBandRegion : int { Shortwave = 0, Longwave = 1 }`.
- Produces `struct RadiationBandInfo { int source_id; RadiationBandRegion region; double lower_wavenumber_cm1; double upper_wavenumber_cm1; }` for the canonical 30 entries, ordered 14 SW then 16 LW.
- Produces `enum class RadiationBandStatus : int { Available = 0, NotInSource = 1 }`.
- Produces `struct RadiationBandComponentAvailability { std::size_t component_index; std::string component_name; std::vector<RadiationBandStatus> status; }`.
- Produces `std::vector<RadiationBandInfo> selectRadiationBands(RadiationBandSet, std::span<const int> band_indices)`; empty selects all 30, otherwise rejects negative, out-of-range, or duplicate indices and preserves caller order.
- Produces `class IOpticsBandCoupling` with exact methods `std::vector<RadiationBandInfo> radiationBandInfo(RadiationBandSet) const`, `std::vector<RadiationBandComponentAvailability> bandAvailability(RadiationBandSet) const`, and `void computeBandOptics(const EnvironmentalStateView&, const View3D<const double>&, RadiationBandSet, std::span<const int> band_indices, View4D<double>&)`. Empty `band_indices` requests all 30; otherwise values are unique, zero-based indices into the canonical metadata and output band order matches the request.
- Produces checked helper `opticsBandCoupling(IAerosolPackage&)` using `dynamic_cast`, matching the existing modal/removal extension-interface pattern.
- Produces `PackageCapability::OpticsBands = 32`; existing bit values stay unchanged.

- [ ] **Step 1: Write metadata and public-header tests.** Add `test_radiation_band_metadata` asserting RRTMG_30 has 30 entries; indices 0–13 are Shortwave and 14–29 Longwave; source IDs are exactly 1–30; each wavenumber interval has finite positive ordered bounds; unknown enum values throw. Add `test_radiation_band_selection` asserting empty selects all, requested `[14, 2]` preserves that order, and negative, out-of-range, and duplicate indices throw. Include `OpticsBandCoupling.hpp` in standalone public-header compilation.

- [ ] **Step 2: Run the new test to verify it fails.**

Run: `cmake --build build --target exaero_radiation_bands_test && ./build/tests/exaero_radiation_bands_test`

Expected: compilation/test fails because the metadata/interface has not yet been implemented.

- [ ] **Step 3: Implement band metadata and interface.** Populate the 30 source bands and exact SW/LW partition from the pinned RRTMG band definition cited in `tools/geosmie_snapshot/PROVENANCE.md`. Keep public declarations Kokkos-free. Add `OpticsBands = 32` capability without changing existing capability values or the base-class vtable.

- [ ] **Step 4: Run metadata and standalone public-header tests.**

Run: `cmake --build build --target exaero_radiation_bands_test exaero_standalone_public_check && ctest --test-dir build -R 'ExaeroRadiationBands' --output-on-failure`

Expected: the radiation-band CTest passes and the standalone public-header target compiles; all band IDs, regions, and physical bounds match the pinned metadata.

- [ ] **Step 5: Commit.**

```bash
git add src/exaero/RadiationBands.hpp src/exaero/OpticsBandCoupling.hpp src/exaero/StateField.hpp tests/test_radiation_bands.cpp tests/standalone_public_check.cpp tests/CMakeLists.txt
git commit -m "feat: define radiation band optics interface"
```

### Task 2: Bulk native RRTMG band optics and atomic availability

**Files:**
- Modify: `src/src_impl/gocart/GocartSpeciesParams.hpp`
- Modify: `src/src_impl/gocart/GocartPackage.hpp`
- Modify: `src/src_impl/gocart/GocartPackage.cpp`
- Modify: `src/src_impl/gocart/GocartSolver.cpp`
- Create: `tests/test_gocart_optics_bands.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `IOpticsBandCoupling`, `RadiationBandSet`, and canonical RRTMG metadata from Task 1.
- Produces: `GocartPackage : public IOpticsBandCoupling`.
- `radiationBandInfo(RRTMG_30)` returns the shared canonical metadata.
- `bandAvailability(RRTMG_30)` returns 30 statuses for each configured GOCART species, preserving species index/name order.
- `computeBandOptics(env, state, RRTMG_30, band_indices, out)` returns `(cell, level, n_requested, NUM_OPTICS)`; empty `band_indices` requests all 30. It uses temporary output/status storage and commits only if every positive-mass species has a source value for every requested cell/level/band.

- [ ] **Step 1: Write focused failing tests.** Add `test_gocart_rrtmg_band_exact_values` comparing exact radius/RH/native band outputs with pinned per-species RRTMG goldens; `test_gocart_band_availability_species` asserting baked DU/SS/SU/BC/OC/NI support and BR unsupported; `test_gocart_band_missing_mass_is_atomic` setting positive BR mass and asserting an error and unchanged output; `test_gocart_band_zero_mass_missing_source` asserting zero BR mass does not block; and `test_gocart_band_coordinate_spaces_are_separate` asserting API index 2 maps to source ID 3 in `computeBandOptics`, while a 550-nm request through `computeOptics` follows the existing physical-wavelength ADT path.

- [ ] **Step 2: Run the focused tests and observe failure.**

Run: `cmake --build build --target exaero_gocart_optics_bands_test && ./build/tests/exaero_gocart_optics_bands_test`

Expected: compile/test fails because the band interface and implementation do not exist.

- [ ] **Step 3: Track spectral axis kind and implement the direct band kernel.** Store whether configured/baked spectral axes are RRTMG IDs or physical wavelengths. For RRTMG_30, pass the native source ID to table interpolation, use table `bext`, `ssa`, `g`, and `bbck` fields, derive scattering and absorption consistently, calculate output AOTs over layer thickness, and never invoke ADT or the legacy RH lookup for missing band values. Leave `computeOptics` wavelength behavior unchanged.

- [ ] **Step 4: Make unsupported positive-mass components fail atomically.** Run into a temporary result plus unsupported flags. Check all cells/levels/bands and throw before copying to caller memory if any unsupported component has positive mass. Ignore unsupported zero-mass components. Ensure `BACKSCATTER_COEFF` and `LIDAR_BACKSCATTER` use source backscatter for the native band path rather than the spectral path's heuristic.

- [ ] **Step 5: Run the focused tests and existing bulk optics tests.**

Run: `cmake --build build --target exaero_gocart_optics_bands_test exaero_test_runner && ctest --test-dir build -R 'GocartOpticsBands|ExaeroUnitTests' --output-on-failure`

Expected: all new band tests and the existing GOCART optics/MIE tests pass; atomic failure preserves the sentinel output bytes.

- [ ] **Step 6: Commit.**

```bash
git add src/src_impl/gocart/GocartSpeciesParams.hpp src/src_impl/gocart/GocartPackage.hpp src/src_impl/gocart/GocartPackage.cpp src/src_impl/gocart/GocartSolver.cpp tests/test_gocart_optics_bands.cpp tests/CMakeLists.txt
git commit -m "feat: compute bulk aerosol optics on RRTMG bands"
```

### Task 3: Reproducible modal SW/LW spectral and band data

**Files:**
- Modify: `tools/modal_optics/manifest.json`
- Modify: `tools/generate_modal_optics.py`
- Modify: `tests/test_modal_optics_generator.py`
- Modify: `tools/modal_optics/README.md`
- Update: `tools/modal_optics/modal_mie_v1.nc`
- Update: `src/src_impl/data/generated/ModalOpticsData.hpp`
- Modify: `tools/modal_optics/mam4xx_manifest.json`
- Update: `tools/modal_optics/mam4xx_mie_v1.nc`
- Update: `src/src_impl/data/generated/Mam4xxOpticsData.hpp`

**Interfaces:**
- Consumes: RRTMG band metadata from Task 1.
- Produces: deterministic generic-modal and MAM4xx material records with physical wavelengths spanning validated SW/LW source coverage and per-material RRTMG band arrays; metadata identifies weighting recipe, spectral domain, input hashes, and unsupported bands.
- Produces `generate_modal_optics.py --check` coverage for the updated checked-in snapshot/header.

- [ ] **Step 1: Add generator tests for source coverage and RRTMG band generation.** Test a synthetic material/source where all spectral nodes cover a band and compare its generated band average against a hand-calculated weighted reference; test a missing extreme-LW source interval marks that band unsupported; test no band is evaluated from its center; test deterministic serialization and provenance fields.

- [ ] **Step 2: Run generator tests and observe failure.**

Run: `python3 -m pytest tests/test_modal_optics_generator.py -q`

Expected: new band-generation tests fail because the manifest/generator only supports point wavelengths.

- [ ] **Step 3: Pin RRTMG band averaging inputs and extend the manifest.** Record the authoritative RRTMG band boundaries/weight source and its revision or content hash. Extend `manifest.json` with a physical spectral grid across validated shortwave and longwave material data, band set `RRTMG_30`, weighting data, and output domain metadata. Reject a material whose source refractive-index range does not cover an advertised band; never extrapolate refractive indices or optical results.

- [ ] **Step 4: Generate band products using the pinned averaging recipe.** Implement deterministic generation of per-band extinction, scattering, absorption, asymmetry, and backscatter from spectral calculations over supported intervals. Handle asymmetry with scattering-weighted averaging. Preserve current dry-mass normalization and apply hygroscopic growth exactly once. Preserve unsupported-band status in the generated artifact.

- [ ] **Step 5: Regenerate and verify snapshots.**

Run: `python3 -m pytest tests/test_modal_optics_generator.py -q && python3 tools/generate_modal_optics.py --source GEOSmie --manifest tools/modal_optics/manifest.json --snapshot tools/modal_optics/modal_mie_v1.nc --header src/src_impl/data/generated/ModalOpticsData.hpp && python3 tools/generate_modal_optics.py --source GEOSmie --manifest tools/modal_optics/mam4xx_manifest.json --snapshot tools/modal_optics/mam4xx_mie_v1.nc --header src/src_impl/data/generated/Mam4xxOpticsData.hpp --namespace mie_mam4xx && python3 tools/generate_modal_optics.py --snapshot tools/modal_optics/modal_mie_v1.nc --header src/src_impl/data/generated/ModalOpticsData.hpp --check && python3 tools/generate_modal_optics.py --snapshot tools/modal_optics/mam4xx_mie_v1.nc --header src/src_impl/data/generated/Mam4xxOpticsData.hpp --namespace mie_mam4xx --check`

Expected: generator tests pass and the deterministic snapshot/header check reports no diff.

- [ ] **Step 6: Commit.**

```bash
git add tools/modal_optics/manifest.json tools/modal_optics/mam4xx_manifest.json tools/generate_modal_optics.py tests/test_modal_optics_generator.py tools/modal_optics/README.md tools/modal_optics/modal_mie_v1.nc tools/modal_optics/mam4xx_mie_v1.nc src/src_impl/data/generated/ModalOpticsData.hpp src/src_impl/data/generated/Mam4xxOpticsData.hpp
git commit -m "feat: generate modal shortwave and longwave optics"
```

### Task 4: Generic modal band evaluation and atomic support

**Files:**
- Modify: `src/src_impl/modal/ModalSolver.hpp`
- Modify: `src/src_impl/modal/ModalSolver.cpp`
- Modify: `src/src_impl/modal/ModalPackage.hpp`
- Modify: `src/src_impl/modal/ModalPackage.cpp`
- Modify: `src/exaero/ModalCoupling.hpp`
- Create: `tests/test_modal_optics_bands.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: shared metadata/API from Task 1 and generated modal band data from Task 3.
- Produces: `ModalPackage : public IOpticsBandCoupling` with component availability keyed by configured mode/material components.
- Produces optional `class IModalBandOpticsProvider` with exact method `ModalOpticalResult evaluateBand(const ModalConfig&, std::size_t mode, const RadiationBandInfo&, double temperature, double rh, std::span<const double> species_mass, std::span<const double> properties) const` and `IModalCoupling::setModalBandOpticsProvider(std::shared_ptr<const IModalBandOpticsProvider>)`; the existing `IModalOpticsProvider` wavelength interface remains unchanged.
- `computeBandOptics(env, state, RRTMG_30, band_indices, out)` emits `(cell, level, n_requested, NUM_OPTICS)` and per-band AOTs with atomic output commit; an empty selection requests all 30.

- [ ] **Step 1: Write focused modal band tests.** Add `test_modal_band_optics_matches_generated_reference` for one-component and internally mixed modes; `test_modal_band_availability_is_per_component`; `test_modal_unsupported_positive_mass_is_atomic`; `test_modal_unsupported_zero_mass_is_ignored`; `test_modal_extreme_longwave_does_not_extrapolate`; and `test_modal_wavelength_provider_remains_unchanged` comparing existing physical-wavelength output before/after band work.

- [ ] **Step 2: Run focused tests and observe failure.**

Run: `cmake --build build --target exaero_modal_optics_bands_test && ./build/tests/exaero_modal_optics_bands_test`

Expected: compile/test fails because modal tables have no band axis or band evaluator.

- [ ] **Step 3: Add generated band values and support metadata to modal tables.** Extend private `modal::Table` and `Solver` upload to carry the RRTMG values/status separately from the physical wavelength axis. Validate exact band count and units during initialization; the two coordinate spaces remain distinct.

- [ ] **Step 4: Implement modal band distribution integration.** Integrate generated component band properties against the configured modal lognormal size distribution using the existing optical size weights. Use `IModalBandOpticsProvider::evaluateBand(const ModalConfig&, std::size_t mode, const RadiationBandInfo&, double temperature, double rh, std::span<const double> species_mass, std::span<const double> properties)` for internal-mixture modes; otherwise report that mode unsupported. Maintain mass-normalized coefficients, scattering-weighted asymmetry, source backscatter, and once-only water growth.

- [ ] **Step 5: Enforce atomic validation and commit.** Use temporary device/host outputs and a per-cell/level/band unsupported flag. Fail if any positive-mass component in any requested cell/level is unsupported; leave the caller output unchanged. An unsupported zero-mass component does not block. Update modal capabilities to advertise band optics only when the configured package has usable data.

- [ ] **Step 6: Run focused and existing modal optics tests.**

Run: `cmake --build build --target exaero_modal_optics_bands_test exaero_modal_test && ctest --test-dir build -R 'ModalOpticsBands|Modal' --output-on-failure`

Expected: new modal band tests and existing modal wavelength/provider/optics tests pass.

- [ ] **Step 7: Commit.**

```bash
git add src/src_impl/modal/ModalSolver.hpp src/src_impl/modal/ModalSolver.cpp src/src_impl/modal/ModalPackage.hpp src/src_impl/modal/ModalPackage.cpp src/exaero/ModalCoupling.hpp tests/test_modal_optics_bands.cpp tests/CMakeLists.txt
git commit -m "feat: compute modal aerosol optics on RRTMG bands"
```

### Task 5: Native MAM4xx band optics delegation

**Files:**
- Modify: `src/src_impl/mam4xx/Mam4xxPackage.hpp`
- Modify: `src/src_impl/mam4xx/Mam4xxPackage.cpp`
- Modify: `tests/test_mam4xx.cpp`

**Interfaces:**
- Consumes: Task 1 band interface and Task 4 generic modal implementation.
- Produces: `Mam4xxPackage : public IOpticsBandCoupling`; metadata and per-mode/material availability forward to the configured modal optics package.
- `computeBandOptics(env,state,set,band_indices,out)` forwards the same interstitial/cloudborne field view used by the existing wavelength optics method.
- `Mam4xxPackage::capabilities()` adds `OpticsBands` only when the configured optics package implements `IOpticsBandCoupling`.

- [ ] **Step 1: Write conditional MAM4xx delegation tests.** Add `test_mam4xx_band_optics_delegates_to_modal_provider` comparing band metadata, availability, and values with the configured component-table modal package; add `test_mam4xx_band_optics_requires_band_capable_tables` asserting unconfigured or unsupported optics throws and preserves output.

- [ ] **Step 2: Run the MAM4xx tests and observe failure.**

Run: `cmake --build build --target exaero_mam4xx_test && ./build/tests/exaero_mam4xx_test`

Expected: compile/test fails because `Mam4xxPackage` does not implement the band interface.

- [ ] **Step 3: Delegate metadata, support, and computation.** Dynamic-cast the configured `optics_` package to `IOpticsBandCoupling`; forward metadata/availability and computation using the same phase-trimmed field view as `computeOptics`. Throw the established fatal error when the selected optics treatment lacks band support. Preserve `IAerosolPackage` vtable and use the existing opt-in optics capability behavior.

- [ ] **Step 4: Run MAM4xx and modal band tests.**

Run: `cmake --build build --target exaero_mam4xx_test exaero_modal_optics_bands_test && ctest --test-dir build -R 'Mam4xx|ModalOpticsBands' --output-on-failure`

Expected: delegation matches direct modal values and invalid/unconfigured paths remain atomic.

- [ ] **Step 5: Commit.**

```bash
git add src/src_impl/mam4xx/Mam4xxPackage.hpp src/src_impl/mam4xx/Mam4xxPackage.cpp tests/test_mam4xx.cpp
git commit -m "feat: expose RRTMG band optics through MAM4xx"
```

### Task 6: C/Fortran integration, documentation, and full verification

**Files:**
- Modify: `src/exaero/IAerosolPackage_C.h`
- Modify: `src/exaero/IAerosolPackage_C.cpp`
- Modify: `src/exaero/exaero_interface.f90`
- Modify: `src/exaero/StateField.hpp`
- Modify: `tests/CMakeLists.txt`
- Create: `tests/test_optics_bands_fortran.f90`
- Modify: `tests/standalone_public_check.cpp`
- Modify: `docs/api-contracts.md`
- Modify: `docs/guides/generic-modal.md`
- Modify: `docs/guides/aerosol-backends.md`
- Modify: `docs/guides/using-exaero-in-ccpp.md`

**Interfaces:**
- Consumes: the completed C++ `IOpticsBandCoupling` interface and package implementations from Tasks 1–4.
- Produces: `exaero_get_radiation_band_info`, `exaero_get_optics_band_availability`, and `exaero_compute_band_optics` with CCPP-style `errmsg`/`errflg`; matching `bind(c)` interfaces and `exaero_RADIATION_BAND_SET_RRTMG_30`, region/status constants, and `exaero_CAP_OPTICS_BANDS = 32` in Fortran/C.
- The compute wrapper accepts the package handle, dimensions, band-set code, environmental/state buffers, and output `(cell, level, 30, NUM_OPTICS)`; it never accepts a wavelength array.
- C metadata signature: `exaero_get_radiation_band_info(pkg, band_set_code, band_capacity, n_bands, n_sw, n_lw, source_ids, regions, lower_wavenumber_cm1, upper_wavenumber_cm1, errmsg, errflg)`.
- C availability signature: `exaero_get_optics_band_availability(pkg, band_set_code, component_capacity, component_name_len, band_capacity, n_components, n_bands, component_indices, component_names, status, errmsg, errflg)`; `status` is column-major `(component_capacity, band_capacity)`.
- C compute signature: `exaero_compute_band_optics(pkg, num_cells, num_levels, num_fields, band_set_code, num_requested_bands, band_indices, temp_ptr, pres_ptr, dens_ptr, rh_ptr, thick_ptr, state_ptr, optics_ptr, errmsg, errflg)`; `num_requested_bands=0` requests all 30, otherwise the unique zero-based selected indices define output axis order.

- [ ] **Step 1: Write C/Fortran wrapper tests.** Assert metadata counts and 14/16 split, source-ID mapping, component support reporting, unsupported-band errors, exact output layout, and atomic output preservation through C and Fortran entry points.

- [ ] **Step 2: Run the wrapper test and observe failure.**

Run: `cmake --build build --target exaero_optics_bands_fortran_test && ctest --test-dir build -R 'OpticsBandsFortran' --output-on-failure`

Expected: compile/test fails because the C and Fortran band entry points do not exist.

- [ ] **Step 3: Implement the C ABI and Fortran bindings.** Validate the band-set code, dimensions, pointers, and buffer extents before dispatch. Use an `IOpticsBandCoupling` dynamic cast so unsupported packages return the established fatal error. Keep output copying atomic when C++ throws.

- [ ] **Step 4: Document both coordinate paths.** Update API contracts and guides with wavelength units versus zero-based RRTMG band indices, band metadata/status queries, SW/LW counts, component support, atomic failure semantics, native versus generated data provenance, spectral-domain limits, and a CCPP/Fortran usage example.

- [ ] **Step 5: Run full verification.**

Run: `python3 -m pytest tests/test_modal_optics_generator.py -q && python3 tools/generate_mie_tables.py --check src/src_impl/data/generated/MieTableData.hpp && cmake --build build && ctest --test-dir build --output-on-failure`

Expected: deterministic generation checks, all C++ tests, all Fortran tests, and CTest pass; the existing physical-wavelength APIs retain their results and the RRTMG_30 APIs pass exact table/generator comparisons.

- [ ] **Step 6: Commit.**

```bash
git add src/exaero/IAerosolPackage_C.h src/exaero/IAerosolPackage_C.cpp src/exaero/exaero_interface.f90 src/exaero/StateField.hpp tests docs/api-contracts.md docs/guides/generic-modal.md docs/guides/aerosol-backends.md docs/guides/using-exaero-in-ccpp.md
git commit -m "feat: expose shortwave and longwave aerosol optics"
```

## Execution Notes

- Task 3's source and weighting audit is a hard correctness gate. If any RRTMG band or material lacks validated optical-constant/weight coverage, record it as unavailable in generated data and test that Task 4 fails atomically for positive mass. Do not fill a gap from a center wavelength or extrapolation.
- The generated modal band path must preserve the existing physical-wavelength provider and default output values at its current five wavelength nodes.
- Re-run snapshot generation only with the pinned, clean GEOSmie revision named in the manifests; generated artifacts must retain exact input hashes and provenance.
- Full CMake/Fortran verification may be environment-limited; report the exact commands and tool availability rather than marking those checks as passed without execution.
