# Modal Diagnostics, Physics, and Host-Removal Adapters Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Match Ex-Aero's generic modal diagnostics and implemented host-facing physical calculations to the corresponding working CAM-SIMA modal behavior, and provide pure tendency adapters for host-computed dry-deposition velocities and wet-scavenging coefficients.

**Architecture:** Keep Ex-Aero's independent mass/number state and host-owned timestep order. Use CAM-SIMA formula outputs as the reference, while implementing behavior through Ex-Aero's modal, cloud/ice, bulk, C, and Fortran APIs. The host computes dry-deposition velocities and species/bin/mode-specific wet-scavenging coefficients and applies accepted tendencies. Ex-Aero preserves bulk species/bin differences and independent modal mass/number moments; components of one uniformly internally mixed mode share its mass coefficient.

**Tech Stack:** C++23, Kokkos, mdspan, Catch2 v2, Fortran C interoperability, Python for reference calculations and report checks.

**Spec:** `docs/superpowers/specs/2026-10-28-cam-sima-aerosol-parity-design.md`

## Global Constraints

- Numerical reference: CAM-SIMA `development` commit `52723a73a2a49e37fdd981f795d08084bb5bfda8`.
- Physics dependency reference: `ESCOMP/atmospheric_physics` commit `6ee9b89afb294408529bea5b72c7bab47c0c55eb`, pinned by CAM-SIMA at `src/physics/ncar_ccpp`. Include `schemes/aerosol/core/wetdep.F90` in the capability crosswalk.
- MKPP/host code retains chemical formation, equilibrium thermodynamics, positivity, accepted/rejected steps, and chemical integration.
- Bulk number remains diagnosed from its configured fixed lognormal distribution; modal number remains prognostic.
- Keep host-owned cloudborne state and explicit transfer acceptance. Diagnostic and removal-tendency helpers leave aerosol state unchanged; `enforceModalBounds` remains an explicit correction API that may change modal number when called by the host.
- The simple dry-deposition input remains one finite nonnegative velocity `[m/s]` per cell, level, and bulk species/bin or modal mode, returning `-velocity/layer_thickness * concentration`. This path applies one population velocity to mass and number and requires finite positive layer thickness.
- Wet-removal inputs are finite nonnegative first-order coefficients `[1/s]` per bulk species/bin or, for modal, separate mass and number coefficients per mode. Ex-Aero maps each coefficient to that species/bin or modal moment without averaging distinct populations, applying updates, or changing output buffers on failure.
- Bulk number tendency follows each species/bin's own mass coefficient and fixed number conversion. Modal number is independent. A uniformly internally mixed mode loses all its constituents together under its common fractional mass coefficient. Its mass and number coefficients can differ because scavenging weights those size-distribution moments differently; both must derive from a consistent host closure. Separate ambient/cloudborne state snapshots and rates prevent implicit double removal.
- The host may use Ex-Aero's size, composition, hygroscopicity, and activation diagnostics to calculate wet-scavenging coefficients from its precipitation and cloud physics. No common rate is imposed across bulk species/bins or modal modes. Arbitrary component-specific fractional scavenging within a uniformly internally mixed mode is outside the modal adapter contract.
- A per-mode coefficient is the result of the host's composition-dependent physics, not a prescribed species-independent constant. CAM-SIMA's `aero_convproc` computes mode hygroscopicity from constituent mass/density/hygroscopicity, with a POM-specific hygroscopicity option affecting primary BC/POM scavenging. Audit Ex-Aero's public composition, effective hygroscopicity, size, density, and activation outputs for the corresponding host inputs. Solubility/scavenging factors must not be inferred from species names or ice-candidate roles.
- CAM-SIMA's property API returns per-bin `scav_diam = dgnum`; its actual wet-deposition core is in the pinned physics dependency. `get_bcscavcoefs(m,...)` returns separate number and volume/mass coefficients per mode; `wetdepa_v2` distinguishes reservoirs and process terms. Retain those representation boundaries, while leaving the collection/precipitation/resuspension scheme external. The local-rate adapter is not numerical parity with the full wet-deposition algorithm.
- CAM-SIMA stores constituent masses and one number field per mode and computes separate mass/number activation fractions (`actfracs`). These support distinct modal moments, but do not establish a wet-deposition rate law. Local sink tendencies do not replace host precipitation transport, evaporation/resuspension, or surface-budget accounting.
- Preserve the existing `IAerosolPackage` vtable and C/Fortran array layout contracts.
- Do not implement CAM-SIMA bulk stubs as if they were working reference behavior.

## Review Focus

- **Empty modes:** zero mass and number must remain empty through bounds, activation, and removal diagnostics. Test in Task 2 and Task 3.
- **Bound endpoints:** exact lower/upper radius limits must not create off-by-one or repeated number corrections. Test in Task 2.
- **Mixed composition:** type weights must use the correct category mass without double-counting overlapping configured roles. Test in Task 5.
- **Rejected host steps:** diagnostic and tendency calls must leave input state unchanged. Test in Tasks 3 and 5.
- **Different removal rates:** preserve distinct bulk species/bin and modal-mode rates, including independent modal mass/number rates. A common mass coefficient within each internally mixed mode preserves component ratios. Test in Task 3.
- **Invalid removal inputs:** negative/nonfinite velocity or scavenging coefficient or concentration, nonpositive/nonfinite thickness on the velocity path, aliased buffers, and nonfinite calculated tendencies must fail atomically. Test in Task 3.
- **Language layouts:** C and Fortran wrappers must preserve zero-based selectors, column-major extents, units, and statuses. Test in Task 6.

---

### Task 1: Pin the CAM-SIMA capability crosswalk and formula fixtures

**Files:**
- Create: `docs/guides/cam-sima-parity.md`
- Reference: `src/aerosol/modal_aerosol_properties_mod.F90` and `src/aerosol/bulk_aerosol_properties_mod.F90` at the pinned CAM-SIMA revision
- Reference: `schemes/aerosol/core/wetdep.F90` at the pinned `ESCOMP/atmospheric_physics` revision
- Reference: `schemes/aerosol/core/aero_convproc.F90` at the pinned physics revision, especially its composition-weighted hygroscopicity and mode/moment activation mapping

**Interfaces:**
- Produces a row for each target behavior: source routine, input/state assumptions, Ex-Aero API, status, and reference case.
- Use statuses `Equivalent`, `Partial`, `Intentional difference`, `Missing`, and `Not applicable`.

- [ ] **Step 1: Record the reference methods and state assumptions** for `apply_number_limits`, `actfracs`, `scav_diam`, `soluble`/`hydrophilic`, `resuspension_resize`, `rebin_bulk_fluxes`, and `icenuc_updates_*`; record only implemented CAM-SIMA bulk behavior and explicitly label its stubs. Crosswalk the dependency's `get_bcscavcoefs`, `init_bcscavcoef`, `wetdepa_v2`, and species-dependent hygroscopicity in `aero_convproc` as implemented host-owned scheme physics rather than absent features or proposed Ex-Aero rate ports. Identify how the host obtains composition and wet-scavenging-related properties through Ex-Aero's current APIs, and record any actual metadata gap separately from the tendency adapter.
- [ ] **Step 2: Add source-pinned reference cases** to the matrix for lower/upper number bounds, modal number/mass activation fractions, per-bin `scav_diam = dgnum`, material redistribution, and deposition rebinning; include concrete inputs, expected outputs, units, and tolerances copied from the pinned formulas. Record the wet-deposition core's separate number/volume moments, reservoirs, and coefficient units. List the host-supplied `/s` wet-rate adapter as an integration addition rather than a direct implementation of CAM-SIMA's precipitation-normalized coefficients and finite-step tendency calculation.
- [ ] **Step 3: Review each Ex-Aero mapping** against `ModalRemoval.hpp`, `CloudIceCoupling.hpp`, `BulkRemoval.hpp`, and their C/Fortran surfaces; record representation differences rather than rewriting APIs to match names.
- [ ] **Step 4: Verify the matrix** against the pinned source URLs and run `git diff --check`; expected result is no unmapped implemented CAM-SIMA modal method in the stated target scope.
- [ ] **Step 5: Commit** `docs/guides/cam-sima-parity.md` as `docs: map CAM-SIMA aerosol parity scope`.

### Task 2: Match modal number bounds and size diagnostics

**Files:**
- Modify if needed: `src/src_impl/modal/ModalRemoval.cpp`
- Modify if needed: `src/src_impl/modal/ModalPackage.cpp`
- Test: `tests/test_modal_removal.cpp`, `tests/test_modal.cpp`
- Reference: matrix rows for `apply_number_limits` and modal size calculations

**Interfaces:**
- Consumes: configured per-mode dry-radius bounds and existing `enforceModalBounds`, `inspectModalBounds`, and `modalScavengingSizes` APIs.
- Produces: documented/tested numerical parity for dry geometric bounds, number correction, and scavenging diameters; no new public signature unless the matrix identifies a true contract gap.

- [ ] **Step 1: Add `test_cam_sima_modal_number_bounds_reference`** in `tests/test_modal_removal.cpp`: configure `sigma_g=1.5`, bounds `[5e-8,2e-7] m`, and mass/number representing `r_g=4e-7 m`; assert enforced number is `8*N_initial`, `delta_number=+7*N_initial`, and every material mass is unchanged. Add below-minimum, exact-boundary, in-range, and empty-mode cases with expected statuses and zero correction for valid/empty states.
- [ ] **Step 2: Add `test_cam_sima_modal_scavenging_sizes_reference`**: configure one material with `kappa=0.5` and `water_uptake=KappaEquilibrium`; assert dry output equals twice the state-derived dry number-median radius and wet output equals dry diameter multiplied by `(1+kappa*RH/(1-RH))^(1/3)` at RH 0, 0.5, and the documented RH cap 0.98. Add `test_modal_composition_changes_host_scavenging_properties` in `tests/test_modal.cpp`: keep mode total dry volume, number, and width fixed while changing volume fractions of materials with kappa `{0.6,0}` from `{0.75,0.25}` to `{0.25,0.75}`. Assert effective kappa changes from `0.45` to `0.15`, dry size stays fixed, and equilibrium wet size differs at nonzero RH. This verifies species-dependent information available to the host without calculating a wet-removal coefficient inside Ex-Aero.
- [ ] **Step 3: Build and run** `cmake --build build --target exaero_modal_removal_test exaero_modal_test && ctest --test-dir build -R 'ExaeroModalRemoval|ExaeroModal' --output-on-failure`; a passing result means the current implementation already matches these cases, while a failure must identify the exact mismatching value before a code change.
- [ ] **Step 4: Correct the smallest formula or unit mismatch** in `ModalRemoval.cpp` or `ModalPackage.cpp`; retain Ex-Aero's radius-based configuration and document its conversion to CAM-SIMA's volume-to-number factors in `docs/guides/cam-sima-parity.md`.
- [ ] **Step 5: Rebuild and rerun the focused tests**; expected result is all reference cases pass and diagnostics leave state unchanged.
- [ ] **Step 6: Commit** source/tests/matrix changes as `fix: match CAM-SIMA modal size diagnostics`.

### Task 3: Return host-computed dry and wet removal tendencies

**Files:**
- Modify: `src/exaero/BulkRemoval.hpp`, `src/src_impl/gocart/BulkRemoval.cpp`
- Modify: `src/exaero/ModalRemoval.hpp`, `src/src_impl/modal/ModalRemoval.cpp`
- Modify: `src/exaero/BulkRemoval_C.h`, `src/src_impl/gocart/BulkRemoval_C.cpp`
- Modify: `src/exaero/ModalRemoval_C.h`, `src/src_impl/modal/ModalRemoval_C.cpp`
- Modify: `src/exaero/exaero_interface.f90`
- Test: `tests/test_bulk_removal.cpp`, `tests/test_modal_removal.cpp`, `tests/test_backend_fortran.f90`, `tests/test_modal_fortran.f90`

**Interfaces:**
- Add C++ entry points `bulkRemovalVelocityTendencies(IAerosolPackage&, const View2D<const double>& layer_thickness, const View3D<const double>& state, const View3D<const double>& velocity, View3D<double>& mass_tendency, View3D<double>& number_tendency)` and `modalRemovalVelocityTendencies(IAerosolPackage&, const View2D<const double>& layer_thickness, const View3D<const double>& state, const View3D<const double>& velocity, View3D<double>& mass_tendency, View3D<double>& number_tendency)`. These retain the previously selected simple velocity contract.
- Add `bulkRemovalRateTendencies(IAerosolPackage&, const View3D<const double>& state, const View3D<const double>& mass_scavenging_rate, View3D<double>& mass_tendency, View3D<double>& number_tendency)` and `modalRemovalRateTendencies(IAerosolPackage&, const View3D<const double>& state, const View3D<const double>& mass_scavenging_rate, const View3D<const double>& number_scavenging_rate, View3D<double>& mass_tendency, View3D<double>& number_tendency)`. Bulk rates are `(cell,level,species)` and independently specified per species/bin. Both modal rate inputs are `(cell,level,mode)`: one common mass coefficient per mode and a separate number coefficient. There is no component-specific modal wet-rate input.
- Add these C entry points (with `nc,nl,nf` bulk extents and `nc,nl,nf,nm` modal extents):
  ```cpp
  void exaero_bulk_removal_velocity_tendencies(
      exaero_package_t pkg, int nc, int nl, int nf, const double *state,
      const double *velocity, const double *layer_thickness,
      double *mass_tendency, double *number_tendency, char *errmsg,
      int *errflg);
  void exaero_modal_removal_velocity_tendencies(
      exaero_package_t pkg, int nc, int nl, int nf, int nm,
      const double *state, const double *velocity,
      const double *layer_thickness, double *mass_tendency,
      double *number_tendency, char *errmsg, int *errflg);
  void exaero_bulk_removal_rate_tendencies(
      exaero_package_t pkg, int nc, int nl, int nf, const double *state,
      const double *mass_scavenging_rate, double *mass_tendency,
      double *number_tendency, char *errmsg, int *errflg);
  void exaero_modal_removal_rate_tendencies(
      exaero_package_t pkg, int nc, int nl, int nf, int nm,
      const double *state, const double *mass_scavenging_rate,
      const double *number_scavenging_rate,
      double *mass_tendency, double *number_tendency, char *errmsg,
      int *errflg);
  ```
  Flat buffers use column-major extents: thickness is `(cell,level)`; bulk arrays are `(cell,level,species)`; modal state and mass output are `(cell,level,nfield)`; modal velocity, both rate inputs, and number output are `(cell,level,mode)`. Number slots in modal mass output are zero. Add matching `bind(c)` interfaces with read-only inputs and inout outputs.
- The velocity adapters return `-velocity/layer_thickness * concentration`. Bulk rate tendencies are `dM_s/dt = -lambda_s*M_s`. Modal tendencies are `dM_ms/dt = -lambda_mass_m*M_ms` for each constituent s of mode m, and `dN_m/dt = -lambda_number_m*N_m`. All coefficients are `[1/s]`. Bulk number tendency is diagnosed from that species/bin's mass tendency under its configured fixed lognormal distribution. Validate dimensions, finite nonnegative inputs, finite positive thickness for the velocity path, and finite calculated tendencies before publishing outputs. Require disjoint output buffers; preserve outputs and state on failure. Zeroing one coefficient disables only its species/bin or modal moment sink. The host applies tendencies only after accepting its timestep.
- These effective `/s` coefficient inputs are an Ex-Aero host-integration proposal. CAM-SIMA's dependency supports the distinct per-mode number/volume structure in `get_bcscavcoefs`; it is not an interchangeable rate API: `wetdepa_v2` accepts precipitation-normalized `/mm` scavenging coefficients and calculates finite-step tendencies with additional inputs. Common constituent-mass removal follows the mode's uniform internal-mixing assumption. Distinct modes, bulk species/bins, and ambient/cloudborne reservoirs retain their own host-computed coefficients. Ex-Aero computes no precipitation, collection law, or size-dependent scavenging integral, and makes no claim that two arbitrary moment coefficients correspond to a valid collection kernel. The host owns that consistency and explicit accepted-step size bounds, precipitation transport, evaporation/resuspension, and surface budgets.

- [ ] **Step 1: Add bulk velocity and species-rate tests**: retain bulk mass `2e-9 kg/m3`, velocity `0.01 m/s`, and thickness `10 m`, giving mass tendency `-2e-12 kg/m3/s` and number tendency `-1e-3*N`. Add two species/bins with masses `{2e-9,3e-9}` and coefficients `{5e-4,2e-4} /s`; assert mass tendencies `{-1e-12,-6e-13}` and species-specific diagnosed number tendencies. Zeroing one coefficient must leave the other species' sink intact. All inputs stay unchanged.
- [ ] **Step 2: Add modal velocity and independent-moment tests**: retain the agreed velocity case. Use material masses `{2e-9,3e-9}`, one mode-mass coefficient `5e-4 /s`, number `5e7`, and independent number coefficient `1e-4 /s`; assert mass tendencies `{-1e-12,-1.5e-12}` and number tendency `-5e3`, with zero number slots in mass output. Assert component mass ratios are preserved. Add two modes with different mass/number coefficient pairs to verify mode/field indexing and ensure rates are not shared across modes. Use separate ambient/cloudborne snapshots and coefficients to check that neither call affects the other reservoir.
- [ ] **Step 3: Add zero and invalid-input tests**: a zero coefficient, or zero velocity with valid thickness, gives zero only for its species/bin or modal moment; other supplied sinks are preserved. Negative/NaN velocity or coefficient, incorrect per-mode rate extents, zero/nonfinite thickness on velocity path, negative/nonfinite concentration, aliased outputs, and finite-input overflow that produces a nonfinite tendency fail and preserve every output sentinel.
- [ ] **Step 4: Add C and Fortran round-trip tests** for both pathways, checking sign, units, zero-based population order, column-major layout, and atomic error status/output preservation.
- [ ] **Step 5: Build the new tests before implementation** using the targets below; expect failure on the missing new API declarations/bindings and confirm the failure is attributable to those new tests.
- [ ] **Step 6: Implement the smallest checked adapters** in the listed C++ sources and C/Fortran wrappers. Apply each bulk coefficient to its species/bin; apply each modal mass coefficient to all constituents of its mode and the number coefficient to that mode's number. Stage validated outputs before publication, and reuse the rate calculation for the simple velocity path without introducing a removal scheme.
- [ ] **Step 7: Build and run** `cmake --build build --target exaero_bulk_removal_test exaero_modal_removal_test exaero_backend_fortran_test exaero_modal_fortran_test && ctest --test-dir build -R 'ExaeroBulkRemoval|ExaeroModalRemoval|ExaeroBackendFortran|ExaeroModalFortran' --output-on-failure`; expected result is identical species/mode/moment-specific tendencies across C++, C, and Fortran with unchanged inputs.
- [ ] **Step 8: Commit** as `feat: return host-computed removal tendencies`.

### Task 4: Match modal resuspension and deposited-population rebinning

**Files:**
- Modify if needed: `src/src_impl/modal/ModalRemoval.cpp`
- Test: `tests/test_modal_removal.cpp`, `tests/test_modal_removal_math.cpp`
- Reference: matrix rows for CAM-SIMA modal `resuspension_resize` and `rebin_bulk_fluxes`

**Interfaces:**
- Consumes: existing route-based `redistributeModalResuspension` and budget-based `rebinModalDeposition` contracts.
- Produces: tested mapping from CAM-SIMA's material-to-larger-mode tendency handling to Ex-Aero's explicit host routes, plus matching lognormal bin/tail budgets where assumptions match.

- [ ] **Step 1: Add `test_cam_sima_resuspension_routes_preserve_material`** in `tests/test_modal_removal.cpp`: configure A routes small→medium and medium→large, set original A rates to 10 in small and 3 in medium, and set unrouted B to 4 in small; assert routed A is 10 in medium and 3 in large (the newly routed 10 is not routed twice), B remains 4 in small, total material rate is unchanged, and all number-rate slots remain zero.
- [ ] **Step 2: Add `test_cam_sima_deposition_rebin_preserves_budget`**: rebin one deposited lognormal mode across edges `[1e-9, 1e-7, 2e-7, 4e-7] m`; assert the sum of bins and both tails equals each input species budget and adding number budget does not alter material outputs.
- [ ] **Step 3: Build and run** `cmake --build build --target exaero_modal_removal_test && ctest --test-dir build -R 'ExaeroModalRemoval' --output-on-failure`; if the cases already pass, record parity and make no code change.
- [ ] **Step 4: Correct route ordering, moment calculation, or units only for a measured mismatch**; preserve Ex-Aero's explicit route graph and integrated-budget inputs.
- [ ] **Step 5: Rebuild and rerun** the same target; expected result is source-once routing and budget closure within the tolerance pinned in Task 1.
- [ ] **Step 6: Commit** as `fix: match modal removal reference budgets`.

### Task 5: Match activation and cloud/ice diagnostic behavior

**Files:**
- Modify if needed: `src/src_impl/CloudIceActivation.cpp`
- Modify if needed: `src/src_impl/modal/ModalCloudIce.cpp`
- Modify if needed: `src/src_impl/gocart/BulkCloudIce.cpp`
- Test: `tests/test_modal_cloud_ice.cpp`, `tests/test_bulk_cloud_ice.cpp`, `tests/test_cloud_ice_config.cpp`
- Reference: matrix rows for `actfracs` and modal ice-nucleation eligibility

**Interfaces:**
- Consumes: existing `modalActivationFractions`, `bulkActivationFractions`, `modalCloudIce`, and `bulkCloudIce` APIs.
- Produces: formula parity for number/mass activation fractions and configurable category/eligibility mapping, without aerosol-state transfer or an Ex-Aero nucleation-rate parameterization.

- [ ] **Step 1: Add `test_cam_sima_modal_activation_reference`** in `tests/test_modal_cloud_ice.cpp`: configure one mode with `r_g=1e-7 m`, `kappa=0.5`, `sigma_g=1.5`, and `T=280 K`; compute critical supersaturation from the configured modal constants, request `{Scrit/2,Scrit,min(2*Scrit,1)}`, and assert at `S=Scrit` that number fraction is `0.5` and mass fraction is `0.5*(1+erf(1.5*sqrt(2)*log(1.5)))`; assert both fractions are monotone over those three values.
- [ ] **Step 2: Add `test_configured_ice_category_weights_match_reference`**: in a 50/50 dry-mass DU/SU mode with dust-only category eligibility and fixed size fraction 0.4, assert dust type weight 0.5, size weight 0.4, eligible number `0.2*N`, and zero eligibility for sulfate; also assert empty and absent cloudborne status values.
- [ ] **Step 3: In the same tests, snapshot ambient/cloudborne inputs** and assert all values are unchanged after activation and cloud/ice diagnostics.
- [ ] **Step 4: Build and run** `cmake --build build --target exaero_modal_cloud_ice_test exaero_bulk_cloud_ice_test exaero_cloud_ice_config_test && ctest --test-dir build -R 'ExaeroModalCloudIce|ExaeroBulkCloudIce|ExaeroCloudIceConfig' --output-on-failure`; passing output establishes current parity, otherwise use the first differing assertion as the fix target.
- [ ] **Step 5: Correct a formula or role mapping only for a measured mismatch**; retain host-owned cloudborne buffers and configuration-driven categories; do not add transfer tendencies.
- [ ] **Step 6: Rebuild and rerun** the focused tests; expected result is matching formulas, explicit unsupported statuses, and unchanged input state.
- [ ] **Step 7: Commit** as `fix: align cloud and ice diagnostics with CAM-SIMA`.

### Task 6: Confirm bulk analogs and C/Fortran parity

**Files:**
- Modify only if needed: `src/src_impl/gocart/BulkRemoval.cpp`, `src/src_impl/gocart/BulkCloudIce.cpp`, `src/src_impl/gocart/BulkRemoval_C.cpp`, `src/src_impl/CloudIceCoupling_C.cpp`, `src/exaero/exaero_interface.f90`
- Test: `tests/test_bulk_removal.cpp`, `tests/test_bulk_cloud_ice.cpp`, `tests/test_cloud_ice_fortran.f90`, `tests/test_backend_fortran.f90`
- Modify: `docs/guides/cam-sima-parity.md`

**Interfaces:**
- Consumes: current mass-only bulk number conversion, physical-radius metadata, scavenging sizes, redistribution/deposition utilities, and cloud/ice helper contracts.
- Produces: a complete status for each applicable bulk analogue; no prognostic bulk number and no implementation of CAM-SIMA bulk stubs.

- [ ] **Step 1: Add `test_bulk_number_and_effective_radius_reference`** in `tests/test_bulk_removal.cpp`: for fixed mass, assert number follows the existing lognormal conversion at `sigma_g=1` and `sigma_g=1.5`, and effective radius follows the configured or derived value without consulting optics tables.
- [ ] **Step 2: Add `test_bulk_cloud_ice_fortran_layout_and_status`** in `tests/test_cloud_ice_fortran.f90`: assert zero-based species selection, column-major `(cell,level,population,12)` output ordering, particles/m³ units, and `Absent` status with a null cloudborne pointer.
- [ ] **Step 3: Build and run** `cmake --build build --target exaero_bulk_removal_test exaero_bulk_cloud_ice_test exaero_cloud_ice_fortran_test exaero_backend_fortran_test && ctest --test-dir build -R 'ExaeroBulkRemoval|ExaeroBulkCloudIce|ExaeroCloudIceFortran|ExaeroBackendFortran' --output-on-failure`; expected result is all analogous supported behaviors pass without introducing prognostic bulk number.
- [ ] **Step 4: Update the capability matrix** so every target bulk behavior has a final status and each intentional bulk/modal difference cites its representation reason and test.
- [ ] **Step 5: Commit** as `test: verify CAM-SIMA bulk analogs and wrappers`.

### Task 7: Full behavior-parity validation

**Files:**
- Test: all targets added above
- Modify: `docs/guides/cam-sima-parity.md`

- [ ] **Step 1: Run all behavior tests** with `ctest --test-dir build --output-on-failure`; expected result is the full supported test suite passes.
- [ ] **Step 2: Configure and run optional MAM4xx tests** in a provisioned build; expected result is the existing adapter tests pass and the parity work does not claim native MAM4xx microphysics not implemented here.
- [ ] **Step 3: Record the tested compiler, Kokkos, yaml-cpp, Fortran, and CAM-SIMA revisions** in the parity matrix; mark any environment-only validation gap explicitly.
- [ ] **Step 4: Commit** the final matrix as `docs: record aerosol parity validation`.
