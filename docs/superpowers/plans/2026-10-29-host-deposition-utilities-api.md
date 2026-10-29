# Common Deposition Utilities Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Supply common resuspension-routing and deposited-budget rebinning calls with explicit field/species selectors.

**Architecture:** Slice 4 uses slice 1 field mapping and retains existing distribution geometry and conservation rules. A common routing layer accepts explicit mass-field routes for either representation; a configuration query translates existing modal routes. Rebinning dispatches to the existing bulk/modal kernels.

**Tech Stack:** C++23 host mdspan and span, C ABI, Fortran bind(c), CMake/CTest.

**Spec:** [Common coupling](../specs/2026-10-29-host-coupling-design.md).

## Global Constraints

- Rates are material returns in kg/m3/s; deposited budgets are kg/m2 plus modal particle number/m2 in its number fields.
- Routes are transformations of host-supplied budgets, not precipitation or resuspension physics.
- Route original sources once; no inferred particle-number return and no implicit number correction.
- Keep existing bulk/modal entry points, configured routes and distribution formulas unchanged.
- Rebin uses dry diameter edges in m, with explicit below/above tails and material-budget closure.

## Review Focus

- Route indices identify state mass fields; number fields and prohibited destinations are invalid.
- A chain A→B→C moves each original input once, without rerouting material newly received at B.
- Duplicate sources, cycles, invalid indices and modal cross-material routes fail atomically.
- Bulk selects one species/bin; modal selects one material across modes. Species indices are not state-field or population indices.
- Deposition mixture geometry uses all relevant mode constituents and number; bins plus tails close the selected material budget.

### Task 1: Add common resuspension routing

**Files:** Create `src/exaero/DepositionCoupling.hpp`, `src/src_impl/DepositionCoupling.cpp`, `tests/test_deposition_coupling.cpp`; modify source/test CMake files and `src/src_impl/HostCoupling.cpp`.

**Interfaces:**
- `ResuspensionFieldRoute { std::size_t source_field, destination_field; }`.
- `configuredResuspensionRoutes(IAerosolPackage&) -> std::vector<ResuspensionFieldRoute>` translates modal configuration routes; bulk returns an empty list, requiring host-supplied routes.
- `redistributeResuspension(package,const View3D<const double>& mass_rates,std::span<const ResuspensionFieldRoute>,View3D<double>& routed_rates)` uses state layout. Modal number slots must be zero; both endpoints must permit the material and identify the same modal material. Bulk routes follow explicit host species/bin mapping without inferring chemical identity from names.

- [ ] **Step 1: Add failing conservation tests** with original returns 10 at A, 3 at B, 4 at unrouted D and routes A→B, B→C. Require B=10, C=3, D=4, total 17 and zero modal number slots. Compare configured modal-route results with the existing helper and bulk explicit routes with the existing bulk helper.
- [ ] **Step 2: Observe missing interface failure; implement staged routing** using metadata. Reject duplicate source routes, self routes, cycles, nonmass endpoints, cross-material modal routes, forbidden constituents, invalid inputs, aliasing and overflow.
- [ ] **Step 3: Test all invalid route classes** and late overflow with output sentinels unchanged; enable `Resuspension` support after implementation.
- [ ] **Step 4: Run `ExaeroDepositionCoupling` and existing bulk/modal removal targets; commit** as `feat: unify aerosol resuspension routing`.

### Task 2: Add common deposition rebinning

**Files:** Modify `src/exaero/DepositionCoupling.hpp`, `src/src_impl/DepositionCoupling.cpp`, `tests/test_deposition_coupling.cpp`, `src/src_impl/HostCoupling.cpp`.

**Interfaces:** `rebinDeposition(package,const View3D<const double>& deposited_budget,std::size_t species_index,const View1D<const double>& diameter_edges,View3D<double>& bin_mass,View3D<double>& tail_mass)`. Outputs use existing bin and two-tail shapes/units; metadata declares species selector scope.

- [ ] **Step 1: Add failing bulk/modal dispatch comparisons** and an independent single-material lognormal CDF case. Require bins plus tails equal the selected input material budget and zero-budget outputs equal zero.
- [ ] **Step 2: Observe missing interface failure; dispatch to existing kernels**, preserving bulk fixed geometry and modal mixture-volume/mass-number geometry. Enable `DepositionRebin` support.
- [ ] **Step 3: Test mixture sensitivity and selector errors.** Changing modal deposited number at fixed mass must change fractions; scaling both mass and number together preserves fractions. Wrong species indices, nonincreasing edges, invalid budgets and aliased outputs fail atomically.
- [ ] **Step 4: Run deposition and removal targets; commit** as `feat: unify aerosol deposition rebinning`.

### Task 3: Deliver C/Fortran interfaces and verification

**Files:** Create `src/exaero/DepositionCoupling_C.h`, `src/src_impl/DepositionCoupling_C.cpp`, `tests/test_deposition_coupling_fortran.f90`; modify source/test CMake files, `src/exaero/exaero_interface.f90`, `tests/standalone_public_check.cpp`, `docs/guides/host-coupling.md`, `docs/guides/cam-sima-parity.md`.

**Interfaces:** `exaero_resuspension_routes` uses count/capacity-based caller buffers for configured source/destination fields; `exaero_redistribute_resuspension` accepts `nc,nl,nf,nroute` and explicit source/destination arrays; `exaero_rebin_deposition` accepts `nc,nl,nf,species_index,nbin` and the existing budget/edge/output buffers. All use the established package/error arguments and zero-based indices.

- [ ] **Step 1: Add failing C/Fortran tests** for both representations, route conservation, separate number slots and selected-species bin/tail closure.
- [ ] **Step 2: Observe missing bindings, implement wrappers and test** insufficient route capacity, wrong field/species indices and invalid edges; preserve all outputs on error.
- [ ] **Step 3: Document explicit routing and selector scope**, emphasizing host-owned return budgets, and update the capability crosswalk without claiming complete wet-deposition physics. Add standalone header coverage.
- [ ] **Step 4: Build/run both configured CTest profiles, run `git diff --check`, and commit** as `feat: expose common deposition utilities across languages`.
