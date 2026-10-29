# Host-Coupling Metadata and Removal Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Establish common population metadata and operation queries, then add representation-aware rate and velocity removal, including constituent-specific modal mass loss.

**Architecture:** This is slice 1 of the [host-coupling roadmap](2026-10-29-host-coupling-roadmap.md). Add metadata and capability queries without changing the package vtable, then dispatch common rate/velocity calls through existing implementations and a constituent-specific modal conversion. Preserve the mode-shared convenience helper. Every operation includes C and Fortran bindings.

**Tech Stack:** C++23, `std::mdspan` backport, existing C ABI wrappers, Fortran `bind(c)`, Catch2 2, CMake/CTest.

**Spec:** [Common coupling](../specs/2026-10-29-host-coupling-design.md) and [removal details](../specs/2026-10-29-representation-aware-removal-design.md).

## Global Constraints

- Rates are finite, nonnegative inverse seconds; concentrations are finite and nonnegative.
- State and mass-tendency layouts remain host-memory column-major.
- Bulk number remains diagnostic and follows its species/bin mass loss.
- Modal number tendencies remain independent and mode-level; modal mass rates apply per constituent field.
- Modal mass-rate values at packed number-field slots must be exactly zero.
- Outputs are staged and published atomically, disjoint from inputs and each other; no state is mutated.
- The common entry point derives representation from `IAerosolPackage::representation()`.
- Existing representation-specific rate and velocity helpers remain behaviorally unchanged.
- Population count `P` means bulk species/bin count or generic modal mode count. Metadata maps every allocated constituent mass slot, including its permission, and the number field where present.
- Metadata is valid for the current initialized configuration; re-query after reinitialization. No name parsing, global registry, new vtable member or `StateField` layout change.
- Operation availability checks actual package support; unimplemented future facade operations and unsupported native MAM4xx helpers report false.

## Review Focus

- A bulk package passed to the common dispatcher must use species/bin shapes and derive number loss consistently from mass loss.
- A modal package must use one mass coefficient per packed constituent field and one number coefficient per mode.
- Missing modal number-rate data or a nonzero coefficient in a modal number field must fail atomically.
- Unsupported package representation or invalid package initialization must fail before publishing outputs.
- Aliasing, invalid values, late invalid entries, and overflow must preserve both output buffers.
- Metadata must distinguish diagnosed bulk number, prognostic modal number, prohibited constituent slots and unsupported helper operations; hosts must not reconstruct offsets.

---

### Task 0: Establish population metadata and capability queries

**Files:**
- Create: `src/exaero/HostCoupling.hpp`, `src/src_impl/HostCoupling.cpp`
- Create: `src/exaero/HostCoupling_C.h`, `src/src_impl/HostCoupling_C.cpp`
- Create: `tests/test_host_coupling.cpp`, `tests/test_host_coupling_fortran.f90`
- Modify: `src/src_impl/CMakeLists.txt`, `tests/CMakeLists.txt`
- Modify: `src/exaero/exaero_interface.f90`, `tests/standalone_public_check.cpp`

**Interfaces:**
- Produce `NumberKind { Diagnosed=0, Prognostic=1 }`, `PopulationInfo`, `ConstituentInfo` and `HostCouplingInfo` with the exact fields in the common design. Use integer indices, `std::string` names, descriptor vectors and `bool permitted`.
- Produce `HostCouplingInfo hostCouplingInfo(IAerosolPackage&)` and `bool supportsHostOperation(IAerosolPackage&, HostOperation)`.
- `HostOperation` enumerates `PopulationProperties`, `CloudIce`, `Activation`, `RemovalRates`, `RemovalVelocity`, `Resuspension`, `DepositionRebin`, `InspectNumberBounds`, `EnforceNumberBounds`. Each flag becomes true only when that callable facade is delivered; bounds flags describe the existing explicit modal helpers.
- C/Fortran provide `exaero_host_coupling_counts`, `exaero_host_population_info`, `exaero_host_constituent_info`, and `exaero_supports_host_operation`: counts first, indexed descriptors next, caller-owned name buffers with explicit byte capacities. Copy enum codes to Fortran constants; insufficient capacity and invalid indices preserve output sentinels and set an error.

- [x] **Step 1: Add failing mapping tests.** Two bulk bins map to populations 0/1 and mass fields 0/1, with number kind diagnosed and number field -1. Two modal modes with two configured materials map mass fields `{0,1,3,4}` and number fields `{2,5}`. Restrict one mode to one material; retain its prohibited allocated mass-field descriptor with `permitted=false`.
- [x] **Step 2: Build the new C++ metadata target and observe the missing-interface failure.**
- [x] **Step 3: Implement metadata from validated configuration and existing state-field metadata.** Reject inconsistent initialization, never split or parse field labels. An initialized legacy package without the required extension returns false for helper support; unsupported metadata queries fail explicitly.
- [x] **Step 4: Add C and Fortran descriptor round-trip tests.** Include reordered configuration, reinitialization, diagnosed/prognostic number distinction, unknown operation codes, invalid indices, too-small name buffers, unsupported native MAM4xx and uninitialized packages. Known unsupported operation queries succeed with false; invalid query arguments set an error.
- [x] **Step 5: Implement bindings and standalone-header coverage; run `ExaeroHostCoupling` and `ExaeroHostCouplingFortran`.** Require correct zero-based indices and no unintended output changes on error.
- [x] **Step 6: Commit** as `feat: expose aerosol population and coupling metadata`.

### Task 1: Add constituent-specific modal rate conversion

**Files:**
- Modify: `src/exaero/ModalRemoval.hpp`
- Modify: `src/src_impl/modal/ModalRemoval.cpp`
- Test: `tests/test_modal_removal.cpp`

**Interfaces:**
- Consumes: `IAerosolPackage`, packed modal state `(cell,level,field)`, configured mode/species field mapping.
- Produces: `modalSpeciesRemovalRateTendencies(package,state,mass_rate_by_field,number_rate_by_mode,mass_tendency,number_tendency)`; mass rates and tendencies match state extent, number rates and tendencies use `(cell,level,mode)`.

- [x] **Step 1: Add the failing behavior test** `modal_species_rate_tendencies_allow_distinct_constituent_losses` to `tests/test_modal_removal.cpp`. In one mode with two populated species masses `2e-9` and `3e-9 kg/m3`, provide rates `5e-4` and `2e-4 /s`, a number `5e7 /m3`, and number rate `1e-4 /s`; assert tendencies `-1e-12`, `-6e-13 kg/m3/s`, and `-5e3 particles/m3/s`. Assert every number slot in mass output is zero.
- [x] **Step 2: Run the focused test and verify it fails** because the species-rate helper is missing.
- [x] **Step 3: Implement `modalSpeciesRemovalRateTendencies`** using the existing validation, staging, overflow, and atomic-publication rules. Require finite nonnegative coefficient values and exact zero in packed number-field mass-rate slots; keep distinct configured constituent rates independent.
- [x] **Step 4: Add and run invalid-input tests** for a nonzero number-slot mass rate, wrong-shaped rates, negative/NaN late entries, aliasing, and output sentinels unchanged on failure.
- [x] **Step 5: Run `ctest --test-dir build -R '^ExaeroModalRemoval$' --output-on-failure`** and require the focused target to pass.
- [x] **Step 6: Commit** as `feat: support constituent-specific modal removal rates`.

### Task 2: Add the common C++ representation-aware API

**Files:**
- Create: `src/exaero/Removal.hpp`
- Create: `src/src_impl/Removal.cpp`
- Modify: `src/src_impl/CMakeLists.txt`
- Modify: `tests/test_bulk_removal.cpp`
- Modify: `tests/test_modal_removal.cpp`

**Interfaces:**
- Consumes: Task 0 metadata, `bulkRemovalRateTendencies` and Task 1's `modalSpeciesRemovalRateTendencies`.
- Produces: `removalRateTendencies(package,state,mass_rate,number_rate_or_null,mass_tendency,number_tendency)`.
- Bulk: mass rate/output and diagnostic number output use `(cell,level,species/bin)`; `number_rate_or_null` must be null; number tendency follows the existing fixed-distribution mass loss.
- Modal: mass rate/output use packed `(cell,level,field)`; modal number rate/tendency use `(cell,level,mode)`; a non-null number-rate view is required.
- Produces: `removalVelocityTendencies(package,thickness,state,velocity,mass_tendency,number_tendency)`; velocity and number output use `(cell,level,P)`, mass output uses state layout. Delegate to existing bulk/modal velocity helpers.
- Both mass-rate and velocity outputs use the metadata field/population mapping; retain separate number output, with zero number slots in packed modal mass output.

- [x] **Step 1: Add failing dispatcher tests** to the existing bulk and modal removal executables. Assert bulk species/bin rates produce matching mass and fixed-distribution number losses, and distinct modal constituent rates plus independent mode-level number rate produce the Task 1 results.
- [x] **Step 2: Build the affected removal targets and verify compilation fails** because the public API/header is absent.
- [x] **Step 3: Implement `removalRateTendencies`** in `src/src_impl/Removal.cpp`, dispatching only on `package.representation()`, validating the representation-specific optional number-rate contract, and delegating to the existing/tested conversion paths.
- [x] **Step 4: Add failing velocity-dispatch tests** for two populations with different velocities and multiple cells/levels. Assert `v=0.01 m/s`, `dz=10 m` gives `0.001 /s`; zero velocity disables only its population. Implement the common velocity dispatcher and enable both removal operation flags.
- [x] **Step 5: Run dispatcher error tests and `ctest --test-dir build -R '^(ExaeroHostCoupling|ExaeroBulkRemoval|ExaeroModalRemoval)$' --output-on-failure`.** Cover a modal package without number rates, bulk with number-rate input, wrong population count, missing extension support, uninitialized package, invalid thickness and aliasing; assert both outputs remain unchanged.
- [x] **Step 6: Commit** as `feat: add representation-aware removal API`.

### Task 3: Expose the common API through C and Fortran

**Files:**
- Create: `src/exaero/Removal_C.h`
- Create: `src/src_impl/Removal_C.cpp`
- Modify: `src/exaero/exaero_interface.f90`
- Modify: `tests/test_removal_fortran.f90`
- Modify: `tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Task 2's C++ `removalRateTendencies`.
- Produces: `exaero_removal_rate_tendencies(pkg,nc,nl,nf,np,state,mass_rate,number_rate,mass_tendency,number_tendency,errmsg,errflg)`.
- Produces: `exaero_removal_velocity_tendencies(pkg,nc,nl,nf,np,state,velocity,thickness,mass_tendency,number_tendency,errmsg,errflg)`.
- The C `number_rate` pointer is null for bulk and required for modal. `nf` is state-field extent; `np` is metadata population count for both representations. Both number outputs use `(nc,nl,np)`.
- The Fortran `bind(c)` declaration uses `type(c_ptr), value` for the optional number-rate pointer so bulk callers pass `c_null_ptr` and modal callers pass a `c_loc` pointer to contiguous mode rates.

- [x] **Step 1: Add failing C/Fortran calls** for rates and velocities, with bulk using `c_null_ptr` number rates and modal using `c_loc` of separate per-mode number-rate storage; assert output values, units, and unchanged input state.
- [x] **Step 2: Build the Fortran test and verify it fails** because the binding is not declared or linked.
- [x] **Step 3: Implement the C wrapper** with null-handle, dimension, and pointer checks; translate exceptions to the established 256-byte error message and error flag without touching outputs on failure.
- [x] **Step 4: Add both `bind(c)` interfaces and run `ExaeroRemovalFortran`** with GNU Fortran runtime checks enabled; require bulk/modal values, mismatched population counts, missing modal number rates and atomic error cases to pass. Add new public headers to `tests/standalone_public_check.cpp`.
- [x] **Step 5: Commit** as `feat: expose removal dispatcher to C and Fortran`.

### Task 4: Document host usage and verify both build profiles

**Files:**
- Modify: `docs/guides/host-removal.md`
- Modify: `docs/guides/cam-sima-parity.md`
- Modify: `docs/superpowers/specs/2026-10-29-representation-aware-removal-design.md`
- Modify: this plan

**Interfaces:**
- Consumes: public C++, C, and Fortran declarations from Tasks 1–3.
- Produces: documented dispatcher shape/unit contracts and the difference between per-constituent modal mass removal and mode-level number removal.

- [x] **Step 1: Update the host guide** with metadata-driven sizing, capability checks, rate and velocity examples for each representation, and the Fortran optional pointer contract. Add `docs/guides/host-coupling.md` as the facade overview linking the later slices.
- [x] **Step 2: Update the CAM crosswalk** to state that Ex-Aero now permits constituent-specific modal mass tendencies while the host still supplies rates and mode-level number removal; do not claim wet-deposition-scheme parity.
- [x] **Step 3: Build and run the generic profile** with the existing configured Debug build and require all CTest tests to pass.
- [x] **Step 4: Build and run the optional MAM4xx profile** with its existing configured Debug build and require all CTest tests to pass.
- [x] **Step 5: Run `git diff --check` and inspect all public declarations/bindings** for matching dimensions and pointer semantics.
- [x] **Step 6: Commit** as `docs: describe representation-aware removal rates`.
