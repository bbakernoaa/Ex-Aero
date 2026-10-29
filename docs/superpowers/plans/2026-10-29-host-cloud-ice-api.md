# Common Cloud/Ice Coupling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose common cloud/ice and activation calls for bulk and generic modal populations.

**Architecture:** Slice 2 of the [roadmap](2026-10-29-host-coupling-roadmap.md) consumes slice 1 metadata and dispatches through existing cloud/ice implementations. Numerical kernels and the existing 12-property contract remain authoritative.

**Tech Stack:** C++23 host mdspan, C ABI, Fortran bind(c), existing CMake/CTest profiles.

**Spec:** [Common coupling](../specs/2026-10-29-host-coupling-design.md).

## Global Constraints

- Finish metadata/removal slice 1 first; use its population count and capabilities.
- Column-major host arrays; unchanged zero-based indices, units, category codes and statuses.
- No activation transfer, cloudborne state creation, ice nucleation integrator or accepted update.
- Preserve every existing public entry point; no package vtable change.

## Review Focus

- Bulk output population count is bins; modal output count is modes, never constituent count.
- Missing cloudborne state returns `Absent`, distinct from a present empty snapshot.
- Ambient and cloudborne snapshots can differ in composition and number.
- Activation mass and number fractions are separate and may differ strongly.
- Invalid categories, unsupported packages, aliasing and late bad inputs preserve every output.

### Task 1: Add common C++ dispatch

**Files:** Modify `src/exaero/CloudIceCoupling.hpp`, create `src/src_impl/CloudIceCoupling.cpp`; modify `src/src_impl/CMakeLists.txt`, `tests/test_bulk_cloud_ice.cpp`, `tests/test_modal_cloud_ice.cpp`.

**Interfaces:**
- `cloudIce(package,env,ambient,cloudborne_or_null,category,preexisting_ice,ambient_out,cloudborne_out,ambient_status,cloudborne_status)` has exactly the existing `bulkCloudIce` argument types and population-based shapes.
- `activationFractions(package,env,ambient,supersaturation,fractions)` has exactly the existing bulk/modal activation argument types; output `(cell,level,P,2*num_S)`.

- [ ] **Step 1: Add failing common-call tests** comparing every returned property/status with the matching existing representation-specific call for two populations and multiple cells/levels. Use different ambient/cloudborne snapshots; verify null cloudborne gives `Absent`, while supplied zero cloudborne gives `Empty`.
- [ ] **Step 2: Build affected targets and observe missing-common-function failure.**
- [ ] **Step 3: Implement dispatch** using slice 1 support queries and existing kernels; advertise `CloudIce` and `Activation` only for supported initialized packages.
- [ ] **Step 4: Add activation and atomic-error tests.** At a configured mode's critical supersaturation, require number fraction 0.5 and the existing independent mass-tail reference. Cover zero-kappa behavior, bad population extents, invalid category/NaN input, unsupported packages and output overlap.
- [ ] **Step 5: Run `ExaeroBulkCloudIce`, `ExaeroModalCloudIce`, and `ExaeroHostCoupling`; commit** as `feat: unify aerosol cloud and activation helpers`.

### Task 2: Deliver bindings and host documentation

**Files:** Modify `src/exaero/CloudIceCoupling_C.h`, `src/src_impl/CloudIceCoupling_C.cpp`, `src/exaero/exaero_interface.f90`, `tests/test_cloud_ice_fortran.f90`, `tests/standalone_public_check.cpp`, `docs/guides/host-coupling.md`.

**Interfaces:** `exaero_cloud_ice` and `exaero_activation_fractions` take explicit `nc,nl,nf,np` counts, then the corresponding existing inputs/outputs and error buffers; activation also takes `nss`. Cloud-state pointers are nullable; caller-provided `np` must agree with metadata.

- [ ] **Step 1: Add failing C/Fortran common-call tests** with two populations and multiple levels; assert slab ordering, category/status codes, zero-based indexing, independent activation moments and unchanged state.
- [ ] **Step 2: Observe missing bindings, then implement wrappers** with the existing 256-byte message/error flag contract and null-cloud handling.
- [ ] **Step 3: Run the cloud/ice Fortran target with runtime checks**, including wrong `np`, missing cloudborne, invalid inputs and all output sentinels preserved on failure.
- [ ] **Step 4: Document common examples and capability checks; add standalone header coverage.**
- [ ] **Step 5: Build/run both configured CTest profiles, run `git diff --check`, and commit** as `feat: expose common cloud coupling to C and Fortran`. Record actual results without claiming a complete cloud microphysics scheme.
