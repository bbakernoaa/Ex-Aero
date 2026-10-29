# Common Population Properties Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Supply hosts with explicitly defined bulk/modal physical properties on a common population axis.

**Architecture:** Slice 3 consumes slice 1 metadata, bulk physical metadata and generic modal properties. Reuse current sizing/growth formulas and expose availability per value; preserve distinct size definitions instead of substituting one for another.

**Tech Stack:** C++23 host mdspan, C ABI, Fortran bind(c), CMake/CTest.

**Spec:** [Common coupling](../specs/2026-10-29-host-coupling-design.md).

## Global Constraints

- All properties use SI units except dimensionless kappa and lognormal width.
- Population order comes from metadata; no dependence on optics/radiation tables.
- Bulk number is diagnostic; modal number comes from the mode's prognostic field.
- Preserve current growth laws/RH caps and physical distribution assumptions; derivatives remain an explicit modal extension.
- No state mutation; values, availability and population-status outputs publish atomically.

## Review Focus

- Caller-anchored scavenging diameter and lognormal median may differ in bulk.
- Missing physical radius bounds are unavailable, rather than physical zero bounds.
- Empty population status is separate from configured/derived property availability.
- Changing modal composition changes kappa and density without changing the metadata mapping.
- Invalid state, environment, property buffers or unsupported package preserves all outputs.

### Task 1: Define and implement common physical properties

**Files:** Create `src/exaero/PopulationProperties.hpp`, `src/src_impl/PopulationProperties.cpp`, `tests/test_population_properties.cpp`; modify source/test CMake files and `src/src_impl/HostCoupling.cpp`.

**Interfaces:**
- `populationProperties(IAerosolPackage&,const EnvironmentalStateView&,const View3D<const double>&,View4D<double>& values,View4D<int>& availability,View3D<int>& population_status)`.
- Values and availability use `(cell,level,P,COUNT)`; population status is `(cell,level,P)` with `Empty`/`Populated` codes.
- Properties, in stable order: `Mass`, `Number`, `DryNumberMedianRadius`, `WetNumberMedianRadius`, `DryScavengingDiameter`, `WetScavengingDiameter`, `DryVolume`, `DrySurfaceArea`, `Kappa`, `DryDensity`, `WetDensity`, `GeometricStdDev`, `LowerDryRadius`, `UpperDryRadius`, `EffectiveDryRadius`.
- Availability codes: `Unavailable`, `Configured`, `Derived`; unavailable numeric slots are zero and must be interpreted with status. Publish these codes and property-unit documentation identically across languages.

- [ ] **Step 1: Add failing analytical tests.** For one matched single-material bulk/modal distribution, verify `V=M/rho`, `N=V/((4*pi/3)*rg^3*exp(4.5*ln(sigma)^2))`, dry area `4*pi*N*rg^2*exp(2*ln(sigma)^2)` and effective dry radius `rg*exp(2.5*ln(sigma)^2)`. Compare bulk median radius from `lognormal_dg/2` against a deliberately different `dry_particle_diameter/2` scavenging radius.
- [ ] **Step 2: Observe missing API failure, then implement** by reusing bulk metadata, fixed number/reference radii and modal property calculations. Wet median/scavenging sizes use their respective dry definitions and the current representation's growth model; wet density includes current aerosol-water volume and the existing water-density constant.
- [ ] **Step 3: Add composition/availability cases.** At fixed dry volume/number, constituent kappas `{0.6,0}` with volume fractions `{0.75,0.25}` yield kappa 0.45; swapping fractions yields 0.15. Missing bounds report unavailable. For empty populations, return zero extensive quantities and configured reference geometry where defined; composition-dependent empty-modal density/kappa remain unavailable. Do not divide zero by zero.
- [ ] **Step 4: Add atomic-error and capability cases** for bad state, NaN RH, wrong population/property dimensions, aliased outputs, uninitialized and unsupported packages; enable `PopulationProperties` support only after implementation.
- [ ] **Step 5: Run `ExaeroPopulationProperties` plus existing bulk/modal/cloud tests and commit** as `feat: add common aerosol population properties`.

### Task 2: Deliver bindings and documentation

**Files:** Create `src/exaero/PopulationProperties_C.h`, `src/src_impl/PopulationProperties_C.cpp`, `tests/test_population_properties_fortran.f90`; modify source/test CMake files, `src/exaero/exaero_interface.f90`, `tests/standalone_public_check.cpp`, `docs/guides/host-coupling.md`.

**Interfaces:** `exaero_population_properties(pkg,nc,nl,nf,np,temperature,rh,state,values,availability,population_status,errmsg,errflg)`; property-axis count comes from a public constant with identical C/Fortran values.

- [ ] **Step 1: Add failing C/Fortran tests** for bulk and modal over multiple cells/levels; verify property slab order, units, availability and empty status independently.
- [ ] **Step 2: Observe missing binding failure, implement wrappers**, then verify wrong dimensions, invalid pointers and late invalid inputs preserve every output buffer.
- [ ] **Step 3: Document exact size definitions, property provenance and empty behavior**, with a metadata-driven host example and modal-only bounds capability query. Add standalone header coverage.
- [ ] **Step 4: Build/run both configured CTest profiles, run `git diff --check`, and commit** as `feat: expose common population properties across languages`.
