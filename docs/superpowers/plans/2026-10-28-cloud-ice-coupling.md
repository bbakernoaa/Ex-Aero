# Cloud and Ice Coupling Helpers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add diagnostic cloud and ice coupling queries for GOCART bulk and generic modal aerosol, with host-owned cloudborne state and configurable ice rules.

**Architecture:** A dependency-free configuration/math header defines categories and size rules. Representation-specific host utilities produce matching per-population outputs without changing packed state or IAerosolPackage's virtual table. C and Fortran wrappers expose the same column-major contract.

**Tech Stack:** C++23, mdspan layout_left, Kokkos-backed existing package properties, yaml-cpp, C ABI, Fortran ISO_C_BINDING.

**Spec:** docs/superpowers/specs/2026-10-27-cloud-ice-coupling-design.md

## Global Constraints

- Mass is kg/m3, number is particles/m3, and radius is metres.
- GOCART fixed number uses lognormal_dg, lognormal_sigma, and dry_density; no optics/radiation-table query.
- Bulk dry number-median radius is lognormal_dg/2; bulk dry/wet scavenging radii use dry_particle_diameter and configured-kappa growth.
- Host owns the cloudborne packed-state buffer. No interstitial-to-cloudborne transfer or aerosol mutation.
- No change to IAerosolPackage's virtual table, packed interstitial state layout, or aggregate computeCCN.
- No implicit category inference from species names and no default CAM-SIMA ice threshold or 1/25 factor.
- All outputs publish atomically after shape, finite-value, range, and alias checks.

## Review Focus

- A role tag without a size rule returns NotConfigured and zero eligible number; test in Task 2.
- A cloudborne pointer omitted by the host returns zeros and Absent even when ambient aerosol exists; test in Task 2.
- Modal mixed dust/non-dust mass uses a dry-mass type fraction and does not reuse number as mass; test in Task 3.
- A point-size bulk distribution (sigma=1) has a deterministic cutoff result at the threshold; test in Task 1.
- Already-cloudborne particles do not enter liquid activation sums; test in Task 4.

---

### Task 1: Configurable physical roles and ice size rules

**Files:**
- Create: src/exaero/CloudIceConfig.hpp (IceCategory, IceSizeRuleKind, IceSizeRule, IcePopulationRule)
- Create: src/exaero/CloudIceMath.hpp (pure lognormal cutoff math)
- Modify: src/exaero/GocartConfig.hpp (per-species ice_categories and ice_rules)
- Modify: src/exaero/ModalConfig.hpp (per-species ice_categories; per-mode ice_rules)
- Modify: src/src_impl/gocart/GocartPackage.cpp and src/src_impl/gocart/GocartPackage.hpp (YAML/structured ingestion and retained validated rules)
- Modify: src/src_impl/modal/ModalPackage.cpp (YAML/structured ingestion and validation)
- Test: tests/test_cloud_ice_config.cpp

**Interfaces:**
- Produce IceCategory {Dust,Sulfate,BlackCarbon,StratosphericSulfate}, IceSizeRuleKind {None,All,FixedFraction,AboveDryDiameter}, IceSizeRule {kind,value}, and IcePopulationRule {category,normal,optional preexisting_ice,heterogeneous_freezing}.
- Produce double iceSizeWeight(const IceSizeRule&, double dry_number_median_diameter, double geometric_std_dev).
- A species carries zero or more explicit category tags; a population carries at most one rule per category. Rule/category mismatches, invalid enum values, nonfinite or out-of-range parameters, and duplicate tags/rules fail initialization.
- YAML uses ice_categories: [dust, sulfate, black_carbon, stratospheric_sulfate] on species and ice_rules: [{category: dust, kind: above_dry_diameter, value: 1.0e-7, preexisting_ice: {kind: all}, heterogeneous_freezing: true}] on a bulk species or modal mode. The names match the structured config fields.

- [ ] **Step 1: Write failing config/math tests.** Test structured and YAML parity, absent rules, duplicate/unknown tags, fixed fraction bounds, sigma=1 cutoff equality, and a lognormal 0.1 micrometre cutoff fraction against the erfc expression.
- [ ] **Step 2: Run targeted test and confirm failure.** Run the test_cloud_ice_config target; expect missing header/API or assertion failure.
- [ ] **Step 3: Add the types, pure math, parsers, and initialization validation.** Retain validated bulk rules in GocartPackage and modal rules in ModalConfig; preserve existing configurations.
- [ ] **Step 4: Run targeted and prior config tests.** Expect new config tests and existing GOCART/modal initialization tests to pass.
- [ ] **Step 5: Commit.** Use feat: configure cloud and ice aerosol roles.

### Task 2: Bulk cloud and ice candidate query

**Files:**
- Create: src/exaero/CloudIceCoupling.hpp (common output indices/status and both representation signatures)
- Create: src/src_impl/gocart/BulkCloudIce.cpp
- Modify: src/src_impl/gocart/GocartPackage.hpp (additive IBulkCloudIceMetadata)
- Modify: src/src_impl/CMakeLists.txt
- Test: tests/test_bulk_cloud_ice.cpp and tests/CMakeLists.txt

**Interfaces:**
- Produce void bulkCloudIce(IAerosolPackage&, const EnvironmentalStateView&, const View3D<const double>& ambient, const View3D<const double>* cloudborne, IceCategory, bool preexisting_ice, View4D<double>& ambient_out, View4D<double>& cloudborne_out, View3D<int>& ambient_status, View3D<int>& cloudborne_status).
- Output shape for each property slab is (cell,level,species,12): mass, number, dry_number_median_radius, dry_scavenging_radius, wet_scavenging_radius, lower_radius, upper_radius, effective_radius, type_weight, size_weight, eligible_number, hetero_freezing_flag. Status shape is (cell,level,species): Absent=-2, Empty=-1, NotConfigured=0, Configured=1. Missing bounds are zero; bulkReferenceRadii separately reports bound availability.
- IBulkCloudIceMetadata returns immutable validated tags/rules by species and reuses IBulkRemovalMetadata for size. A cloudborne bulk mass field uses the same fixed mass-to-number conversion as ambient.

- [ ] **Step 1: Write failing bulk tests.** Test ambient number parity with bulkFixedNumber, distinct number-median/scavenging radii, matching and nonmatching type weights, no-rule status, ice cutoff, explicit heterogeneous-freezing flag, host cloudborne input, omitted cloudborne input, and atomic rejection of malformed arrays.
- [ ] **Step 2: Run targeted test and confirm failure.** Expect missing bulkCloudIce/API.
- [ ] **Step 3: Implement bulkCloudIce and metadata.** Reuse removal math and radius semantics; stage results before copying to outputs.
- [ ] **Step 4: Run targeted and bulk removal tests.** Expect both to pass.
- [ ] **Step 5: Commit.** Use feat: expose bulk cloud and ice candidates.

### Task 3: Generic modal cloud and ice candidate query

**Files:**
- Create: src/src_impl/modal/ModalCloudIce.cpp
- Modify: src/src_impl/CMakeLists.txt
- Test: tests/test_modal_cloud_ice.cpp and tests/CMakeLists.txt

**Interfaces:**
- Produce void modalCloudIce(IAerosolPackage&, const EnvironmentalStateView&, const View3D<const double>& ambient, const View3D<const double>* cloudborne, IceCategory, bool preexisting_ice, View4D<double>& ambient_out, View4D<double>& cloudborne_out, View3D<int>& ambient_status, View3D<int>& cloudborne_status).
- Same output indices/status as Task 2, but population axis is mode. Number is the packed prognostic number; type weight is category-tagged dry mass divided by total mode dry mass. Radius is the current modal lognormal number-median radius. The host cloudborne buffer has exactly the same packed field order.

- [ ] **Step 1: Write failing modal tests.** Test pure and mixed-composition modes, empty mode, current number/radius parity with modalProperties, cutoff-tail rule and preexisting-ice override, cloudborne-only input, disallowed species, and atomic shape/NaN rejection.
- [ ] **Step 2: Run targeted test and confirm failure.** Expect missing modalCloudIce/API.
- [ ] **Step 3: Implement modalCloudIce.** Keep the existing interstitial layout and calculate from the supplied population snapshot without applying tendencies.
- [ ] **Step 4: Run targeted, modal removal, and modal core tests.** Expect pass.
- [ ] **Step 5: Commit.** Use feat: expose modal cloud and ice candidates.

### Task 4: Liquid activation fractions, cross-language bindings, and guide

**Files:**
- Modify: src/exaero/CloudIceCoupling.hpp
- Create: src/src_impl/gocart/BulkCloudIce_C.cpp and src/src_impl/modal/ModalCloudIce_C.cpp
- Create: src/exaero/CloudIceCoupling_C.h
- Modify: src/exaero/exaero_interface.f90 and src/src_impl/CMakeLists.txt
- Create: docs/guides/cloud-ice-coupling.md
- Modify: docs/index.md and docs/guides/aerosol-backends.md
- Test: tests/test_cloud_ice_activation.cpp, tests/test_cloud_ice_c.cpp, and tests/standalone_public_check.cpp

**Interfaces:**
- Produce void bulkActivationFractions(IAerosolPackage&, const EnvironmentalStateView&, const View3D<const double>& ambient, const View1D<const double>& supersaturation, View4D<double>& fractions) and modalActivationFractions with the same signature. Fractions shape is (cell,level,population,2*nsupersaturation), number then mass for each supersaturation.
- Bulk fractions are the existing critical-supersaturation 0/1 decision; modal number uses the existing lognormal erfc CCN fraction and mass uses the corresponding third moment. Only interstitial aerosol is an activation input.
- C wrappers expose both cloud/ice and activation functions with zero-based indices, nullable cloudborne state, fixed property/status constants, and 256-byte error buffer; Fortran ISO_C_BINDING mirrors them.
- Name the C entry points exaero_bulk_cloud_ice, exaero_modal_cloud_ice, exaero_bulk_activation_fractions, and exaero_modal_activation_fractions. Each accepts the package, nc/nl/nf extents, input pointers, output pointers, errmsg, and errflg; cloud/ice also accepts category and preexisting_ice integer flags.

- [ ] **Step 1: Write failing activation and C ABI tests.** Test sum(number*fraction) against computeCCN at several supersaturations, modal third-moment mass fraction, monotonicity, nonhygroscopic/empty behavior, cloudborne exclusion, C extent/error handling, and C/Fortran declarations.
- [ ] **Step 2: Run targeted test and confirm failure.** Expect missing activation functions or C symbols.
- [ ] **Step 3: Implement activation fractions, C wrappers, Fortran declarations, and guide.** Preserve existing computeCCN and avoid host state mutation.
- [ ] **Step 4: Run targeted tests, standalone public-header check, all runnable regression tests, and formatting/static checks.** Expect pass; record full-build dependencies if unavailable.
- [ ] **Step 5: Commit.** Use feat: expose cloud and ice coupling across C and Fortran.

## Completion

Run git diff --check, the project test suite when the toolchain is present, standalone public/API checks otherwise, and inspect the final branch diff against the design. Report numerical semantics, tests run, and any unavailable dependencies.
