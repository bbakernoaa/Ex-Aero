# Cloud and ice coupling helpers for GOCART bulk and generic modal aerosol

**Status:** proposed; scientific and interface review pending.

## Intent

Give host cloud and ice microphysics schemes the aerosol quantities and eligibility information they need for both Ex-Aero representations. GOCART bulk should offer the useful cloud/ice queries exposed by CAM-SIMA's bulk aerosol interface, and generic modal should offer equivalent queries with its prognostic particle number and internally mixed composition. Ex-Aero remains the aerosol property and diagnostic boundary; the host chooses nucleation, activation, freezing, and transfer rates and owns cloudborne buffers. This preserves the existing division with MKPP, where chemistry and precursor formation occur upstream.

## Goals and constraints

- Expose ambient and optional host-supplied cloudborne mass and number, CCN activation fractions, dry and wet size, ice candidate number, aerosol-type and size weights, and heterogeneous-freezing eligibility for each configured bulk species or modal mode.
- Keep the outputs diagnostic. No aerosol state updates, interstitial-to-cloudborne transfer, cloud droplet or ice crystal tendencies, or new nucleation parameterization occur in Ex-Aero.
- Use Ex-Aero's existing concentration units: kg/m3, particles/m3, and metres. A host that needs CAM-SIMA's particles/cm3 divides by 1e6; no air-density conversion is needed because Ex-Aero state is already per air volume.
- Preserve existing GOCART mass-to-number conversion from configured lognormal diameter, width, and density. The physical lower/upper/effective radii remain size metadata and are not read from radiation or GEOSmie tables. Modal number comes from its prognostic number field.
- Preserve the existing aggregate `computeCCN` result. New per-population activation diagnostics must sum back to that result within numerical tolerance.
- Leave the `IAerosolPackage` virtual table and the packed interstitial state layout unchanged. Add optional C++ surfaces, C ABI wrappers, and Fortran bindings consistent with the removal helpers.
- Validate finite, nonnegative state and finite environmental inputs, exact extents, valid category/rule configuration, and output ranges before publishing results. Diagnostics should be deterministic and should not require global gathering or MPI communication.

## Approaches considered

1. **Configurable diagnostic helpers — selected.** Explicit aerosol roles and size rules determine ice eligibility; per-population outputs use the existing Ex-Aero state and sizing. This works for GOCART and arbitrary generic modal configurations without embedding a CAM climate aerosol profile.
2. **Literal CAM-SIMA/BAM bulk emulation — rejected for this increment.** Copy BAM's constant `1/25` ice size weight, sulfate scaling, and zero cloudborne return. It improves BAM-specific numerical parity, but those constants are not general GOCART physics and CAM-SIMA bulk has unimplemented activation and transfer methods. A named compatibility profile can be considered separately after scientific review.

## Host-facing behavior

One common conceptual query applies to both representations, with population axis equal to configured bulk species or modal modes. The implementation may use separate bulk and modal entry points to preserve type safety, but the output fields and units should match. Accept an interstitial packed state and an optional *host-owned* cloudborne packed state of the same shape. For bulk, every slot is species mass and particle number is diagnosed from that population's mass. For generic modal, the input follows the existing per-mode species-mass and number layout, so cloudborne number must be supplied in its number slots. An omitted cloudborne state yields zero cloudborne mass and number with a distinct `Absent` status; it is never inferred from ambient state.

Return, per cell/level/population, total dry mass, number, dry number-median radius, wet radius, and configured/derived physical size bounds where available. Reuse `bulkFixedNumber`, `bulkReferenceRadii`, `bulkScavengingSizes`, `modalProperties`, and `modalScavengingSizes` math rather than creating a second interpretation. A zero-mass bulk population has zero diagnosed number; an empty modal population has zero number and its existing reference-radius status. Expose a separate populated/absent status so zero is not confused with missing host cloudborne input.

For a queried ice aerosol category, report separately: `type_weight` in [0,1], `size_weight` in [0,1], `eligible_number = number * type_weight * size_weight` in particles/m3, and a heterogeneous-freezing eligibility flag. Bulk type weight is 1 when the configured species carries that category and 0 otherwise. Modal type weight is the matching category's dry-mass fraction within the mode; a mode with no dry mass has weight zero. Each category is queried independently, including a distinct stratospheric-sulfate category, so overlap is visible to the caller rather than silently summed twice. The host decides how to combine category-specific candidates.

Ice category tags belong to physical species configuration, not to optical tables or name-string guesses. Add an explicit, validated cloud/ice role or category list to both structured and YAML configuration paths. Add explicit ice size rules scoped to bulk species or modal mode and queried category. The supported generic rules are `All`, `None`, `FixedFraction(f)` and `AboveDryDiameter(d)`; the last integrates the configured lognormal number distribution above a dry diameter cutoff. A rule may declare a separate weight for pre-existing ice; otherwise the same weight applies. No category or no rule means zero ice candidates and a status that says `NotConfigured`. Configured lower/upper bulk radii are metadata and do not silently truncate this lognormal integral. Add a separate explicit heterogeneous-freezing eligibility flag scoped to population and category, default false; neither a dust/black-carbon tag nor a nonzero ice size weight turns it on implicitly. The host can therefore opt into a documented ice scheme without Ex-Aero selecting a CAM-SIMA threshold, `1/25` constant, or heterogeneous-freezing rule by default. The optional CAM-SIMA-like rules (for example sulfate Aitken diameter above 0.1 µm) must be entered explicitly and reviewed by the host team.

For liquid cloud activation, report number and mass fractions at each requested supersaturation and temperature **for interstitial aerosol only**. Already cloudborne aerosol is excluded from this activation query to prevent double counting; its mass and number remain available through the separate cloudborne population query. Generic modal uses its present lognormal kappa-Koehler number fraction and the corresponding third-moment mass fraction (the CAM-SIMA modal `actfracs` form). GOCART bulk uses the current `computeCCN` all-or-none critical-supersaturation decision based on configured `dry_particle_diameter`, so its reported number and mass fractions are each 0 or 1. This preserves current bulk CCN sums; adopting continuous lognormal bulk activation would be a separate behavior change. A non-hygroscopic or empty population returns both fractions zero. Fractions are diagnostics only; the host computes activated numbers/mass and any transfer budget.

The host may use existing wet sizing, water, density, and dry-volume diagnostics to drive coating or water-activity schemes. This increment does not copy CAM-SIMA's coating-fraction or water-activity formulas, which are specific scientific parameterizations and require a separate review. It also does not implement CAM-SIMA bulk stubs such as `actfracs`, `hetfrz_species`, and `update_bin` literally; instead it supplies the useful information those call sites need under explicit Ex-Aero configuration.

## Interfaces and integration

- Add a standalone cloud/ice coupling header and implementation for each representation, following `BulkRemoval` and `ModalRemoval`. Where an initialized package is required, check its representation/coupling interface and fail clearly for the wrong backend.
- Expose stable category and rule codes and zero-based population indices through C ABI wrappers; mirror them in `exaero_interface.f90`. Accept caller-owned column-major arrays and a nullable cloudborne-state pointer. Document exact extents, units, and status codes in each language.
- Make configuration ingestion identical for YAML and structured callers. Existing configurations remain valid and produce unchanged diagnostics; no ice category is inferred from names such as `DU` or `SO4`.
- Keep `computeCCN` on its current path and test the sum of the new per-population number fractions times population number against it. Shared arithmetic is preferred where practical to prevent drift.
- Report configured cloud/ice rules and their source in the guide. Distinguish a generic mathematical size rule from a scientifically calibrated ice-nucleation parameterization.

## Verification and acceptance

- Bulk and modal queries produce the same fields and units, with finite values and weights/fractions in [0,1]. Category rules cover matching, nonmatching, empty, mixed-composition, cutoff-tail, and pre-existing-ice cases.
- Bulk fixed number matches the existing lognormal mass-to-number calculation and never consults optics or radiation tables. Modal number matches the packed prognostic field. Host cloudborne buffers are read without mutation; omission gives zero plus `Absent`.
- Bulk per-species and modal per-mode interstitial activated number sums reproduce existing `computeCCN` across multiple supersaturations, including kappa zero and empty populations. Cloudborne particles do not enter those sums. Modal mass fractions follow the lognormal third moment, and all activation fractions are monotone with supersaturation.
- C++, C, and Fortran paths agree on dimensions, ordering, statuses, and units. Invalid input fails without partial output or changes to state.
- Existing removal, sizing, chemistry, optics, and CCN tests continue to pass. If the full project toolchain is unavailable, record the missing dependencies and run the standalone public/API checks that can execute.

## Decisions left to the builder

Choose exact class/function names and how much computation to share across the existing Kokkos kernels and additive helper implementations. Preserve the specified numerical semantics and cross-language contract. Use an isolated feature branch/worktree when implementation begins.

## Reference points

- CAM-SIMA [bulk state](https://github.com/ESCOMP/CAM-SIMA/blob/development/src/aerosol/bulk_aerosol_state_mod.F90): fixed `1/25` ice size weight, species match weight, no cloudborne state update.
- CAM-SIMA [bulk properties](https://github.com/ESCOMP/CAM-SIMA/blob/development/src/aerosol/bulk_aerosol_properties_mod.F90): several activation/ice methods are stubs or no-ops.
- CAM-SIMA [modal state](https://github.com/ESCOMP/CAM-SIMA/blob/development/src/aerosol/modal_aerosol_state_mod.F90) and [modal properties](https://github.com/ESCOMP/CAM-SIMA/blob/development/src/aerosol/modal_aerosol_properties_mod.F90): lognormal ice size and activation fractions, type weighting, and transfer policies.
- Ex-Aero `BulkRemoval.hpp`, `ModalRemoval.hpp`, `ModalCoupling.hpp`, `GocartConfig.hpp`, and `ModalConfig.hpp` define the existing additive API and physical metadata surfaces.
