# Representation-aware aerosol removal design

This is the removal contract within the broader
[common host-coupling design](2026-10-29-host-coupling-design.md). The expanded
plan adds metadata and capability queries first, then delivers common cloud/ice,
property and deposition utilities in subsequent slices.

## Intent

Give host aerosol schemes one Ex-Aero tendency API that can express different
mass removal rates for bulk species/bins and for individual modal constituents.
Modal particle-number removal remains a separate mode-level tendency because
number describes particles, not individual material constituents.

## Current behavior

Bulk removal accepts one mass-scavenging coefficient per species/bin and derives
the diagnostic number tendency from that bin's fixed lognormal distribution.
Generic modal removal accepts one mass coefficient per mode, applies it to each
constituent in that mode, and accepts an independent number coefficient per
mode. This cannot express tracer-dependent uptake within a mode.

The initialized package already reports whether its representation is bulk or
modal. The API can therefore dispatch from package metadata without a second
caller-supplied representation selector.

## Goals

- Provide a common C++ entry point, with matching C and Fortran bindings, for
  host-supplied first-order removal coefficients.
- For bulk packages, accept mass coefficients and return mass and diagnostic
  number tendencies per species/bin.
- For generic modal packages, accept a mass coefficient per constituent slot
  and a separate number coefficient per mode; return constituent mass
  tendencies and mode-level number tendencies.
- Preserve atomic output behavior, validation, column-major host-memory
  layouts, and no state mutation.
- Retain the current mode-shared modal helper as a convenience for whole-particle
  removal, including the existing velocity adapter.

## Non-goals

- Implement collection efficiencies, precipitation processing, activation,
  cloudborne transfer, resuspension, time integration, or accepted-state
  updates.
- Infer modal number removal from species-specific mass loss. The host remains
  responsible for providing a physically consistent mode-level number rate.
- Change the native MAM4xx state layout or its APIs.

## Approaches considered

1. Add a sibling modal species-rate helper. This preserves the current split
   between bulk and modal entry points, but duplicates the public capability
   and leaves hosts to select functions themselves.
2. Replace the current modal helper with per-constituent mass rates. This
   simplifies the modal surface but breaks existing callers that provide one
   whole-particle mass rate per mode.
3. Add one representation-aware rate entry point, dispatched from the package,
   and retain the current mode-shared modal helper as a convenience. This gives
   hosts one general path, avoids a redundant selector argument, and preserves
   whole-particle callers. **Recommended.**

## Proposed API contract

The common rate function accepts the initialized package, packed state, mass
rate array, an optional number-rate view, mass tendency output, and number
tendency output. Package representation determines the expected shapes and
whether the number-rate view is present:

- Bulk: state, mass rates, and mass tendencies use `(cell, level, species/bin)`.
  Number rates are not supplied; diagnostic number tendency follows the same
  fractional loss as the corresponding species/bin mass.
- Generic modal: state, mass rates, and mass tendencies use the packed
  `(cell, level, field)` layout. Rates in number-field slots must be zero. A
  number-rate view is required and uses `(cell, level, mode)`; number tendency
  has the same extent. Each constituent mass tendency uses that constituent's
  own coefficient.

All coefficients are finite, nonnegative inverse seconds. Outputs remain
disjoint from inputs and each other. Failure leaves outputs unchanged. The
existing representation-specific functions remain available; the common
function is an additional host-facing path.

## Implementation details resolved in the plan

- The public C++ entry point is `removalRateTendencies`; it derives bulk or modal
  behavior from the initialized package and takes a nullable number-rate view
  (null for bulk, required for modal).
- The C entry point is `exaero_removal_rate_tendencies`, with state-field and
  population extents supplied by the caller; its number-rate pointer is null for bulk
  and required for modal. Fortran passes `c_null_ptr` or a contiguous array's
  `c_loc` as appropriate.
- The common velocity entry point uses one velocity per population and the same
  metadata-derived output shapes. Population count equals species/bin count for
  bulk and mode count for modal; it is never zero for an initialized package.
- Existing bulk/modal APIs remain behaviorally unchanged. The common function
  delegates to those tested paths, adding a modal constituent-rate conversion.
- Tests cover distinct constituent mass losses in one mode, independent mode
  number loss, bulk diagnostic number behavior, C/Fortran calls, invalid shapes,
  zero number-field rates, aliasing, overflow, and atomic failure.
