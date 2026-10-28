# Shortwave and Longwave Aerosol Optics Design

## Intent

Ex-Aero should serve both radiation schemes that request a named band set and
consumers that request optical properties at physical wavelengths. Both bulk
GOCART and generic modal representations should support shortwave and
longwave radiation where validated data are available. The two coordinate
systems must remain explicit: an RRTMG band index is not a wavelength.

## Why

The public `computeOptics` contract takes wavelengths in metres and returns
extinction, scattering, asymmetry, backscatter, lidar backscatter, and column
optical depths. The baked GEOSmie bulk tables instead store 30 RRTMG band
averages indexed 1–30 (14 SW and 16 LW). The bulk optical kernel currently
passes physical wavelengths to this band-indexed axis, so normal radiation
wavelength requests do not use those native tables.

Generic modal optics has a physical-wavelength provider and generated starter
tables at 355, 550, 870, 1064, and 1640 nm. Those tables cover visible/NIR SW,
not the full SW range or thermal IR. The guides explicitly describe that
limit. The native GEOSmie monochromatic files cover 250 nm–40 µm but are not
part of the default baked snapshot. CAM-SIMA's aerosol optics contract
reinforces the radiation split: SW consumers use extinction/absorption/SSA/g,
while LW consumers request absorption.

## Goals

- Preserve existing arbitrary physical-wavelength `computeOptics` behavior
  and ABI.
- Add a separately named band-indexed operation for a named band set, starting
  with the source-backed RRTMG 30-band set: 14 shortwave followed by 16
  longwave bands.
- Make band metadata, ordering, and SW/LW partition queryable and stable at
  the C++/C/Fortran boundary.
- Use native GEOSmie band averages for the bulk species that have them.
- Extend generic modal optics to provide both wavelength and SW/LW band
  products using documented material and spectral assumptions.
- Keep backend differences explicit: a missing band product reports a
  structured unsupported/error outcome; it never falls back to a band
  midpoint, ADT, or unrelated coordinate.
- Preserve output units, aerosol state layout, and the current optics output
  slots where their meanings remain valid.
- Record data provenance and spectral/band domains for each product.

## Non-goals

- Replacing or changing the host radiation transfer solver.
- Adding cloud or gas optics, radiative heating tendencies, or aerosol state
  mutation.
- Treating the 30 RRTMG band IDs as physical wavelengths.
- Claiming all spectral points or materials are supported when source optical
  constants do not cover them.
- Adding a new optical treatment for brown carbon until its data source is
  validated.
- Changing the existing aerosol microphysics, hygroscopic growth, or modal
  mixing assumptions beyond what is required to evaluate the new optical
  products.

## Approaches considered

1. **Convert everything to a continuous spectrum and derive radiation bands.**
   This gives one conceptual table, but it needs validated spectral weighting
   and optical-constant coverage for every band. It cannot reproduce the
   existing GEOSmie native RRTMG averages by treating each band as a point.
2. **Maintain two explicit products and interfaces.** Keep physical-wavelength
   queries and add a band-set-indexed API. Use native GEOSmie RRTMG averages
   where supplied; generate modal band products using a documented,
   validated RRTMG band-averaging recipe. Preserve domain/provenance and report
   unsupported products explicitly. **Selected.** This matches both host
   workflows and preserves the intended meaning of the existing data.
3. **Use only radiation bands.** This simplifies radiation coupling but drops
   general wavelength diagnostics, photolysis use, and existing provider
   behavior. It does not meet the request for both products.

## Design

### Separate coordinate spaces

`computeOptics` remains the physical-wavelength operation. Add an additive
band operation that takes a `RadiationBandSet` identifier and returns the
complete band axis in canonical order. The first supported identifier is
`RRTMG_30`; its entries are zero-based at the API boundary and map explicitly
to source IDs 1–30. Metadata reports 14 SW bands, 16 LW bands, the source
convention/version, and the band edges/centers where available. No implicit
conversion between wavelength and band index is allowed.

Both operations return the existing volume-property slots and column
integrals. The band operation additionally identifies the band set so a
buffer cannot be interpreted without its coordinate metadata. SW callers may
consume extinction, scattering/absorption, SSA, and asymmetry. LW callers may
consume absorption (derived consistently from extinction minus scattering
when source fields require it). Ex-Aero does not perform gas optics or
radiative transfer.

### Bulk backend

For the RRTMG band operation, consume the pinned GEOSmie band-averaged tables
directly for DU, SS, SU, BC, OC, and NI. Respect their radius and RH
coordinates and documented interpolation rules. BR remains unsupported unless
an explicit source supplies matching RRTMG band data. The existing physical
wavelength operation uses only curves that declare physical wavelength
coordinates; when a request is outside a curve's domain, keep the established
explicit fallback/error semantics and never reinterpret band IDs as metres.

The bulk backend's treatment of species without native RRTMG tables must be
explicit in per-band status output. The legacy RH-only lookup and ADT path
remain available to the physical-wavelength operation; they are not an
automatic substitute for missing RRTMG data.

### Generic modal backend

Retain the `IModalOpticsProvider` wavelength interface. Extend the generated
modal optical data workflow to cover an approved physical wavelength domain
across SW and LW, using pinned refractive-index inputs and provenance. Do not
extrapolate beyond the material source domain.

Generate `RRTMG_30` modal band products with a documented band-averaging
procedure and weighting data. The generator must validate that every band's
spectral support is covered for every material advertised as supported. A
band may be marked unsupported if its required source/weighting range is
missing; it must not be populated from its center wavelength. Preserve current
component and internal-mixture semantics, dry-mass normalization, and
single-application humidification. Native MAM4xx uses the same band metadata
contract; its availability is based on its own material/source coverage, not
assumed from the generic modal tables.

The exact spectral domain, wavelength grid/quadrature, and RRTMG band-weight
source for the regenerated modal snapshots are builder decisions only after
verifying the pinned optical-constant data and upstream band-averaging
procedure. These choices must be recorded in the snapshot manifest and
provenance. If a band cannot be justified from those inputs, stop and report
that gap rather than introduce a surrogate value.

### Public API and host contract

- Add Kokkos-free band-set and band metadata definitions in the public C++
  interface.
- Add an additive virtual/interface path that does not change the existing
  `IAerosolPackage` vtable; use the established separate extension-interface
  pattern for modal coupling where appropriate.
- Add C and Fortran entry points and constants with explicit band-set ID,
  band count, SW/LW counts, output layout, statuses, and error handling.
- Keep the existing `computeOptics` signature and output layout unchanged.
- Update host integration guidance to show physical wavelengths and radiation
  bands as separate calls, including buffer extents and status handling.

### Validation

- Verify all 30 RRTMG metadata entries map to source IDs 1–30, with a 14/16
  SW/LW split and monotone/valid band boundaries.
- Compare bulk band results against direct pinned-table values at exact RH and
  radius nodes, and test interpolation between nodes.
- Test a backend/species without native band data and ensure it reports
  unsupported without using ADT, a band center, or zero as a silent result.
- Test physical wavelength queries continue to use physical-coordinate data
  and reject/flag band-index coordinates as incompatible.
- Test modal generated spectral coverage and band integration against
  independent generator reference values for representative absorbing,
  scattering, hygroscopic, and coarse particles.
- Check mass-normalized extinction/scattering/absorption relationships,
  single-scattering albedo, asymmetry bounds, finite outputs, RH behavior,
  and per-band column AOT.
- Exercise C++ and C/Fortran APIs for both paths and verify dimensions,
  ordering, band metadata, unsupported statuses, and error propagation.
- Run deterministic snapshot generation checks and existing bulk/modal optics
  tests.

## Constraints and assumptions

- Source-backed initial radiation band set is RRTMG-30 (14 SW + 16 LW), from
  the pinned GEOSmie table snapshot. Other band sets require their own IDs,
  metadata, and data; they cannot reuse RRTMG indices implicitly.
- Bulk GEOSmie tables are band averages, not monochromatic samples.
- Existing generic modal wavelength nodes are starter data, not enough for
  complete SW/LW spectral coverage.
- The RRTMG band averaging weights and optical-constant coverage must be
  verified before claiming modal coverage. No spectral extrapolation or
  band-center approximation is accepted.
- CAM-SIMA reference: `aerosol_optics_mod.F90` defines distinct SW and LW
  property calls; `aerosol_optics_core.F90` passes SW band count, LW band
  count, and the separate SW/LW refractive-index arrays.

## Builder decisions

- Determine the canonical RRTMG band metadata source and preserve its exact
  band ordering and edges.
- Validate which pinned refractive-index inputs and spectral tables support
  all RRTMG bands, particularly the longest-wave LW interval.
- Specify whether modal band coefficients are precomputed at generation time
  or integrated from a generated physical-wavelength grid at runtime, choosing
  the method that is reproducible and keeps the runtime path performant.
- Use explicit availability status for partial band sets and preserve atomic
  output behavior on invalid requests.
