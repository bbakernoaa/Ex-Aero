# Generic modal aerosols and MKPP coupling

Explicit removal-coupling utilities and dry size/number bounds are documented
in [modal removal and bounds](modal-removal.md).

Select `createPackage("modal")` or `exaero_create_package("modal", ...)` and
initialize with YAML or `initializeModal(package, ModalConfig)`. The backend uses
the existing Kokkos and yaml-cpp dependencies; it does not require MAM4xx or EKAT.
It is a generic, experimental two-moment model inspired by modal aerosol methods,
with explicit process equations and limits below. It does not reproduce MAM4 or
MADE as a complete scientific model.

Editable MAM4 mode/material presets and optional native coagulation/water
processes are documented in [the MAM4 profile guide](modal-mam4.md). The generic
profile remains independent of MAM4xx; native process choices require it.

## State, sizing, and configuration

Each mode packs all configured species masses in species order, followed by its
particle number. Query `stateFields()` before allocating arrays. Names are
`<species>_<mode>_i` and `num_<mode>_i`. Mass is kg/m3 of air and number is
particles/m3, including on the C/Fortran boundaries. Arrays are host memory,
`layout_left`, `(cell,level,field)`. Kernels stage to the selected Kokkos space and
return synchronously. Release packages before finalizing the environment.

With the default generic process choices, every mode specifies a reference **dry geometric radius** in metres and a fixed
geometric standard deviation greater than one. Its evolving radius is

\[
 V=\sum_s M_s/\rho_s,\qquad
 r_g=\left[\frac{3V}{4\pi N\exp(4.5\ln^2\sigma_g)}\right]^{1/3}.
\]

The reference radius is used for emitted-number inference and the empty-mode
property convention; it never overrides a populated mode's mass/number sizing.
An empty mode has zero mass and number. The caller must remove residual number
when MKPP evaporates a mode completely; EX-aero rejects ghost particles.
Molecular weights are **kg/mol**, unlike some legacy GOCART settings.

`species` within a mode lists permitted materials; omission permits all. All
species retain packed slots, including disallowed zero slots. Coagulation and
transfer destinations must permit every material carried by their sources.
This allows a coagulated particle to retain its coatings.

See [three-mode](../../examples/modal_three_modes.yaml) and
[four-mode](../../examples/modal_four_modes.yaml) examples, including dust and
ammonium nitrate. The sample organic/dust molecular weights are host mechanism
surrogates, not material constants. The sample nucleation coefficient is an
illustration, not an atmospheric calibration.

## Processes implemented in EX-aero

- **Coagulation:** Fuchs flux-matching Brownian collision kernel, with
  Cunningham-corrected particle diffusion and Sutherland viscosity, integrated
  with eight-point Gauss-Hermite quadrature in each participant's log radius.
  The number and material-volume-weighted kernels are separate. Intermode
  destinations are required for every unordered pair; self collisions retain
  material and reduce number. Explicit substeps limit fractional net depletion
  to 10%; exceeding 100000 substeps or numerical underflow fails atomically.
- **Growth transfer:** one outgoing route per source, with an explicit dry-radius
  threshold. When a mode's geometric radius exceeds that threshold, transfer its
  number tail and third-moment material tail above the threshold, then reclose
  at each mode's configured width. Apply routes in topological order. This is a
  declared modal reclosure approximation; repeated calls are not a sectional
  solution of the size-distribution equation.
- **Default water:** volume-weighted kappa and equilibrium growth
  \(r_w/r_d=[1+\kappa RH/(1-RH)]^{1/3}\), with RH capped at 0.98 for sizing.
  This neglects the Kelvin term in subsaturated sizing and does not implement
  deliquescence/efflorescence hysteresis. Water mass and wet density use the
  same volume budget.
- **Diagnostics/CCN:** dry mass, independent number, wet surface area, water,
  dry-geometric PM1, wet-aerodynamic PM2.5/PM10, Stokes settling, surface/column
  integrals, and the lognormal kappa-Kohler CCN spectrum. Settling/cuts assume
  spheres and omit slip/shape corrections. CCN is a diagnostic, not an activation
  tendency. Level zero is the surface.
- **Emissions:** already-speciated mass/number rates, or area fluxes divided by
  positive layer thickness. `emittedNumber(mode, masses)` supplies a linear
  conversion from consistent species masses/rates/fluxes to particle number.

Nucleation, condensation/evaporation, and chemistry are integrated by MKPP.
Cloudborne evolution, aging, cloud activation tendencies, and deposition process
rates are not implemented in this backend. Explicit deposited-budget rebinning
is available through the removal coupling utilities linked above.

## Optical tables and mixing treatments

The selected treatment is per mode:

- `component_tables` (default): sum constituent extinction/scattering and
  scattering-weight the asymmetry. This is an **external-mixture optical
  approximation**, even though the mass/number state is internally mixed.
- `internal_mixture`: install a C++ `IModalOpticsProvider`. It may use a
  composition-aware table or a Mie solver, must validate its domain, and returns
  complete mode integrals. Without a provider the operation fails explicitly.
  Built-in component tables are not substituted for this treatment.

The new baked sources use the Mie engine and refractive-index inputs from
`bbakernoaa/GEOSmie` revision
`c535388356c7e8982ead98deeebff33985254f9c`. They are particle-resolved,
dry-mass-normalized, and independent of modal spread:

| Source | Material approximation | Density kg/m3 | Kappa |
| --- | --- | ---: | ---: |
| `MODAL_SU` | OPAC sulfate/H2SO4-equivalent optics | 1700 | 0.507 |
| `MODAL_BC` | OPAC soot | 1000 | 0 |
| `MODAL_OC` | OPAC water-soluble organic surrogate | 1800 | 0.14 |
| `MODAL_SS` | OPAC sea salt | 2200 | 1.16 |
| `MODAL_DU` | Spherical mineral dust, `ri-dust.wsv` | 2650 | 0 |
| `MODAL_NI` | Ammonium-nitrate-equivalent material, `ri-nitrate.wsv` | 1725 | 0.67 |

The nitrate normalization is the **complete dry salt**, not nitrate-ion mass.
MKPP mappings must respect that composition basis. Dust is spherical Mie dust;
these tables do not model dust nonsphericity. The chosen kappas are explicit
recipe assumptions, rather than the original bulk tables' growth curves.

Bind with `optics: {source: MODAL_SU}`. Authoritative generated metadata supplies
the table semantics and physical wavelength axis; material density and kappa
must match the recipe. `modalCoupling(package).opticalTables()` exposes source
versions, labels, semantics, and citations. Packages own immutable table snapshots.

Coverage: dry particle radii 0.1 nm–50 um; RH nodes 0, 0.5, 0.8, 0.9, 0.95, 0.98;
wavelengths 355, 550, 870, 1064, and 1640 nm. RH interpolation is linear and
spectral interpolation is linear in log wavelength. These starter visible/NIR
curves are **not RRTMG band averages** or a full thermal-infrared dataset.

Native RRTMG calculations use the separate C++ `IOpticsBandCoupling` surface;
they do not pass band IDs through `computeOptics`. Query `bandAvailability()`
for each mode/material before requesting bands. The committed generic-modal and
MAM4xx snapshots are still wavelength-only and therefore report
`NotInSource` for RRTMG bands. A component table can supply a separate
`@axis lambda band_index` coordinate with native `bext`, `bsca`, `bbck`, and `g`
values, or an internally mixed mode can register an
`IModalBandOpticsProvider`. The package never substitutes wavelength nodes or
band centers for missing source data, and positive-mass unsupported requests
leave the caller's output unchanged. Modal band tables are not yet included in
the checked-in generated snapshots.
The smallest tabulated radii extend interpolation support and should not be
interpreted as validation of continuum sphere physics at molecular scales.

Particle-resolved coefficients are integrated analytically against lognormal
mass over every linear radius interpolation interval. At least `1-1e-8` of dry
mass must lie inside the table radius domain; otherwise the query fails.
No optical size extrapolation or silent treatment fallback occurs. Fine/PM2.5
integrals use the same interval integration; aerodynamic cuts use the modal
wet-density approximation. Hemispheric backscatter uses a Henyey-Greenstein
closure, while lidar backscatter uses the source's 180-degree coefficient.

The existing loader, source bindings, overrides, and provenance are reused.
Legacy bulk curves require an explicit `kind` (`particle_resolved` or
`lognormal_average`) and physical `wavelength_nodes`. A lognormal average also
requires a matching `geometric_std_dev`; it is not integrated a second time.
The old RRTMG coordinates are band IDs and must never be passed as metres.
Duplicate size coordinates representing different populations cannot be treated
as a continuous size curve. Humidification is applied once in the optical recipe;
EX-aero does not multiply the coefficients by an additional growth factor.

Balanced condensation/evaporation RHS terms, analytic coupled Jacobians,
coating-driven aging and post-acceptance renaming are documented in the
[chemical coupling guide](modal-chemistry.md). These are additive C++/C/Fortran
interfaces; MKPP remains the chemical integrator.

## MKPP contract and validation

Include `ModalCoupling.hpp` and obtain the separate `IModalCoupling` interface
through `modalCoupling(package)`. This leaves the `IAerosolPackage` vtable and
existing C entry points unchanged. New C functions are in `ModalCoupling_C.h`
and mirrored in `exaero_interface.f90`.

`modalProperties` returns `(cell,level,mode,9)` with dry/wet radius, number,
mass, dry volume, kappa, wet surface area, wet density, and water. Optional
derivatives are `(cell,level,mode,property*field_count+field)`, holding RH fixed.
Native MAM4xx water rejects wet-property Jacobian requests; dry uptake
derivatives remain available. Empty-mode radius derivatives are defined as zero; new-particle seeding is a
separate source operation. `modalUptake` returns dry-geometry Fuchs-Sutugin
transport coefficients [s^-1], with optional analytic state derivatives.
MKPP supplies diffusivity [m2/s], molecular weight [kg/mol], accommodation,
equilibrium concentrations, reaction laws, and the implicit integration.

`ModalMechanismMap.hpp` resolves field names to mechanism indices and explicit
conversion factors. Mass kg/m3 converts to molecules/cm3 using
`6.02214076e17 / molecular_weight`; particle number converts using `1e-6`.
Apply inverse factors to returned states/tendencies and the corresponding
chain-rule factors to Jacobians. Every aerosol field must map exactly once.
Other MKPP concentration conventions require an explicit host conversion.

`nucleationTerms(H2SO4_kg_per_m3)` evaluates `J=A*C^p`, where C is H2SO4
molecules/m3. It returns number production and balanced seed-mass/precursor-loss
rates without modifying gas or aerosol. The empirical helper supports pure dry
H2SO4-equivalent seeds; mixed-cluster chemistry belongs to MKPP. Seed mass uses
the receiving mode's spread and configured seed geometric radius. MKPP must
limit depletion through its coupled solver and include nucleation number as a
separate prognostic variable.

The host advances MKPP chemistry, advances EX-aero coagulation/transfer, then
refreshes properties. Geometry must be reevaluated at MKPP stage states if its
Jacobian is to include geometry feedback. Freezing geometry is a host numerical
approximation; the supplied derivatives do not by themselves implement an MKPP
Jacobian. No MKPP repository changes or unified coagulation/chemistry solver
are included here.

Build with CMake/C++23, Kokkos, yaml-cpp, Fortran, and Python NumPy/netCDF4/xarray.
Tests additionally use Catch2 v2; the Python generator tests use pytest and Dask.
Run `ctest --test-dir build --output-on-failure` and
`python -m pytest tests/test_modal_optics_generator.py`.
Checks cover dynamic layouts, analytical moments/derivatives, independent
Brownian integration, conservation, timestep refinement, optical moment
integration, Mie/Rayleigh limits, unit budgets, C/Fortran boundaries, and exact
snapshot-to-header serialization. These are limited software/science checks,
not production atmospheric validation. GPU execution needs separate validation.

References: Fuchs (1964), *The Mechanics of Aerosols*; Seinfeld & Pandis (2016),
*Atmospheric Chemistry and Physics*; Petters & Kreidenweis (2007),
doi:10.5194/acp-7-1961-2007; Hess et al. (1998), OPAC,
doi:10.1175/1520-0477(1998)079<0831:OPAACT>2.0.CO;2; and the pinned
[GEOSmie documentation](https://github.com/bbakernoaa/GEOSmie/tree/c535388356c7e8982ead98deeebff33985254f9c).

The optional [MAM4xx adapter](aerosol-backends.md#mam4xx-component-optics) uses
this same component integrator with a separate GEOSmie recipe matching native
MAM4xx densities and kappa. The `MAM4_*` curves are also available to structured
generic-modal callers that use those exact material properties.
