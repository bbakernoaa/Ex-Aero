# Ex-Aero aerosol parity with CAM-SIMA — design

## Intent

Bring Ex-Aero's GOCART bulk and generic modal aerosol surfaces close to the
working diagnostic and modal-physics behavior in CAM-SIMA, while retaining the
Ex-Aero/MKPP split: MKPP/host code owns chemical formation, equilibrium
thermodynamics, and integration; Ex-Aero supplies aerosol geometry, diagnostics,
and tendency adapters for host-computed rates. The comparison target is CAM-SIMA `development` at
`52723a73a2a49e37fdd981f795d08084bb5bfda8`. The Ex-Aero starting point is
`codex/native-aerosol-optics` at `5b347cf` plus this design.
CAM-SIMA's pinned `src/physics/ncar_ccpp` submodule is
`ESCOMP/atmospheric_physics` at `6ee9b89afb294408529bea5b72c7bab47c0c55eb`;
its aerosol wet-deposition core is part of the comparison scope.

## Goals

- Map every in-scope CAM-SIMA aerosol capability to an Ex-Aero equivalent,
  an intentional representation difference, or a documented remaining gap.
- Match the numerical behavior of implemented CAM-SIMA modal formulas where
  the state definitions and physical assumptions are comparable.
- Add reference cases and cross-interface tests before extending behavior.
- Preserve Ex-Aero's bulk/modal representations and host-owned step ordering.
- Let host schemes compute dry-deposition velocities and wet-scavenging rate
  coefficients and apply accepted updates; Ex-Aero may translate supplied
  species/bin-specific bulk rates and modal mass/number-moment rates into tendencies
  without mutating aerosol state.
- Complete useful shortwave and longwave band optics for modal data wherever
  the pinned source refractive-index data provide valid coverage.
- Keep C++, C, and Fortran APIs consistent without changing the existing
  `IAerosolPackage` virtual interface.

## Non-goals

- Move chemical reaction formation, nucleation parameterizations, equilibrium
  thermodynamics, or the implicit chemical integrator from MKPP/host code into
  Ex-Aero.
- Modify MKPP or automatically call ISORROPIA-Lite from an MKPP-generated
  mechanism in this workstream.
- Port CAM-SIMA stubs as if they were working features. Its bulk `actfracs`,
  `scav_diam`, `resuspension_resize`, and `rebin_bulk_fluxes` contain explicit
  unimplemented behavior and do not define a parity target.
- Turn Ex-Aero's diagnostic cloud/ice helpers into an autonomous cloud
  microphysics scheme. The host continues to own cloudborne state and accepted
  transfer tendencies.
- Calculate dry-deposition velocity, wet-scavenging or collection rates, choose
  active removal levels, apply a removal tendency, or advance aerosol state in
  an Ex-Aero removal scheme.
- Claim full CAM-SIMA model parity, production readiness, or GPU parity from
  unit tests alone.

## Constraints and current state

- CAM-SIMA modal behavior is the numerical reference only for formulas and
  semantics it actually implements. Reference revision and all intentional
  deviations must be recorded alongside tests.
- Ex-Aero's bulk GOCART number is diagnosed from its fixed configured
  lognormal distribution; bulk number is not a prognostic field. Modal number
  remains an independent prognostic field.
- Existing Ex-Aero modal helpers include number/radius bounds, activation
  fractions, scavenging sizes, removal-budget redistribution/rebinning,
  cloud/ice diagnostics, and an RRTMG band interface. These require a behavior
  audit against CAM-SIMA, not an assumption of parity from matching names.
- In the pinned CAM-SIMA aerosol-property API, modal `scav_diam(bin_ndx)`
  returns `dgnum(bin_ndx)`, one representative diameter per modal bin. The
  property API itself does not apply wet removal. The actual wet-deposition
  core is in the pinned physics submodule,
  `schemes/aerosol/core/wetdep.F90`; inspecting `src/aerosol` alone does not
  establish the full CAM-SIMA wet-deposition behavior.
- That core's
  [get_bcscavcoefs](https://github.com/ESCOMP/atmospheric_physics/blob/6ee9b89afb294408529bea5b72c7bab47c0c55eb/schemes/aerosol/core/wetdep.F90#L805-L862)
  returns separate number and volume/mass below-cloud coefficients for mode m,
  with a shared size-dependent collection kernel integrated against the
  different distribution moments. Its `wetdepa_v2` distinguishes interstitial
  and cloudborne aerosol, in-cloud/below-cloud and convective/stratiform terms,
  and resuspension. These are host-scheme responsibilities for Ex-Aero.
  `wetdepa_v2` takes a precipitation-normalized `scavcoef` in `/mm` and a
  timestep; Ex-Aero's proposed local-rate adapter takes effective first-order
  coefficients in `/s`. The host must supply the appropriate effective rates;
  raw CAM-SIMA scavenging coefficients are not interchangeable with them.
- CAM-SIMA groups constituent masses and one number field within each mode
  ([modal state accessors](https://github.com/ESCOMP/CAM-SIMA/blob/52723a73a2a49e37fdd981f795d08084bb5bfda8/src/aerosol/modal_aerosol_state_mod.F90#L163-L229))
  and provides separate number and mass activation fractions
  ([actfracs](https://github.com/ESCOMP/CAM-SIMA/blob/52723a73a2a49e37fdd981f795d08084bb5bfda8/src/aerosol/modal_aerosol_properties_mod.F90#L497-L516)).
  These support the mode/moment distinction; activation fractions are not
  themselves wet-deposition coefficients.
- Mode-specific removal is species-dependent. The pinned physics core
  [aero_convproc](https://github.com/ESCOMP/atmospheric_physics/blob/6ee9b89afb294408529bea5b72c7bab47c0c55eb/schemes/aerosol/core/aero_convproc.F90#L1445-L1479)
  forms volume-weighted mode hygroscopicity from constituent density,
  hygroscopicity, and mass, including a configurable POM hygroscopicity
  treatment explicitly intended to affect primary BC/POM wet scavenging.
  CAM-SIMA's modal `soluble()` also distinguishes primary-carbon modes.
  A shared mass coefficient inside a uniformly mixed mode does not mean
  different modes or total chemical species have the same removal rate.
  Species occupy different modes and reservoirs; their total losses reflect
  those populations' composition-dependent removal and their own abundances.
  Ex-Aero must supply current composition and effective size/hygroscopicity/
  density information so the host can calculate these differing coefficients.
  No composition-independent or globally fixed mode coefficient is implied.
- The agreed simple dry-deposition adapter accepts one finite nonnegative
  velocity `[m/s]` per cell, level, and bulk species/bin or modal mode, returning
  `-velocity/layer_thickness * concentration`. It requires finite positive
  layer thickness and applies that population's velocity to its mass and
  number. This convenience path does not resolve different distribution moments.
- The proposed first-order wet-removal adapter preserves host-supplied
  coefficients `[1/s]` per configured bulk species/bin. Modal accepts one mass
  coefficient and an independent number coefficient per mode, cell, and level.
  Bulk tendencies are `dM_s/dt = -lambda_s*M_s`; modal tendencies are
  `dM_ms/dt = -lambda_mass_m*M_ms` for every constituent s in mode m, and
  `dN_m/dt = -lambda_number_m*N_m`. Distinct species/bins and modes may have
  distinct rates. No global wet-removal coefficient is imposed.
- Wet deposition removes particles. With uniform composition within an
  internally mixed mode, all material components share the same fractional
  mass-removal coefficient, which the adapter applies to all components of
  that mode. Size-selective scavenging can give different mass and number
  coefficients; the host supplies effective coefficients derived consistently
  from its scavenging closure. Arbitrary chemical-component coefficients
  within one uniformly mixed mode are outside this adapter's contract.
  Distinct particle populations should be represented distinctly.
- Host solubility and activation assumptions are distinct from ice-candidate
  categories. The host wet scheme owns those factors; Ex-Aero must not infer
  sulfate/BC/POM scavenging rates from names or reuse ice eligibility as a
  wet-scavenging efficiency. Crosswalk the existing public property/configuration
  surfaces to these needs before adding any missing metadata utility.
- Both adapters validate before publishing outputs and leave input state
  unchanged. The host sets an individual input rate to zero where it wants no
  sink and applies returned tendencies only after accepting its step. Ambient
  and cloudborne reservoirs use separate snapshots and rate inputs; a rate
  must never be applied to both implicitly.
- Ex-Aero's existing dry/wet size and cloud/ice activation diagnostics can help
  a host compute wet-scavenging coefficients. Precipitation, cloudborne state,
  in-cloud/below-cloud collection rates, and the resulting coefficients remain
  host-owned; the tendency adapter does not infer wet-deposition rates.
  A local first-order sink is only one part of wet deposition: precipitation
  transport, evaporation/resuspension, and surface-budget accounting remain
  part of the host scheme. The adapter is not a complete wet-deposition model.
- Bulk number remains diagnostic: each species/bin's number tendency follows
  its own mass-removal coefficient and configured fixed lognormal conversion;
  no number field is added to bulk state. Modal number tendency uses the
  independently supplied number coefficient and prognostic number field.
  Material mass tendencies are returned per existing mass slot.
- Existing resuspension redistribution and deposition rebinning helpers only
  transform host-supplied returned rates or deposited budgets; they do not
  calculate a removal rate, flux, or state update.
- GOCART has bundled native RRTMG band data. Generic modal and MAM4xx snapshots
  remain wavelength-only at the starting commit, so their native band
  capability is not available from the bundled data yet.
- Modal optical source coverage must remain explicit and must not be filled by
  wavelength substitution, band-center substitution, or extrapolation.
- MKPP keeps ownership of condensable concentrations, equilibrium
  concentrations, reaction laws, positivity, accepted/rejected steps, and
  integration. Ex-Aero returns terms and Jacobians only.

## Approaches considered

1. **Direct port of CAM-SIMA methods.** This would reproduce method structure
   and formulas quickly, but risks importing CAM-specific state assumptions
   into Ex-Aero's generic package and host-owned interfaces.
2. **Parity-tested slices — selected.** Pin CAM-SIMA as a numerical reference,
   write small reference cases, then implement or adjust one capability group
   at a time through Ex-Aero's existing C++/C/Fortran surfaces. This preserves
   the requested MKPP boundary and makes differences visible and reviewable.
3. **Optics-first.** Completing modal RRTMG data would close the most visible
   current gap, but would postpone missing modal diagnostic and process
   behavior and would not establish a parity definition for those features.

The selected approach is the parity-tested slice sequence below. The user
specified that CAM-SIMA's implemented formulas are the numerical reference;
the implementer may choose specific code structure, but must document and test
any numerical or state-transition deviation.

## Proposed parity slices

### 1. Capability matrix and reference harness

- Inventory CAM-SIMA's aerosol diagnostics, modal numerical methods, cloud/ice
  hooks, and SW/LW optical outputs against Ex-Aero's current public APIs.
- For each row, record source path and routine, state assumptions, Ex-Aero
  entry point, status (`equivalent`, `partial`, `intentional difference`,
  `missing`, or `not applicable`), and the evidence needed to close it.
- Capture deterministic reference cases from the pinned CAM-SIMA formulas for
  the in-scope methods. Keep a path to running a direct CAM-SIMA comparison
  when its full build environment is available.
- Specify per-method tolerances and edge cases before changing the method.

### 2. Modal diagnostic and numerical parity

Audit and close the applicable behavior in these groups:

- Lognormal geometry and lower/upper number limits, including empty modes and
  preservation of material mass.
- Number and mass activation fractions and dry/wet scavenging size semantics.
- Ice-nucleation eligibility and species/mode classification, without adding
  an Ex-Aero nucleation rate parameterization that belongs to MKPP.
- Resuspension redistribution and deposited-flux rebinning, preserving material
  and number budgets and making the host's rate/budget units explicit.
- Modal cloud/ice diagnostics using host-owned cloudborne state, including
  activation and candidate fractions where CAM-SIMA provides a comparable
  defined formula.
- Host-removal tendency adapters: accept either an externally computed
  population dry-deposition velocity per cell and level, or first-order
  wet-scavenging coefficients per bulk species/bin or modal mass/number moment,
  then return material-mass and number tendencies. The velocity path scales
  by layer thickness; the rate path uses each supplied `[1/s]` coefficient
  directly. These adapters do not calculate rates or update state. Bulk number
  tendency is diagnostic; modal number uses its own supplied coefficient and
  prognostic slot. Host precipitation/cloud physics may use Ex-Aero's size,
  composition, hygroscopicity, and activation diagnostics to form these rates.
  The comparison records CAM-SIMA's per-mode number/volume coefficient
  structure and distinct reservoirs, while keeping its collection kernel,
  precipitation/timestep treatment, and resuspension physics external.
  Ex-Aero's local-rate adapter is a host-integration addition, not an exact
  port of CAM-SIMA's complete wet-deposition algorithm.

Where CAM-SIMA's behavior relies on a fixed mode taxonomy, map that taxonomy
to explicit Ex-Aero configuration rather than inferring roles from names.
Document any behavior that cannot be represented by Ex-Aero's independent
mass/number state without changing the user's selected architecture.

### 3. Modal shortwave/longwave optics data

- Keep the separate RRTMG-30 API and the pinned 100-sample uniform-wavenumber
  recipe, with explicit support metadata for each material/band.
- Make regeneration practical using the pinned GEOSmie implementation (for
  example, optimize or safely parallelize the exact calculation); retain
  deterministic artifacts and validate optimized results against scalar
  GEOSmie reference calculations.
- Regenerate generic modal and MAM4xx artifacts only after the data checks
  pass. Check every generated header against its NetCDF snapshot.
- Compare integrated modal coefficients against the reference optical
  calculations and verify no source-domain extrapolation. Uncovered
  material/band pairs remain explicitly unavailable.

### 4. Bulk analogs and interface consistency

- Check that each applicable modal capability has a bulk counterpart where
  the representation permits it, without inventing prognostic bulk number or
  transplanting modal equations into mass-only state.
- Keep zero-based C/C++ selectors, Fortran layout, units, status codes, and
  error behavior aligned across language wrappers.
- Review MAM4xx forwarding only for shared diagnostics/optics contracts; do
  not expand its native microphysics scope as a side effect of generic-modal
  parity.

### 5. Integrated validation and parity report

- Run focused numerical reference tests, conservation/budget tests, atomic
  failure tests, and C/C++/Fortran boundary tests for each slice.
- Run the complete supported CMake test suite with and without the optional
  MAM4xx backend in a provisioned build environment.
- Publish an updated matrix marking each feature complete, partial, or
  intentionally different, with the tested CAM-SIMA revision and evidence.

## Acceptance criteria

- Every in-scope implemented CAM-SIMA modal capability is either matched in
  Ex-Aero within its stated tolerance or listed as an intentional difference
  with a scientific/API reason and a test protecting that choice.
- No unimplemented CAM-SIMA bulk stub is presented as a feature requirement.
- Ex-Aero bulk and modal diagnostics have representation-appropriate contracts
  and pass parity or invariant tests through C++, C, and Fortran where exposed.
- Supplied zero velocity or scavenging coefficient returns zero tendencies;
  negative/nonfinite rates, nonpositive/nonfinite layer thickness on the
  velocity path, or nonfinite calculated tendency fails atomically without
  changing output buffers. Valid velocity-path tendencies equal
  `-velocity/layer_thickness * state`; valid wet-scavenging tendencies equal
  `-coefficient * state`. Input state remains unchanged. The wet adapter's
  coefficient semantics are validated as an Ex-Aero host contract, not claimed
  as CAM-SIMA numerical parity.
- Different supplied bulk-species/bin and modal-mode coefficients produce
  their corresponding distinct tendencies without cross-talk. Independent
  modal mass and number coefficients are preserved; the common mass
  coefficient within each internally mixed mode preserves its component ratios.
- All source-supported RRTMG modal material/band pairs have deterministic
  generated native data; unsupported pairs remain visible in availability
  results and fail positive-mass requests atomically.
- MKPP chemistry ownership and accepted/rejected integration semantics remain
  unchanged.
- Full-build evidence is recorded separately from syntax and focused unit
  checks; the roadmap is not marked complete on partial-environment checks.

## Open implementation choices

- The parity matrix will identify the exact CAM-SIMA modal routines and
  equations to use; the builder may group related methods to keep commits
  reviewable.
- Reference fixtures should favor source-derived numeric outputs and a direct
  optional CAM-SIMA comparison over copied implementation code.
- Performance work on GEOSmie must preserve its scientific implementation and
  output definition; the builder should profile before selecting a batching or
  parallel execution strategy.
