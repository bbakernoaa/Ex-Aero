# Ex-Aero aerosol parity with CAM-SIMA — design

## Intent

Bring Ex-Aero's GOCART bulk and generic modal aerosol surfaces close to the
working diagnostic and modal-physics behavior in CAM-SIMA, while retaining the
Ex-Aero/MKPP split: MKPP/host code owns chemical formation, equilibrium
thermodynamics, and integration; Ex-Aero supplies aerosol geometry, process
terms, and diagnostics. The comparison target is CAM-SIMA `development` at
`52723a73a2a49e37fdd981f795d08084bb5bfda8`. The Ex-Aero starting point is
`codex/native-aerosol-optics` at `5b347cf` plus this design.

## Goals

- Map every in-scope CAM-SIMA aerosol capability to an Ex-Aero equivalent,
  an intentional representation difference, or a documented remaining gap.
- Match the numerical behavior of implemented CAM-SIMA modal formulas where
  the state definitions and physical assumptions are comparable.
- Add reference cases and cross-interface tests before extending behavior.
- Preserve Ex-Aero's bulk/modal representations and host-owned step ordering.
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
