# Common host-coupling execution record

## Delivered scope

All four [roadmap slices](2026-10-29-host-coupling-roadmap.md) are implemented:
population/constituent metadata and operation capabilities, common rates and
velocities, cloud/ice and activation, physical population properties, and
resuspension/deposition utilities. Every operation has C++, C and Fortran
surfaces. Host/MKPP ownership and the diagnosed bulk / prognostic modal number
contracts are preserved. See [usage and verification](../../guides/host-coupling.md).

Implementation commits before final review: `3214390`, `d0ebb05`, `0f5320f`,
`0331f33`, `47542ad`, `390d5bf`; planning commit `fe83d4a`. Work is local on
`codex/native-aerosol-optics`, in the pre-existing linked worktree.

## Verification evidence

Tests were written before each new surface: missing C++ headers/functions and
missing Fortran bindings failed first, followed by passing implementations.
Regression tests independently check analytical mass/number losses, lognormal
moments, composition kappas 0.45/0.15, explicit radius definitions, route chains,
CDF bin/tail closure, and atomic error behavior. Configuration reordering and
initialized native-backend rejection are exercised.

Before review: generic Debug 29/29 and optional MAM4xx Debug 31/31 passed.
The independent reviewer confirmed those counts against logs. GNU C++/Fortran
12.2, Kokkos 5.1.1 Serial; Fortran uses `-fcheck=all -fbacktrace`.
Eleven facade/removal implementation translation units pass
`-Wall -Wextra -Wpedantic -Werror` syntax checks. Public C headers pass strict
C11 compilation; new public headers build in the standalone compile gate.

After review, all three finding classes have regression tests that failed
before the fixes and pass after them. Post-fix complete suites pass: generic **29/29** (13.62s), optional MAM4xx
**31/31** (28.74s). Strict checks for the eleven facade/removal implementation
translation units also pass after the fixes. Commands were
`cmake --build build -j 2 && ctest --test-dir build --output-on-failure` and
the same commands with `build-mam4xx`. No GPU, MPI, full CAM model or coupled
host-rate-law validation is claimed.

## Independent review and fix pass

A fresh reviewer examined `53a833f..390d5bf`, read-only, against both specs and
all four plans. It found no Critical issues, two Important issues, and one
Minor dimension issue. The executor treated the dimension rejection as a
supported-input correctness issue and fixed it in the same pass.

- Environmental buffer overlap: common cloud/activation forwarding omitted
  legacy environment/output checks; properties omitted optional environment
  buffers. New regression tests reproduced accepted aliases. The facade now
  protects every supplied environmental view against every output; optional
  pressure/density/thickness remain optional. C error paths preserve buffers.
- Failed bulk initialization: species names were visible before all GOCART
  validation completed. Fresh/previously initialized bulk-package failures now
  hide partial field metadata until a successful initialization. Tests cover
  late YAML/config errors, malformed YAML, operation support and recovery.
- Rebin dimensions: output size incorrectly included the population count.
  A dedicated checked rebin-output extent accepts a 100-bin, 300,000-cell
  configuration without allocating its grid, while rejecting real overflow.

No second review was dispatched; the required fix pass is verified through
failing-then-passing regressions and complete suites. No code findings deferred.

## Rulings made during execution

1. Use one umbrella ledger and manually mirror task records because cloud skill
   scripts are resources rather than local executables. Cost if wrong: execution
   tooling only; no scientific or API effect.
2. Use a dedicated initialized-package dispatcher test executable instead of
   duplicating its fixtures in old bulk/modal tests. Cost: tests have a different
   file location; both representations still execute.
3. Commit the common dispatcher and bindings together because their integration
   test exercises C++ and C under one nullable-number-rate contract. Cost: coarser
   commit granularity; no runtime behavior change.
4. Put common cloud C wrappers in a separate implementation file to avoid
   anonymous-namespace collisions and preserve old wrappers. Cost: one extra
   source entry; the public interface is unchanged.
5. Reuse a multi-population Fortran driver for cloud, property and deposition
   assertions. Cost: a shared fixture; per-slice assertions remain distinct.

6. Failed GOCART initialization invalidates exposed helper metadata until a
   successful retry because its initialization mutates members before validation
   completes. Generic modal initialization remains transactional. Cost: hosts
   must retry failed bulk initialization before calling common helpers.

## Rulings on areas the reviewer set aside

- Full CAM wet-deposition equivalence, precipitation and accepted integration
  remain host-owned by explicit user choice. Cost: host integration/rate laws
  still need their own validation.
- Existing activation thresholds and cloud/ice weights are preserved, rather
  than recalibrated. Cost: their pre-existing scientific limitations remain.
- GPU/MPI and production throughput are outside the Serial verification used
  here. Cost: those execution environments remain unverified for this facade.
- Native water-uptake science is delegated to the existing implementation;
  native-package facade helpers remain explicitly unsupported. Cost: no new
  independent validation of that underlying scientific model.
- Deliberately inconsistent third-party extension metadata is not a supported
  package contract; real built-in initialization failures were fixed. Cost:
  third-party implementations must obey their metadata contracts.
- Descriptor scalar/name output aliasing has no numerical-array non-overlap
  guarantee. Cost: callers must provide separate descriptor storage; atomic
  error behavior and name-capacity checks remain tested.
- Fabricated C++ views beyond addressable storage are outside the valid view
  contract. Cost: callers remain responsible for real buffer storage; practical
  dimension overflow guards and regressions are covered.
