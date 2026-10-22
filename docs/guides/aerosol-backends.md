# Bulk and modal aerosol backends

The separate configurable `modal` backend is documented in the
[generic modal guide](generic-modal.md), including its generated GEOSmie optics
and MKPP contract. The fixed MAM4xx backend described below has a separate
opt-in adapter for tables matching its native material properties.

EX-aero can select `gocart` (bulk) or the optional, experimental `mam4xx`
(modal) package at runtime. Both use the same column-major
`(cell, level, field)` boundary. GOCART remains the default build and the
existing `exaero_create_gocart_package()` entry point remains available.
Rebuild consuming C++ applications after upgrading the interface library.

## Build

All backends require Kokkos 5.1.1 or newer; CI pins Kokkos 5.1.1 at revision
`267ebc25fc5c8b96bb321f34d78f97b1c30f8830`. CMake rejects older discovered
packages, and compilation also checks the version for host-provided targets.
GOCART uses Kokkos and yaml-cpp. MAM4xx additionally
needs EKAT built against **the same Kokkos**, with its Kokkos utilities enabled.
The pinned MAM4xx revision is `c994f6889f270c3da4aca8a4703a1ebe30f92f69`.
The EKAT revision used by that upstream checkout is
`007bdcff6cb26948517d2d26f3d7ca005363d119`.

```sh
git submodule update --init --recursive
cmake -S . -B build -DEXAERO_ENABLE_MAM4XX=ON \
  -DCMAKE_PREFIX_PATH="${KOKKOS_PREFIX};${EKAT_PREFIX};${YAML_CPP_PREFIX}"
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

The optional build compiles the pinned upstream library sources, reusing the
host's Kokkos and EKAT targets. It does not invoke the upstream standalone
driver, fetch a second Kokkos, or download validation datasets at configure
time. An existing parent-provided `mam4xx` target is also accepted; its
precision must be double. With `EXAERO_ENABLE_MAM4XX=OFF`, neither EKAT nor
initialized MAM4xx sources are required.

## Select a package

```cpp
#include <exaero/PackageFactory.hpp>

auto package = exaero::createPackage("mam4xx");
package->initialize("backend: mam4xx");
const auto fields = package->stateFields();
```

For C/Fortran hosts, use `exaero_create_package("gocart", ...)` or
`exaero_create_package("mam4xx", ...)`. A disabled or unknown backend returns
a null handle and an error; it never falls back to GOCART. Query
`exaero_package_available`, `exaero_get_representation`,
`exaero_get_capabilities`, `exaero_get_state_field_count`, and
`exaero_get_state_field_info` before allocating buffers or calling optional
operations. `exaero_execute_microphysics` advances either backend through the
generic interface. These routines are mirrored in `exaero_interface.f90`.

## Packed modal state and units

GOCART keeps one mass concentration per configured bulk species. MAM4xx has
64 fields: 32 interstitial fields followed by 32 cloudborne fields. Each
phase has four modes in upstream order: accumulation, Aitken, coarse, primary
carbon. Each mode has seven mass fields in **global AeroId order** (SOA, SO4,
POM, BC, NaCl, DST, MOM), then one independently prognosed number field.
This is not the mode-local species ordering used by some other upstream
MAM4xx processes.

Names are `<species>_<mode>_<i|c>` and `num_<mode>_<i|c>`, for example
`SO4_accumulation_i` and `num_aitken_c`. Use the name queries or the constants
and index functions in `exaero/Mam4xxConfig.hpp`.

Mass is **kg/m³** and number is **particles/m³**, including cloudborne fields
expressed per volume of air. Number must never be treated as another mass
species. The full species axis preserves material transferred by coagulation
into coated primary-carbon particles. Hosts should initialize species absent
from a mode to zero. A populated mode must have both positive total mass and
positive number; an empty mode has both zero. Invalid shapes, negative or
nonfinite values, and inconsistent mass/number pairs fail before modifying
the caller's state. Zero timestep leaves the state unchanged.

The public views refer to **host memory**. The MAM4xx solver stages them to
the configured Kokkos execution space and copies results back synchronously.
Initialize the EX-aero/Kokkos environment before computing and release all
packages before finalizing it. GPU execution has not yet been validated.

## Supported physics and current limits

| Surface | GOCART bulk | Experimental MAM4xx modal |
| --- | --- | --- |
| Microphysics | Existing bulk step | Upstream clear-air coagulation |
| Wet sizing | Existing GOCART sizing | Upstream MAM4xx water uptake and hysteresis |
| Diagnostics | Existing diagnostics | Both phases: mass, number, water, area, size cuts, settling, columns |
| CCN spectrum | Existing bulk activation | Modal lognormal kappa-Kohler integral |
| Emissions | Configured raw-species mapping | Already-speciated 64-field mass/number rates |
| Optics | Existing GOCART/GEOSmie optics | Opt-in GEOSmie component-table approximation |

MAM4xx initialization accepts an empty YAML mapping, `backend: mam4xx`, and
the optional `optics` mapping described below; unknown settings are rejected. The timestep updates
interstitial coagulation only. Cloudborne state is carried unchanged;
diagnostics include both phases. Gas chemistry, nucleation, condensation,
aging, mode renaming, cloud activation tendencies, and wet/dry deposition
are not yet connected. This is not a complete MAM4 atmosphere driver or
production scientific validation.

Dry geometric sizes use each mode's mass-derived volume, independent number,
and upstream lognormal width. Wet sizes use upstream material properties and
the water-uptake kernel, with RH capped at 0.98 as in upstream MAM4xx.
PM1 is a dry geometric diameter cut. PM2.5/PM10 use wet aerodynamic diameter
with a spherical density correction; slip and nonspherical corrections are
omitted. Settling uses mass-weighted lognormal Stokes velocity and Sutherland
air viscosity, omitting slip correction. Column and surface diagnostics are
written at level 0, which is assumed to be the surface.

CCN uses the lognormal kappa-Kohler integral underlying MAM4's
`ndrop::ccncalc`, extended to the caller's supersaturation axis. It returns
particles/m³, includes both phases, and does not impose upstream's
fixed-axis number clipping. This is a spectrum diagnostic, not a cloud
droplet tendency. Optical properties use the explicit approximation below.

For modal emissions, input/output field order is identical to the state.
Concentration rates are kg/m³/s or particles/m³/s according to field metadata.
Area fluxes are kg/m²/s or particles/m²/s and are divided by positive layer
thickness. This does not split a bulk inventory into modes or infer emitted
number. Hosts must supply both mass and number consistently.

## MAM4xx component optics

Enable the shared generic-modal integrator explicitly:

```yaml
backend: mam4xx
optics:
  treatment: component_tables
  include_cloudborne: false
  # Optional source aliases; names use the public global species order.
  sources: {SOA: MAM4_SOA, SO4: MAM4_SO4, POM: MAM4_POM, BC: MAM4_BC,
            NaCl: MAM4_NACL, DST: MAM4_DST, MOM: MAM4_MOM}
```

The adapter reads native MAM4xx densities, hygroscopicities and fixed mode
widths, and integrates the dedicated GEOSmie particle tables over each mode's
mass/number-derived dry lognormal distribution. It shares the existing 64-field
state without repacking or modifying it. Native microphysics and diagnostics
continue to use MAM4xx. Reinitializing without `optics` disables optical support;
failed initialization leaves the previous configuration intact.

The seven `MAM4_*` curves are embedded in `Mam4xxOpticsData.hpp`, independently
of the generic `MODAL_*` defaults and legacy GOCART tables. Their
[recipe](../../tools/modal_optics/mam4xx_manifest.json) pins GEOSmie and MAM4xx,
records density and kappa, and **labels SOA, POM and MOM as configurable OPAC
`waso00` refractive-index surrogates**. Change those recipe entries and regenerate
for alternative organic materials. These surrogates do not constitute validated
MAM4 organic optical constants. Dust is spherical. Native MAM4 has no nitrate
state field; nitrate remains available to the generic modal backend.

This is a component-table optical approximation: components share their mode's
dry radius distribution, and their extinction/scattering are added by dry mass.
It does not reproduce the native internally mixed MAM optical calculation.
Each table includes its material's equilibrium kappa growth and mixing with
water exactly once. Native MAM water diagnostics instead use mixture kappa and
hysteresis, so their wet sizes need not match this optical approximation.
No diagnostic wet radius or water mass is applied again in optics.

Default optics include interstitial aerosol only. `include_cloudborne: true`
adds the carried cloudborne particle distributions as ambient aerosol; it does
not compute cloud-droplet optics. Tables cover dry radii 0.1 nm–200 micrometres,
RH 0–0.98 and five wavelengths 355–1640 nm. The integrator rejects an unsupported
RH/wavelength or a distribution with more than 1e-8 of its dry mass outside the
radius range. It does not extrapolate. These are visible/NIR starter tables,
not full RRTMG bands. The generic guide documents output slots and integration.

Optional `sources` entries select alternate curve labels; an optional `data_file`
is merged through the existing EX-aero MIE loader before resolving them. Every
curve must carry recognized particle/distribution semantics, physical wavelength
nodes and compatible density/kappa metadata. Legacy bulk tables cannot be
silently substituted. For the initial embedded set, configure refractive indices
in the offline recipe and regenerate the snapshot/header together.

## Validation

The modal tests compare dry-air coagulation with the native pinned kernel
over multiple columns and levels, check mass conservation and number
evolution, preserve cloudborne state, verify column integrals and CCN bounds,
and exercise C/Fortran selection and failure paths. Existing GOCART tests
run in both build configurations. These are software and limited physics
checks, not a substitute for reference-case scientific validation.
