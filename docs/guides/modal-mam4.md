# Configurable modal aerosols and the MAM4 reference profile

The `modal` backend now has editable MAM4 mode/material presets and independent
process choices. This is the first validated step toward a configurable MAM4
representation. The native reference profile matches the **clear-air coagulation
and water diagnostics currently exposed by EX-aero's MAM4xx wrapper**. It is not
a complete MAM4 atmospheric integration.

## Select a profile

A general four-mode setup, available with or without the optional MAM4xx build:

```yaml
backend: modal
preset: mam4
processes:
  coagulation: brownian
  water_uptake: kappa_equilibrium
```

This expands to accumulation, Aitken, coarse and primary-carbon modes, with
native reference geometric radii and widths, and SOA/SO4/POM/BC/NaCl/dust/marine
organic materials. The preset's Brownian routes explicitly select destinations
for all six intermode pairs. These generic routes differ from the native MAM4
coagulation parameterization, which uses three selected mode pairs and its own
analytic update. No growth-transfer or nucleation calibration is inferred.

For the native reference processes, build with `EXAERO_ENABLE_MAM4XX=ON`:

```yaml
backend: modal
preset: mam4xx_reference
```

Equivalent process choices on the `mam4` layout are:

```yaml
backend: modal
preset: mam4
processes:
  coagulation: mam4xx
  water_uptake: mam4xx
```

The profile uses the pinned MAM4xx coagulation kernel, including its native
number/concentration conversion and wet-size assumptions. Cloudborne fields in
the private native bridge are zero: the generic public state remains **32
interstitial fields**. There is no implicit cloudborne state or cloud model.

Native coagulation requires native mode order, widths, reference radii,
water thresholds, material order, densities and kappa. Incompatible changes
fail at initialization. It uses fixed native routes, so explicit
`coagulation_routes` are rejected for this choice. Switching to `brownian`
restores arbitrary configurations and requires explicit routes for every pair.
`processes.coagulation: none` disables coagulation without disabling explicitly
configured growth transfer. Native coagulation, when enabled, runs before
optional generic growth-tail transfer; that combination is a declared hybrid.

## Edit the representation

Change individual mode parameters without repeating a preset's mode list:

```yaml
backend: modal
preset: mam4
mode_overrides:
  accumulation: {geometric_radius: 6e-8, geometric_std_dev: 1.7}
  aitken: {crystallization_rh: 0.2, deliquescence_rh: 0.7}
```

Water thresholds affect the native water choice. Native coagulation rejects
changes incompatible with its fixed widths/thresholds. Unknown mode names or
settings fail explicitly.

C++ callers can obtain a complete editable configuration:

```cpp
#include <exaero/ModalPresets.hpp>
#include <exaero/ModalCoupling.hpp>
#include <exaero/PackageFactory.hpp>

auto config = exaero::mam4ModalConfig();
config.modes[0].geometric_std_dev = 1.7;
config.coagulation = false;
config.coagulation_routes.clear();
config.modes.push_back({"extra", 4e-8, 1.4, {}});
auto package = exaero::createPackage("modal");
exaero::initializeModal(*package, config);
```

`mam4xxReferenceModalConfig()` returns the native reference choices. In YAML,
explicit `species` or `modes` sequences replace the corresponding preset list.
Replacing modes clears inherited routes; callers must supply routes naming the
new modes. Explicit route lists replace inherited route lists. Unknown settings
and conflicting legacy `coagulation`/`processes.coagulation` settings fail.

Each mode retains every configured material slot. The presets permit coatings
in every mode, including primary carbon, so coagulation cannot drop incoming
sulfate/organic mass before an aging operation is connected. They do not impose
the native model's mode-local emission/species inventory. Hosts should speciate
emissions explicitly. Query field metadata before allocating state arrays.

Preset molecular weights are explicit host chemistry surrogates in kg/mol.
MOM uses 0.250093 kg/mol; native coagulation/water do not use molecular weights.
MKPP must supply a consistent mechanism and mass basis. The generic empirical
pure-H2SO4 nucleation helper is not automatically calibrated or substituted for
native MAM4 nucleation; sulfate's molecular-weight convention differs.

## Physical constants

`constants` configures `gravity` [m/s2], `water_molecular_weight` [kg/mol],
`gas_constant` [J/mol/K], and `surface_tension` [N/m]. All must be finite and
positive. Defaults preserve the existing generic diagnostic/transport constants.
Both MAM4 presets explicitly retain the pinned native constants for settling
and CCN, including native water molecular weight and gas constant. These values
are part of process-parity comparisons. The native water/coagulation kernels
retain their own upstream constants; changing this block affects generic
settling, CCN and dry uptake, not those native kernels.

## Water uptake independently of mode count

`water_uptake: mam4xx` calls the native scalar Kohler solver for each generic
mode, using its volume-weighted mixture kappa and mean-volume dry radius.
It supports arbitrary material/mode counts and configurable
`crystallization_rh`/`deliquescence_rh` per mode (defaults 0.35/0.8). The valid
range is `0 <= crystallization < deliquescence <= 1`. RH is capped at 0.98.
The native linear transition prescription is used; it does not introduce a
prognostic wet/dry history flag. Native solver constants and its 30-micrometre
wet mean-radius cap are retained. There is no second radius-growth correction.

Native wet-property Jacobians are currently unavailable and requests fail
before outputs change. Dry-geometry uptake Jacobians remain available for MKPP.
`kappa_equilibrium` retains the existing analytic property and uptake Jacobians.
The native option requires the MAM4xx dependency even with generic coagulation;
unavailable processes fail explicitly rather than falling back.

## Optical and chemistry scope

Both profiles bind the embedded `MAM4_*` GEOSmie component tables. Density/kappa
must match each recipe. Changing optical material properties requires compatible
tables or removing/replacing their bindings. SOA, POM and MOM remain labeled OPAC
organic surrogates. Components include their own equilibrium humidification;
the native mixture-water diagnostics may therefore differ from optical growth.
This is a component-table optical approximation, not native internally mixed
MAM4 optics. Native property Jacobian limitations do not affect optics queries.

MKPP continues to own gas chemistry and condensation/nucleation integration.
Aging, native growth renaming, cloud activation tendencies, scavenging, and
deposition remain unconnected. The next scientific work is to connect these
processes with explicit host gas/cloud inputs and independent reference cases;
none are silently approximated by enabling a preset. Existing generic growth
transfer, CCN diagnostics, emissions and coupling interfaces remain available.

## Validation and build

Tests compare all standard diagnostic slots over multiple columns/levels and
RH regimes, and compare coagulation state updates against the fixed MAM4xx
wrapper. Additional tests exercise arbitrary mode/material counts, dry-uptake derivative finite differences, conservation,
C ABI preset initialization, process guards and transactional failures.

With Kokkos 5.1.1, yaml-cpp, Catch2 v2 and Python data dependencies installed:

```sh
cmake -S . -B build -DEXAERO_BUILD_CHEMISTRY=OFF \
  -DEXAERO_ENABLE_MAM4XX=ON \
  -DCMAKE_PREFIX_PATH="${KOKKOS_PREFIX};${EKAT_PREFIX};${YAML_CPP_PREFIX};${CATCH2_PREFIX}"
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

Also build with `EXAERO_ENABLE_MAM4XX=OFF` to test the dependency-free generic
profile and explicit rejection of unavailable native processes. These software
and limited process-parity checks do not establish full MAM4 scientific
validation. GPU execution still requires separate validation.
