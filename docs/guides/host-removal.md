# Host-owned aerosol removal

Ex-Aero translates **host-computed** removal inputs into signed concentration
tendencies. The host owns precipitation, collection efficiencies, deposition
velocities, integration, positivity, accepted/rejected steps, and state updates.
The adapters perform no implicit modal bounds correction.

Use the C++ functions in `BulkRemoval.hpp` / `ModalRemoval.hpp`, their
`exaero_*_removal_*_tendencies` C equivalents, or the matching
`exaero_interface` Fortran bindings. All views are synchronous host-memory,
column-major arrays. These generic modal helpers do not use the native
`mam4xx` package's different state layout.

## Diagnose properties and calculate host coefficients

Bulk physical size, density, kappa and fixed-distribution number are exposed
through `bulkSpeciesSize`, `bulkFixedNumber`, `bulkReferenceRadii` and
`bulkScavengingSizes`, independently of radiation tables. For modal aerosol,
`modalConfig` exposes constituents and mode width; `modalProperties` supplies
current effective kappa, wet density, volume, number, and dry/wet radius.
`modalScavengingSizes` returns current number-median diameters. Cloud/ice and
activation helpers retain separate ambient and host-owned cloudborne snapshots.

The host combines those properties with its meteorology to compute removal
coefficients. A per-mode coefficient depends on **current species composition**
and size; it is not necessarily a fixed mode parameter. Ice-candidate roles are
not substitutes for host solubility/activation/scavenging factors.

## First-order rate inputs

`bulkRemovalRateTendencies(package, state, mass_rate, mass_sink, number_sink)`
accepts nonnegative `/s` coefficients per species/bin and returns
`dM_s/dt=-lambda_s*M_s`. Diagnostic number loss uses that species' configured
mean lognormal particle mass,
`rho*pi*Dg^3*exp(4.5*log(sigma)^2)/6`. There is no bulk number state field.

`modalRemovalRateTendencies(package, state, mass_rate, number_rate,
mass_sink, number_sink)` accepts independent mass and number coefficients per
mode. It returns `dM_ms/dt=-lambda_mass_m*M_ms` for every material s in mode m
and `dN_m/dt=-lambda_number_m*N_m`. Uniformly mixed particles lose their
constituents together, preserving component mass ratios. Different modes can
have different coefficients. Mass and number moments may differ because of
size-selective scavenging; the host derives a physically consistent pair.

| Buffer | Bulk extent | Modal extent |
| --- | --- | --- |
| State | `(cell,level,species/bin)` | Existing packed `(cell,level,field)` |
| Mass coefficient, /s | `(cell,level,species/bin)` | `(cell,level,mode)` |
| Number coefficient, /s | Not supplied; fixed distribution | `(cell,level,mode)` |
| Mass tendency, kg/m3/s | State extent | State extent; number slots zero |
| Number tendency, particles/m3/s | State extent | `(cell,level,mode)` |

For masses `{2e-9,3e-9} kg/m3`, number `5e7 particles/m3`, mass coefficient
`5e-4 /s` and number coefficient `1e-4 /s`, a mode returns mass tendencies
`{-1e-12,-1.5e-12} kg/m3/s` and number tendency `-5e3 particles/m3/s`.
Add the separate number result only to the mode's number field: the packed
mass output contains zero in that slot.

Call separately with appropriate rates for ambient/cloudborne snapshots.
Bulk number in either snapshot follows the same configured fixed distribution
as the existing bulk cloud/ice diagnostics. CAM-SIMA's precipitation-normalized
scavenging coefficients in `/mm` cannot be passed directly as `/s` inputs here.
The host must derive the effective local sink rates. Precipitation transport,
evaporation/resuspension and surface budgets remain part of the host scheme.
Existing route/rebin helpers transform host-supplied returns/deposited budgets.

## Velocity inputs

`bulkRemovalVelocityTendencies(package, thickness, state, velocity,
mass_sink, number_sink)` and its modal counterpart use `-v/dz` times each
population's mass and number. Velocity is nonnegative m/s per species/bin or
mode; thickness is positive m with extent `(cell,level)`. For example,
`0.01 m/s` over `10 m` gives `0.001 /s`. This simple path applies the same
effective velocity to both moments. Hosts needing distinct moment losses use
the rate API. The host supplies zero velocity at inactive levels.

## Errors and accepted steps

Both outputs must be disjoint from all inputs and each other. Wrong shapes,
nonfinite/negative coefficients or concentrations, invalid modal mass/number
pairs or material permissions, nonpositive/nonfinite thickness, and nonfinite
calculated tendencies fail without changing either output buffer or state.
C++ throws a `FATAL ERROR:` exception; C/Fortran set nonzero `errflg` and
write a 256-byte `errmsg`. Zeroing a coefficient disables only that sink.

The host may reject a step or integrate accepted tendencies. Its explicit
`enforceModalBounds` call can then correct number and return a separate
correction budget; no removal helper adjusts number or budgets implicitly.

## Common representation-aware entry points

The [common host-coupling API](host-coupling.md) provides metadata-driven rate
and velocity adapters for both representations. Its modal rate path accepts a
mass coefficient for each constituent field, plus a separate number coefficient
per mode. The existing `modalRemovalRateTendencies` remains the convenient
mode-shared whole-particle path. All precipitation physics, solubility factors
and accepted updates remain host-owned.
