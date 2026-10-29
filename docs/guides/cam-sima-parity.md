# CAM-SIMA aerosol behavior crosswalk

Reference: [CAM-SIMA 52723a73](https://github.com/ESCOMP/CAM-SIMA/tree/52723a73a2a49e37fdd981f795d08084bb5bfda8/src/aerosol)
and its pinned [atmospheric_physics 6ee9b89a](https://github.com/ESCOMP/atmospheric_physics/tree/6ee9b89afb294408529bea5b72c7bab47c0c55eb).
Ex-Aero baseline: `5b347cf`, with the subsequent parity design. The cases below
are derived from the pinned source equations, not from an executed CAM-SIMA
model. Numerical tolerances are relative `1e-12` unless specified otherwise;
zero-valued budgets and reserved number slots are checked exactly.

| Capability / source routine | Ex-Aero surface | Status and representation contract |
| --- | --- | --- |
| Modal `get`, `species_type`, mode width and constituent properties | `modalConfig`, `stateFields`, `modalProperties` | Equivalent property access with explicitly configured material names, density and kappa. Species mass slots and independent number are packed per mode. |
| Modal `amcube` / mode geometry | `modalProperties`, `modalRadius` | Equivalent lognormal volume/number relation. Ex-Aero uses radius in m; CAM mode configuration uses diameter in m. |
| Modal `apply_number_limits` | `inspectModalBounds`, `enforceModalBounds` | Equivalent volume-derived number interval when diameter/radius bounds are converted. Explicit host call corrects number and preserves every material mass. Empty modes stay empty. |
| Modal `actfracs` | `modalActivationFractions` | Equivalent separate zeroth-moment number and third-moment mass tails. Host supplies supersaturation; activation and transfer integration remain external. |
| Modal `scav_diam` | `modalScavengingSizes` | Intentional difference: CAM returns configured `dgnum(m)`; Ex-Aero returns current dry/wet number-median diameters from state. The CAM accessor's cm comment conflicts with the consuming wetdep code's explicit m-to-cm conversion; its actual property is in m. |
| `aero_convproc` composition-weighted hygroscopicity | `modalConfig`, `modalProperties` kappa, wet size and density | Equivalent volume-weighted constituent kappa for the generic model. CAM's POM-specific override is a host closure, not an Ex-Aero default. |
| Modal `soluble` / `hydrophilic` | configured kappa, mode composition and host wet-removal policy | Intentional difference: CAM tests fixed mode taxonomy (`primary_carbon` / `accum`). Ex-Aero does not infer wet-removal efficiency from names or reuse ice roles. Hosts can access mode/material metadata and set their own activation/solubility factors. |
| `wetdep::get_bcscavcoefs`, `init_bcscavcoef` | size/composition diagnostics for a host scheme | Intentional difference in ownership: implemented CAM physics stays external. Distinct per-mode number/volume moments are preserved by Ex-Aero's host-rate adapter contract. |
| `wetdep::wetdepa_v2` | `bulkRemovalRateTendencies`, `modalRemovalRateTendencies` | Implemented host integration addition, not an equivalent wet-deposition scheme. Host supplies effective nonnegative /s sink coefficients; CAM's precipitation-normalized /mm coefficients are not interchangeable. Precipitation transport, step limiting, evaporation/resuspension and surface budgets stay host-owned. |
| Host dry-deposition velocity | `bulkRemovalVelocityTendencies`, `modalRemovalVelocityTendencies` | Implemented host integration addition: one velocity per species/bin or mode; return `-v/dz * concentration`. Host chooses active levels and applies accepted tendencies. |
| Ambient/cloudborne state getters | `bulkCloudIce`, `modalCloudIce`, host-owned snapshots | Equivalent reservoir separation with different Ex-Aero units/layout. Missing cloudborne input returns explicit `Absent`; no inferred transfer or double removal. |
| Modal `icenuc_updates_*`, `hetfrz_species`, category/size weights | configured `IceCategory` / `IcePopulationRule`, cloud/ice diagnostics | Intentional difference: configurable eligibility and available number, without CAM-specific mode-name rules or an Ex-Aero ice-nucleation/transfer integrator. |
| Modal `resuspension_resize` | `redistributeModalResuspension` | Partial / configurable generalization: CAM sends each material's returns to its largest participating mode. Direct routes reproduce that mapping; other explicit route graphs are intentional host choices. Read original source rates once; infer no particle number. |
| Modal `rebin_bulk_fluxes` | `rebinModalDeposition` | Equivalent single-material lognormal mass fractions; intentional extension for mixture-volume geometry, selected configured material, stable tails and explicit below/above budgets. CAM's type-density approximation is not copied into mixtures. Both use deposited mass and number to infer distribution size. |
| Bulk fixed number and size metadata | `bulkFixedNumber`, `bulkReferenceRadii` | Equivalent fixed-lognormal mass-to-number conversion when physical parameters agree. Number remains diagnostic and radii are independent of optical tables. Zero-kappa material contributes to ordinary number/PM/surface diagnostics while remaining inactive in the existing bulk kappa activation model. |
| Bulk `apply_number_limits`, `icenuc_updates_*` | mass-only bulk state, diagnostic cloud/ice helpers | Not applicable as prognostic-number correction / aerosol transfer. CAM bulk number-limit routine is a no-op and ice updates return false. Ex-Aero's configured diagnostic candidates are a host-facing extension. |
| Bulk `actfracs`, `scav_diam`, `resuspension_resize`, `rebin_bulk_fluxes` | bulk activation/size/route/rebin utilities | Ex-Aero extensions: CAM bulk methods are explicitly unimplemented; they are not numeric parity targets. |
| Modal SW/LW native band datasets | `modalBandOptics`, availability APIs | Partial; the independent optics-data plan supplies source-supported bands. This behavior work does not claim new spectral coverage. |

## Source-derived numerical cases

For width `sigma`, the mean particle volume is
`(4*pi/3)*r_g^3*exp(4.5*log(sigma)^2)`. At fixed dry volume, clamping a
mode from `r_g=4e-7 m` to upper radius `2e-7 m` increases number by exactly
a factor of eight. Below-minimum, exact-boundary, in-range and empty cases
must preserve mass and report the original bound status.

CAM `actfracs` uses `x=(2/3)*log(Scrit/S)/(sqrt(2)*log(sigma))` and
`y=x-1.5*sqrt(2)*log(sigma)`, then `fn=erfc(x)/2`, `fm=erfc(y)/2`.
At `S=Scrit`, `fn=0.5` and
`fm=(1+erf(1.5*sqrt(2)*log(sigma)))/2`; these are different moments.

For kappa equilibrium sizing, dry/wet diameter ratio is
`(1+kappa*RH/(1-RH))^(1/3)` before the documented model RH cap. At fixed
total dry volume/number, changing constituent kappa `{0.6,0}` volume fractions
from `{0.75,0.25}` to `{0.25,0.75}` changes effective kappa from `0.45` to
`0.15`, preserving dry size and changing wet size. The host can therefore
compute composition-dependent collection rates; per-mode input does not mean
a composition-independent coefficient.

For single-material deposited mass `M`, number `N` and density `rho`, CAM's
volume-median diameter is `(6*M/(pi*rho*N))^(1/3)*exp(1.5*log(sigma)^2)`.
Integrating the third-moment lognormal tail matches Ex-Aero's radius-based
fractions. Bins plus both tails sum to the input material budget. Changing
number at fixed mass changes bin fractions; scaling both mass and number by
the same factor preserves fractions and scales the integrated mass budgets.

The new local sink contract is checked independently of CAM's rate law:
bulk masses `{2e-9,3e-9} kg/m3` and coefficients `{5e-4,2e-4} /s` produce
`{-1e-12,-6e-13} kg/m3/s`. Within one modal population, masses
`{2e-9,3e-9}`, common mass coefficient `5e-4 /s`, number `5e7 /m3` and
number coefficient `1e-4 /s` produce `{-1e-12,-1.5e-12} kg/m3/s` and
`-5e3 /m3/s`. Material ratios are preserved while the two moments can differ.
A host velocity `0.01 m/s` over a `10 m` layer corresponds to `0.001 /s`.

## Validation record

The generic Debug profile passes all 22 CTest tests with GNU C++/Fortran 12.2,
Kokkos 5.1.1 Serial, yaml-cpp 0.8.0, Catch2 2.13.10 and Python 3.11.2.
The optional Debug profile passes all 24 tests, including native MAM4xx and
the generic/native process comparison, using MAM4xx
`c994f6889f270c3da4aca8a4703a1ebe30f92f69` and its EKAT
`007bdcff6cb26948517d2d26f3d7ca005363d119` dependency (2.1.0).
Both profiles enable Fortran runtime checks (`-fcheck=all -fbacktrace`). The
four changed removal implementation translation units also pass strict
`-Wall -Wextra -Wpedantic -Werror` compilation. The optional public band
forwarding case runs in a separate Kokkos-free test translation unit.
The new source-derived cases pass without changes to modal bounds, equilibrium
wet sizing, activation fractions, resuspension routing, or rebinning formulas.
The adapter tests cover independent population/moment rates, empty populations,
zero sinks, invalid late-array inputs, aliases, overflow and atomic outputs.
The dedicated Fortran removal test runs both backends over two cells and three
levels; cloud/ice Fortran tests check multiple species/property slabs.

The first full run exposed a misplaced zero-kappa guard in bulk diagnostics
and its absence from aggregate CCN activation. Correcting its placement retains
hydrophobic material in diagnostics while preserving zero-kappa nonactivation.
Existing test doubles were updated for the preceding band-optics interface;
an old optics test was aligned with the source-declared physical-wavelength axis.

No full CAM-SIMA model execution, GPU validation, or wet-deposition-rate
equivalence is claimed by these source-derived cases. Optical dataset coverage
is unchanged and remains a separate workstream. See [host removal](host-removal.md)
for units, buffer layouts and the host integration boundary.
