# Common aerosol host-coupling design

## Intent and agreed direction

Provide a common host interface for bulk and generic modal aerosol helpers.
The user requested population metadata and removal first, followed by common
cloud/ice, activation, physical-property, resuspension and deposition-rebinning
operations. The initialized package selects representation internally. Existing
diagnostics, CCN, emissions and optics interfaces already provide common entry
points and are reused.

The [removal design](2026-10-29-representation-aware-removal-design.md) remains
the detailed mass/number-rate contract within this broader design. No full CAM
wet-deposition algorithm or numerical equivalence is implied: CAM's tracer
solubility factors motivate the flexibility, while the host supplies the rates.

## Population and state-field contract

- A bulk population is one configured species/bin, with diagnosed particle
  number and one mass state field.
- A generic modal population is one configured mode, with one prognostic
  number field and a mass field for each configured material. Each material's
  permission in that mode is explicit, including allocated but prohibited slots.
- Species indices are package-local. Bulk bin names do not imply a universal
  chemical taxonomy; modal material indices refer to configuration order.
- Metadata maps every state field to its population, quantity, constituent
  identity where applicable, and permission. Hosts never parse field names or
  reconstruct packed offsets. Indices are zero-based across C++/C/Fortran.
- Queries return value snapshots valid for the current initialized configuration;
  hosts re-query after initialization changes. No shared mutable global registry.

Use additive `HostCoupling.hpp` metadata types and free functions, without
changing `IAerosolPackage`'s vtable or the existing `StateField` layout:

- `PopulationInfo`: `index`, `name`, `number_kind`, `number_field` (-1 if diagnosed).
- `ConstituentInfo`: `population_index`, `species_index`, `mass_field`, `name`,
  `permitted`. Include every allocated constituent mass slot.
- `HostCouplingInfo`: representation, field count, population descriptors and
  constituent descriptors.
- `hostCouplingInfo(package)` and `supportsHostOperation(package, operation)`.

Operation queries cover property evaluation, cloud/ice diagnostics, activation,
rate/velocity removal, resuspension, deposition rebinning and modal bound
inspection/enforcement. They check actual extension support and initialization;
representation alone is insufficient. Unsupported operations return false and
calls fail explicitly without modifying outputs. Support reports executable
package behavior, separately from state-dependent per-property availability.

## Shared operations and shapes

All arrays are synchronous column-major host-memory views. Let `P` be the
metadata population count and `F` the state-field count.

| Operation | Shared contract |
| --- | --- |
| Removal rates | Mass coefficients/output `(cell,level,F)`; bulk diagnostic number loss per population; modal number coefficients/output `(cell,level,P)`. Bulk number-rate input absent; modal required. |
| Removal velocity | Velocity `(cell,level,P)`, thickness `(cell,level)`; one particle velocity per population maps to its mass constituents and number. |
| Cloud/ice | Existing 12-property slabs `(cell,level,P,12)` and population status; separate ambient and optional host-owned cloudborne snapshots. |
| Activation | Existing output `(cell,level,P,2*num_S)`; number/mass fractions remain distinct. |
| Properties | `(cell,level,P,property_count)` with per-value availability and separate empty/populated status; exact radius definitions and SI units declared per property. |
| Resuspension | Packed material rates `(cell,level,F)`; routes identify source/destination mass fields; read original sources once and preserve mass, with no inferred particle number. |
| Deposition rebinning | Packed deposited budgets plus package-local species index and dry diameter edges; bins and explicit below/above tails. Modal material budgets aggregate across modes; bulk selection identifies one configured bin/species. |

The property set includes mass, number, dry number-median radius, wet
number-median radius, existing dry/wet scavenging diameters, dry volume, dry
surface area, kappa, dry/wet density, lognormal width, dry lower/upper radii,
and dry effective radius. Distinct size definitions remain distinct outputs:
bulk caller-anchored scavenging diameter need not equal its lognormal median.
Configured, derived and unavailable values have explicit status; an unavailable
slot's numerical zero must never be interpreted without that status. Empty
population status is separate from property availability. Reuse current growth
laws and RH caps, independent of optical tables.

## Representation-specific operations

Modal number-bound inspection and enforcement remain explicit extensions.
Bulk physical size bounds are properties, not an invitation to create or clamp
a bulk number prognostic. Modal uptake derivatives, nucleation terms and
chemistry/MKPP mapping retain their current extension interfaces. Native MAM4xx
keeps its distinct layout: only advertise common helpers when a real adapter
exists; unsupported helpers fail atomically.

## Compatibility, ownership and implementation approach

Alternatives considered were replacing the package interface with a large new
virtual interface, extending only removal, or adding a common facade over the
existing implementations. The chosen direction is the additive facade, with
metadata and capability queries first. It supplies consistent host calls while
preserving existing public functions and numerical kernels.

Every delivered operation has C++, C and Fortran interfaces, matching error
semantics, dimensions and index conventions. No hidden state update, timestep
integration, positivity adjustment or bounds correction is introduced.
Precipitation transport, collection laws, species solubility factors, chemistry,
cloudborne prognostics and accepted-state updates remain host/MKPP-owned.
Independent constituent mass rates may change modal composition; the host
supplies a consistent particle-number rate. Whole-particle removal can still use
equal constituent mass rates and the existing mode-shared convenience helper.

Numerical calls validate complete inputs and stage all outputs; invalid shapes,
NaN/Inf, aliasing and overflow leave outputs unchanged. C reports the established
256-byte message and error flag; Fortran uses interoperable wrappers. Metadata
queries do no I/O and introduce no process-global state.

## Delivery and review

The [roadmap](../plans/2026-10-29-host-coupling-roadmap.md) divides work into four
independently testable slices. Each slice includes bindings, documentation and
verification before moving on. This revision changes plans only; implementation
has not begun.
