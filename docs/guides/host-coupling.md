# Common host aerosol coupling

Initialize a package, then query `hostCouplingInfo(package)` and
`supportsHostOperation(package, operation)` from `HostCoupling.hpp`. Metadata
maps every constituent mass field to a population and package-local species,
including allocated but prohibited modal slots. Re-query after reinitialization.
No host should reconstruct packed offsets or parse field names.

A bulk population is one configured species/bin, with diagnostic number and
`number_field=-1`. A modal population is a mode with constituent mass fields
and one prognostic number field. Species indices, population indices and state
field indices are distinct and zero-based in C++, C and Fortran.

## Removal

Use `removalRateTendencies(package,state,mass_rate,number_rate,mass_sink,number_sink)`
from `Removal.hpp`. State, mass rates and mass sinks use `(cell,level,field)`;
number sinks use `(cell,level,population)`. For bulk, pass `nullptr` for the
number-rate view: diagnostic number loss follows mass loss. For modal, pass a
pointer to a `(cell,level,population)` number-rate view and set all packed
number-field mass coefficients to zero. Each constituent mass coefficient can
vary independently. Equal constituent coefficients represent whole-particle
mass removal; the host supplies physically consistent particle-number loss.

`removalVelocityTendencies(package,thickness,state,velocity,mass_sink,number_sink)`
uses one velocity per population. Velocity is m/s, thickness m, coefficients
inverse seconds, and outputs kg/m3/s and particles/m3/s. The host sets zero
velocity at levels where removal is inactive. See [host removal](host-removal.md)
for the external precipitation and integration responsibilities.

The matching C names are `exaero_removal_rate_tendencies` and
`exaero_removal_velocity_tendencies`. Supply `nc,nl,nf,np` from the host grid and
metadata. In Fortran `use exaero_interface`; pass `c_null_ptr` for bulk number
rates or `c_loc(number_rate)` for modal rates in a contiguous `target` array.

## Metadata and error contract

The C/Fortran metadata surface provides `exaero_host_coupling_counts`, indexed
`exaero_host_population_info` / `exaero_host_constituent_info`, and
`exaero_supports_host_operation`. Names use caller-owned buffers and explicit
byte capacities including the terminating NUL. Query counts before sizing arrays.
Known unsupported operation queries succeed with false; invalid operation codes
or descriptor indices set an error and preserve descriptor outputs.

Every numerical facade uses synchronous column-major host arrays, validates the
complete request, requires outputs disjoint from inputs and each other, and
publishes outputs atomically. No call applies a timestep or mutates aerosol
state. C/Fortran error messages need 256 bytes; `errflg=0` indicates success.
Existing representation-specific APIs remain available. Modal bounds and
chemistry/uptake derivatives retain their explicit interfaces. Native MAM4xx
helper support is queried separately from its modal representation label.

## Cloud and activation coupling

`cloudIce` has the same arguments and property/status definitions as the existing
bulk/modal calls in `CloudIceCoupling.hpp`, using the metadata population axis.
The common `activationFractions` returns separate number/mass fractions for
host-supplied supersaturations. Its output is `(cell,level,population,2*num_S)`.
`exaero_cloud_ice` and `exaero_activation_fractions` expose these operations to
C/Fortran with explicit `nf,np` dimensions. A null cloudborne pointer produces
`Absent`; a supplied empty reservoir produces `Empty`. Both reservoir snapshots
are host-owned. No aerosol transfer is inferred or applied.

## Physical population properties

`populationProperties` returns values and per-value availability with shape
`(cell,level,population,15)`, plus a separate `(cell,level,population)` empty/
populated status. Use `PopulationProperty` or matching `exaero_POP_*` constants.
The C/Fortran entry point is `exaero_population_properties`.

| Properties | Units / definition |
| --- | --- |
| Mass, Number | kg/m3 air; particles/m3 air |
| Dry/WetNumberMedianRadius | m; lognormal number-median radius |
| Dry/WetScavengingDiameter | m; existing representation's scavenging size |
| DryVolume | m3 dry aerosol / m3 air |
| DrySurfaceArea | m2 dry aerosol / m3 air; lognormal second moment |
| Kappa, GeometricStdDev | dimensionless; modal kappa is volume-weighted |
| Dry/WetDensity | kg/m3 of dry/wet particle volume, respectively |
| Lower/UpperDryRadius | m; configured physical radius bounds, when present |
| EffectiveDryRadius | m; configured bulk value or lognormal third/second-moment ratio |

Bulk scavenging diameter uses `dry_particle_diameter`; its number-median radius
uses `lognormal_dg/2`, so these can differ. Bulk kappa sizing retains its RH cap
of 0.99. Modal properties retain the configured existing water-uptake model
(default kappa sizing RH cap 0.98); the facade does not change those models.
The common surface area is **dry**; the older modal property surface reports
area at the wet radius.

Availability is `Unavailable`, `Configured`, or `Derived` (including values
read/derived from current state). An unavailable numerical slot is zero. Empty
populations have zero extensive quantities and retain configured dry/reference
geometry and width. Bulk configuration still defines kappa and wet reference
properties. Empty modal mixture kappa, density and wet geometry are unavailable.
Neither missing bounds nor an empty population triggers number correction.
