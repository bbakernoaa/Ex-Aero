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
