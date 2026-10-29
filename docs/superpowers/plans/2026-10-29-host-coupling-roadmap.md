# Common Aerosol Host-Coupling Roadmap

**Goal:** Give hosts one representation-aware interface for applicable bulk and
generic modal helpers, using the package's representation and explicit metadata.

**Design:** [Common host coupling](../specs/2026-10-29-host-coupling-design.md).

## Delivery order

| Slice | Deliverable | Detailed plan |
| --- | --- | --- |
| 1 | Population/constituent mapping, operation capabilities, constituent-specific modal removal, common rate and velocity adapters | [Metadata and removal](2026-10-29-representation-aware-removal.md) |
| 2 | Common cloud/ice and activation dispatch using the same population axis | [Cloud/ice](2026-10-29-host-cloud-ice-api.md) |
| 3 | Common physical population properties with explicit units, size definitions and availability | [Properties](2026-10-29-host-population-properties-api.md) |
| 4 | Common resuspension routing and deposition rebinning using explicit field/species selectors | [Deposition utilities](2026-10-29-host-deposition-utilities-api.md) |

Execute slice 1 first. Slices 2–4 consume its metadata and capability contract;
each delivers complete C++, C and Fortran calls, tests and a host usage guide.
Use native inline execution for continuity unless the user requests otherwise.
The execution skill supplies a final whole-slice review. No implementation,
branch publication or integration is included in this documentation revision.

## Shared decisions

- The initialized package selects bulk/modal behavior; no separate selector can
  disagree with it. Capability checks also inspect the actual extension support.
- Hosts query populations and field mappings once per initialized configuration,
  then size arrays from counts; field names are labels, not encoded APIs.
- Bulk number is diagnosed; modal number is one field per mode. A constituent
  never gains its own modal particle-number prognostic.
- Supported operations, empty populations and unavailable properties have
  distinct meanings and explicit codes.
- Existing bulk/modal helpers and common optics/diagnostic APIs remain available;
  the facade reuses their calculations and preserves their contracts.
- Modal bounds and chemistry/uptake derivatives stay explicit extensions. A
  capability query reports their applicability without manufacturing bulk analogs.
- Each numerical call is synchronous, uses column-major host memory, and either
  publishes complete outputs or preserves them on error.
- Host/MKPP ownership of chemistry, removal physics and accepted updates continues.

## Completion criteria

For each slice, demonstrate two bulk species/bins and two modal populations
with at least two configured materials, multiple cells/levels, nontrivial state
field mapping and errors that preserve sentinels. Tests compare facade outputs
with existing kernels where applicable, plus independent analytical cases for
new behavior. C and Fortran must agree with C++ units and indices.

Run the existing generic and optional MAM4xx configured CTest profiles after
each completed slice; optional builds must retain their existing functionality
and report unsupported facade helpers explicitly. Add public headers to the
standalone compile gate. Do not infer new scientific-scheme parity from facade
tests or include optical dataset generation in these slices.
