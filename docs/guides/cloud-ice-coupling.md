# Cloud and ice coupling

Ex-Aero supplies diagnostic aerosol inputs to a host cloud and ice microphysics scheme. The host owns cloudborne aerosol buffers and decides activation, freezing, nucleation, and transfer tendencies. These helpers never change aerosol state. Chemical precursor formation remains upstream in MKPP.

The GOCART bulk and generic modal backends expose matching C++ functions, `bulkCloudIce` and `modalCloudIce`, with C/Fortran equivalents. Their ambient and optional cloudborne output slabs are `(cell,level,population,12)`: population means species for bulk and mode for modal. Each slot is:

| Slot | Quantity | Units |
| --- | --- | --- |
| 0, 1 | Dry mass, number | kg/m3, particles/m3 |
| 2 | Dry number-median radius | m |
| 3, 4 | Dry and wet scavenging radii | m |
| 5, 6, 7 | Lower, upper, effective dry radius | m |
| 8, 9 | Type and size weights for the queried ice category | fractions |
| 10 | Eligible number = number × type weight × size weight | particles/m3 |
| 11 | Configured heterogeneous-freezing eligibility | 0 or 1 |

The status slab is `(cell,level,population)`: `Absent=-2` for an omitted cloudborne buffer, `Empty=-1` for an empty population, `NotConfigured=0` for a populated population without an ice rule, and `Configured=1` when a rule is present. Missing physical lower/upper bounds are zero; `bulkReferenceRadii` separately reports whether they were configured. An omitted cloudborne input produces zero cloudborne outputs with `Absent` status. A bulk cloudborne buffer contains mass fields and gets the same fixed number conversion as ambient bulk mass. A modal cloudborne buffer uses the same packed mass and independent number field order as interstitial modal state.

Bulk number comes from configured `lognormal_dg`, `lognormal_sigma`, and density. Its number-median radius is `lognormal_dg/2`. Its scavenging radii use `dry_particle_diameter` and configured-kappa wet growth, so the two dry radii can differ. The lower/upper/effective radii are physical metadata and do not change fixed number or truncate the lognormal ice cutoff. None of these values depend on radiation or GEOSmie tables. Modal number is prognostic; its type weight is matching category dry mass divided by total mode dry mass.

Ice categories are `dust`, `sulfate`, `black_carbon`, and `stratospheric_sulfate`. Add `ice_categories` to each species and `ice_rules` to each bulk species or modal mode. For a bulk species:

```yaml
ice_categories: [dust]
ice_rules:
  - category: dust
    kind: above_dry_diameter
    value: 1.0e-7
    preexisting_ice: {kind: all}
    heterogeneous_freezing: true
```

The supported rule kinds are `none`, `all`, `fixed_fraction` (value in [0,1]), and `above_dry_diameter` (positive threshold in metres). A missing role or rule gives zero candidates; no species name supplies a role implicitly. A role does not imply heterogeneous-freezing eligibility; set the flag explicitly. Stratospheric sulfate is a separate query, so a host combining overlapping sulfate categories must avoid double counting. The examples are configuration mechanics, not a calibrated ice-nucleation scheme.

`bulkActivationFractions` and `modalActivationFractions` return interstitial-only liquid activation fractions in `(cell,level,population,2*nss)`; slots `2*j` and `2*j+1` are number and mass fractions at supersaturation `j`. Bulk preserves its current 0/1 critical-supersaturation decision. Modal uses lognormal number and third-moment mass fractions. Multiplying number by the number fraction and summing populations reproduces Ex-Aero's aggregate CCN spectrum. Already cloudborne aerosol is excluded from this activation query.

The C functions in `CloudIceCoupling_C.h` use zero-based indices, caller-owned column-major arrays, a nullable cloudborne pointer, and a 256-byte error buffer. The Fortran `exaero_interface` bindings use `c_loc(cloud_state)` for an optional cloudborne buffer or `c_null_ptr` when absent. For CAM-SIMA-style number density in particles/cm3, divide Ex-Aero's particles/m3 by 1e6; no air-density multiplication is needed. CAM-SIMA bulk's fixed `1/25` ice weight and unimplemented activation/transfer methods are not defaults here.
