# Design: Structured (No-YAML) GOCART Package Configuration

- **Date:** 2026-09-10
- **Status:** Approved (user sign-off on sections 1–2, pending written-spec review)
- **Author:** Copilot + Barry (NOAA NWS OMD)
- **Relates to:** ADR-001 (compilation isolation), ADR-003 (curve mapping model),
  feature `001-geosmie-lut-attributes`, `docs/api-contracts.md` §4

## 1. Problem

`GocartPackage::initialize(const std::string& config_yaml)` parses YAML *text* (in
memory — never a file on disk; `YAML::Load`, not `YAML::LoadFile`). Callers that
already hold typed configuration in C++ (CCPP-ish drivers, embedded callers, tests)
must serialize it to YAML first. The typed curve/activation setters
(`setSpeciesCurveConfig`, `setAttributeActivation`) exist but are only usable
*after* `initialize()`, because `initialize()` owns the fixed precedence order
(reset → runtime file → config) and the one-shot solver-state upload.

**Goal:** let a C++ caller configure an entire GOCART package directly from
in-memory structs — everything the YAML can express — with zero YAML text.

## 2. Decisions (user-approved)

| # | Decision |
|---|----------|
| D1 | One-shot struct overload (not a builder, not setter+finalize). |
| D2 | Full coverage: species scalars, legacy RH lookups, activation, curve configs, emissions mappings. YAML becomes fully optional. |
| D3 | C++ first. C/Fortran bindings keep taking the YAML string; a flat C-ABI struct surface is explicitly out of scope. |
| D4 | Approach A — `ConfigReader` seam: all orchestration lives in one `initializeImpl(reader)`; YAML and struct paths are two readers over it. |

## 3. Public surface

### 3.1 New header `src/exaero/GocartConfig.hpp`

Public, Kokkos-free, mdspan-free (ADR-001 / SC-002 style, same constraints as
`AttributeQuery.hpp`). Presence of optional blocks is `std::optional` (C++23);
counts/lengths are always data (no fixed-size arrays cross the boundary).
`SpeciesCurveConfig` is **reused** from `AttributeQuery.hpp` — no parallel curve type.

```cpp
namespace exaero {

/// Legacy band-integrated RH optics lookup (Mode B), the four parallel lists
/// currently under has_optics_lookup/rh_bins/ext_lookup/ssa_lookup/asm_lookup.
struct GocartLegacyOpticsLookup {
    std::vector<double> rh;    // RH axis [fraction], size >= 2
    std::vector<double> ext;   // extinction, same size
    std::vector<double> ssa;   // single-scattering albedo, same size
    std::vector<double> asm_;  // asymmetry factor, same size ('asm' is a MS keyword)
};

/// One species: the eight scalars currently parsed from species[i], plus
/// optional lookup / MIE-table blocks.
struct GocartSpeciesConfig {
    std::string name;                  // required (YAML default "Species_i" not adopted)
    double dry_density            = 0; // [kg/m3]
    double molecular_weight       = 0; // [g/mol]
    double dry_particle_diameter  = 0; // [m]
    double hygroscopicity         = 0; // kappa
    double lognormal_sigma        = 0; // GSD
    double lognormal_dg           = 0; // GMD [m]
    double refractive_index_real  = 0; // n
    double refractive_index_imag  = 0; // k

    std::optional<GocartLegacyOpticsLookup> optics_lookup; // engaged => Mode B
    std::optional<SpeciesCurveConfig>       mie_table;     // engaged => curve binding
};

/// Top-level activation block (categories/species/data_file).
struct GocartActivationConfig {
    std::vector<AttributeCategory> categories;   // extra bits OR-ed onto the default mask
    std::vector<std::string> species;            // empty => all species
    std::string data_file;                       // empty => no runtime MIE file
};

/// One raw CECE emissions species and its split mappings (emissions_mapping[s]).
struct GocartEmissionsMappingConfig {
    std::string raw_name;                        // parity note: YAML parses but ignores
    struct Mapping {
        std::string target_species;
        double mass_split_fraction    = 0;
        bool   is_modal_mode          = false;
        double emitted_particle_diameter = 0;    // used iff is_modal_mode
        double lognormal_sigma        = 0;       // used iff is_modal_mode
    };
    std::vector<Mapping> mappings;
};

/// The whole package configuration — exactly what the YAML document expresses.
struct GocartConfig {
    std::vector<GocartSpeciesConfig> species;    // required non-empty (parity with YAML)
    std::optional<GocartActivationConfig> activation;
    std::vector<GocartEmissionsMappingConfig> emissions_mapping;
};

} // namespace exaero
```

### 3.2 One new interface method

`IAerosolPackage` gains an overload with a default that throws, following the same
additive pattern as `setSpeciesCurveConfig` (contract §2):

```cpp
virtual void initialize(const GocartConfig& config) {
    throw std::logic_error("EX-aero Error: this package does not support structured configuration");
}
```

`GocartPackage` overrides it. Existing packages (MAM4xx wrapper, test doubles)
compile and behave unchanged.

## 4. Internal seam (Approach A)

All real orchestration moves verbatim out of `initialize(const std::string&)` into
one private method in `GocartPackage.cpp`:

```
initialize(string)       -> YamlReader{YAML::Node}     -> initializeImpl(reader)
initialize(GocartConfig) -> StructReader{const GocartConfig&} -> initializeImpl(reader)
```

- `PackageConfigReader`, `YamlReader`, `StructReader` are **impl-private**
  (declared in `GocartPackage.cpp` or a `src/src_impl/gocart/` header — not public).
- `initializeImpl` keeps the exact current order and semantics:
  1. allocate `GocartSpeciesParams` + names,
  2. `MieTableStore::reset_to_baked_in()`,
  3. activation (`setAttributeActivation`, optional runtime file),
  4. curve configs (`setSpeciesCurveConfig`),
  5. micro / spectral / pmom pool wiring,
  6. emissions mappings,
  7. `create_solver_state` single H2D upload.
- **Precedence guarantee unchanged:** config > runtime file > baked-in, enforced by
  construction order, identical on both paths.
- `YamlReader` preserves the *exact* existing fail-fast messages ("missing 'species'
  field", "sets has_optics_lookup but omits a required lookup list", "lookup lists
  have mismatched lengths", "override names unknown attribute", ...).
- `StructReader` is direct field access. Struct-side fail-fast (same "EX-aero
  Error:" style, species name included): empty `species` vector; engaged
  `optics_lookup` with < 2 points or mismatched list lengths; engaged `mie_table`
  with non-empty `overrides` but empty `radius_nodes`; override `values.size()`
  != `radius_nodes.size()`; unknown attribute names (already validated inside
  `SpeciesCurveConfig::Override` resolution).

### 4.1 Reader interface (semantic, minimal)

```cpp
struct PackageConfigReader {
    virtual ~PackageConfigReader() = default;
    virtual int numSpecies() const = 0;
    virtual GocartSpeciesScalars speciesScalars(int i) const = 0; // 8 doubles + name
    virtual bool hasOpticsLookup(int i) const = 0;
    virtual LegacyLookup opticsLookup(int i) const = 0;           // 4 vectors
    virtual bool hasMieTable(int i) const = 0;
    virtual SpeciesCurveConfig mieTable(int i) const = 0;         // species_name filled
    // activation
    virtual bool hasActivation() const = 0;
    virtual int activationExtraMask() const = 0;                  // category bits
    virtual std::vector<std::string> activationSpecies() const = 0;
    virtual std::string activationDataFile() const = 0;
    // emissions
    virtual int numEmissionsMappings() const = 0;
    virtual EmissionsMappingEntry emissionsMapping(int s) const = 0;
};
```

(Exact member types are impl detail; the point is `initializeImpl` never touches
YAML nodes or `GocartConfig` fields directly.)

## 5. Behavior parity notes

- **Empty species:** YAML keeps throwing on a missing/empty `species` node; the
  struct path likewise requires non-empty `species` (parity chosen over
  convenience; a 0-species package is meaningless for GOCART).
- **`raw_name`** in emissions mappings is parsed but unused by today's YAML path;
  the struct keeps the field for parity and documents it as reserved.
- **`raw_cece_index`** remains sequential (`s`), matching current behavior; no new
  semantics introduced.
- **Species default naming:** YAML's `"Species_" + i` fallback is *not* adopted in
  the struct path — `name` is required there (typed callers should name things).

## 6. Testing

1. **Equivalence gate (the core test):** build one `GocartConfig` and the
   byte-equivalent YAML string; initialize two packages; assert bit-identical
   results for `queryAttribute` / `computeAttributes` / microphysics diagnostics /
   optics across a fixed probe set. Additions: curve-bound species, overrides,
   runtime-file activation, emissions mappings (all in one config).
2. **Fail-fast tests:** each struct-side violation in §4 throws with "EX-aero
   Error:" and the species name in the message.
3. **Regression:** all existing YAML tests and the four ctest suites stay green
   unchanged (backward compatibility proof).
4. **Unsupported package:** a package without the override throws `logic_error`
   from the default interface method.

## 7. Out of scope

- C / Fortran bindings for structured config (flat C-ABI struct surface deferred).
- Any change to curve semantics, resampling, interpolation, or the store.
- Builder / fluent APIs; setter-before-initialize lifecycles.
- Removal or deprecation of the YAML entry point (stays canonical).

## 8. Files touched

| File | Change |
|------|--------|
| `src/exaero/GocartConfig.hpp` | new public header |
| `src/exaero/IAerosolPackage.hpp` | additive virtual overload (needs `<stdexcept>`) |
| `src/src_impl/gocart/GocartPackage.hpp` | override declaration + private `initializeImpl` |
| `src/src_impl/gocart/GocartPackage.cpp` | reader seam, both entry points |
| `tests/test_main.cpp` | equivalence + fail-fast tests |
| `docs/api-contracts.md` | §4.x structured-config surface |
