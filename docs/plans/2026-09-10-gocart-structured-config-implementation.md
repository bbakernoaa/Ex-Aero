# Structured (No-YAML) GOCART Configuration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a C++ caller configure an entire GOCART package from in-memory structs (`GocartConfig`) with zero YAML text, while the existing YAML string entry point stays canonical and byte-for-byte behavior-compatible.

**Architecture:** Approach A ("ConfigReader seam", design spec `docs/plans/2026-09-10-gocart-structured-config-design.md`): all orchestration in `initialize()` moves verbatim into one private `initializeImpl(PackageConfigReader&)`. `YamlReader` (wraps `YAML::Node`) and `StructReader` (wraps `const GocartConfig&`) are impl-private readers over it. The public surface is one new Kokkos-free header plus one additive virtual overload on `IAerosolPackage` that defaults to throwing.

**Tech Stack:** C++23, CMake, yaml-cpp, Kokkos 5.1 (impl only — never in public headers), hand-rolled assert-style test runner (`tests/test_main.cpp`).

## Global Constraints

- C++23; every public class/struct/method has Doxygen `///` docs with `@brief` (`.github/instructions/cpp.instructions.md` §5).
- Public headers under `src/exaero/` MUST NOT include `Kokkos_Core.hpp` or define mdspan namespaces beyond the existing `IAerosolPackage.hpp` pattern (ADR-001, SC-002). `GocartConfig.hpp` includes only `<optional>`, `<string>`, `<vector>`, `<exaero/AttributeQuery.hpp>`.
- Fail fast, loudly: every validation error throws `std::runtime_error` with an `"EX-aero Error: "` prefix and the species name (never a silent fallback).
- New structured-config tests use explicit throw-based gates — the `gate(bool, const char*)` helper defined in Task 1 that throws `std::runtime_error("FATAL ERROR: ...")` — never bare `assert`, which compiles out under NDEBUG and would make Release runs assert nothing (T038 precedent). Existing assert-based tests are left as-is.
- Precedence is invariant on BOTH paths: `reset_to_baked_in → runtime file → config curves` (config > file > baked-in).
- Axis lengths are data, never literals (ADR-003): no hardcoded bin/RH/band counts in new code.
- Overloads: `initialize(const std::string&)` and `initialize(const GocartConfig&)` are unambiguous (a `GocartConfig` never converts to `std::string`).
- Build: `cmake --build build -j8` (Debug; `build-rel/` is the Release twin). Tests: `./build/tests/exaero_test_runner` then `ctest --test-dir build`.
- `specs/` and `docs/superpowers/` are git-ignored — NEVER force-add; tracked docs go in `docs/`.
- Do not commit `docs/guides/using-exaero-in-ccpp.md` (pre-existing uncommitted user edit).

---

### Task 1: Public config header + additive interface overload

**Files:**
- Create: `src/exaero/GocartConfig.hpp`
- Modify: `src/exaero/IAerosolPackage.hpp` (add include + virtual overload after line ~49 `virtual void initialize(const std::string& config_yaml) = 0;`)
- Test: `tests/test_main.cpp` (new `test_structured_config_unsupported_package()`)

**Interfaces:**
- Consumes: `exaero::SpeciesCurveConfig`, `exaero::AttributeCategory` from `src/exaero/AttributeQuery.hpp` (unchanged).
- Produces (later tasks rely on these exact names/types, all in `namespace exaero`):
  - `struct GocartLegacyOpticsLookup { std::vector<double> rh, ext, ssa, asm_; };`
  - `struct GocartSpeciesConfig { std::string name; double dry_density, molecular_weight, dry_particle_diameter, hygroscopicity, lognormal_sigma, lognormal_dg, refractive_index_real, refractive_index_imag; std::optional<GocartLegacyOpticsLookup> optics_lookup; std::optional<SpeciesCurveConfig> mie_table; };` (all doubles default 0)
  - `struct GocartActivationConfig { std::vector<AttributeCategory> categories; std::vector<std::string> species; std::string data_file; };`
  - `struct GocartEmissionsMappingConfig { std::string raw_name; struct Mapping { std::string target_species; double mass_split_fraction; bool is_modal_mode; double emitted_particle_diameter; double lognormal_sigma; }; std::vector<Mapping> mappings; };`
  - `struct GocartConfig { std::vector<GocartSpeciesConfig> species; std::optional<GocartActivationConfig> activation; std::vector<GocartEmissionsMappingConfig> emissions_mapping; };`
  - `virtual void IAerosolPackage::initialize(const GocartConfig&)` — default throws `std::logic_error("EX-aero Error: this package does not support structured configuration")`.

- [X] **Step 1: Write the failing test**

Add to `tests/test_main.cpp` (near the other Mie tests; register in `main()` after `test_mie_hot_path_overhead(...)`):

```cpp
// A package that does not override the structured overload must throw
// std::logic_error (additive-surface pattern, contract §2).
struct UnsupportedPkg : public exaero::IAerosolPackage {
    void initialize(const std::string&) override {}
    void executeMicrophysics(const exaero::EnvironmentalStateView&,
                             exaero::View3D<double>&, double) override {}
    void computeDerivedDiagnostics(const exaero::EnvironmentalStateView&,
                                   const exaero::View3D<const double>&,
                                   exaero::View3D<double>&) override {}
    void computeEmissions(const exaero::EnvironmentalStateView&,
                          const exaero::EmissionsInputView&,
                          exaero::View3D<double>&) override {}
    void computeOptics(const exaero::EnvironmentalStateView&,
                       const exaero::View3D<const double>&,
                       const exaero::View1D<const double>&,
                       exaero::View4D<double>&) override {}
    void computeCCN(const exaero::EnvironmentalStateView&,
                    const exaero::View3D<const double>&,
                    const exaero::View1D<const double>&,
                    exaero::View4D<double>&) override {}
    int getSpeciesIndex(const std::string&) const override { return -1; }
    std::string getSpeciesName(int) const override { return {}; }
};

// Explicit throw-based gate: bare assert() compiles out under NDEBUG and would
// make Release runs a no-op (T038 precedent). ALL structured-config tests use it.
static void gate(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(std::string("FATAL ERROR: structured-config test: ") + what);
}

void test_structured_config_unsupported_package() {
    UnsupportedPkg pkg;
    exaero::GocartConfig cfg;
    cfg.species.push_back(exaero::GocartSpeciesConfig{});
    bool threw = false;
    try { pkg.initialize(cfg); }
    catch (const std::logic_error& e) {
        threw = std::string(e.what()).find("does not support structured configuration") != std::string::npos;
    }
    gate(threw, "default structured initialize must throw logic_error with the contract message");
    std::cout << "Structured config unsupported-package guard: PASS" << std::endl;
}
```

Add `#include <exaero/GocartConfig.hpp>` at the top of `test_main.cpp` (with the other exaero includes).

- [X] **Step 2: Run test to verify it fails to compile**

Run: `cmake --build build -j8 2>&1 | tail -20`
Expected: FAIL — `exaero/GocartConfig.hpp: No such file or directory`.

- [X] **Step 3: Write `src/exaero/GocartConfig.hpp`**

Full content (Doxygen on every type/member; mirrors design spec §3.1 exactly):

```cpp
#pragma once
// Structured (no-YAML) configuration surface for the GOCART package.
//
// MUST NOT include <Kokkos_Core.hpp> (ADR-001): public, zero-dependency header.
// Every optional block is std::optional; all array lengths are data — no
// fixed-size shapes cross the public boundary. Species curve bindings reuse
// SpeciesCurveConfig from AttributeQuery.hpp (ADR-003); there is no parallel
// curve type.
#include <exaero/AttributeQuery.hpp>
#include <optional>
#include <string>
#include <vector>

namespace exaero {

/// @brief Legacy band-integrated RH optics lookup (Mode B): four parallel
/// lists indexed by the rh axis. Engaged optics_lookup replaces the ADT
/// analytical solver for that species; rh must hold >= 2 strictly ordered points.
struct GocartLegacyOpticsLookup {
    std::vector<double> rh;   ///< RH axis [fraction 0-1], size >= 2
    std::vector<double> ext;  ///< Extinction lookup, same length as rh
    std::vector<double> ssa;  ///< Single-scattering albedo lookup, same length as rh.
    std::vector<double> asm_; ///< Asymmetry-factor lookup, same length as rh ('asm' is reserved on MSVC).
};

/// @brief One aerosol species: the eight physical scalars plus optional
/// legacy-lookup and MIE-curve blocks.
struct GocartSpeciesConfig {
    std::string name;                    ///< Species label (required; e.g. "DU").
    double dry_density = 0;              ///< [kg m^-3]
    double molecular_weight = 0;         ///< [g mol^-1]
    double dry_particle_diameter = 0;    ///< [m]
    double hygroscopicity = 0;           ///< kappa (Kohler), dimensionless
    double lognormal_sigma = 0;          ///< Geometric standard deviation
    double lognormal_dg = 0;             ///< Geometric mean diameter [m]
    double refractive_index_real = 0;    ///< n
    double refractive_index_imag = 0;    ///< k

    std::optional<GocartLegacyOpticsLookup> optics_lookup; ///< Engaged => Mode B RH lookup.
    std::optional<SpeciesCurveConfig> mie_table;           ///< Engaged => curve binding (ADR-003 R9).
};

/// @brief Top-level attribute activation block (FR-010/FR-012 parity with YAML
/// `activation:`). Absent activation engages the defaults: microphysical +
/// spectral-optical categories, all species, no runtime file.
struct GocartActivationConfig {
    std::vector<AttributeCategory> categories; ///< Extra category bits OR-ed onto the default mask.
    std::vector<std::string> species;          ///< Activated species; empty => all.
    std::string data_file;                     ///< Optional runtime MIE table file; empty => none.
};

/// @brief Emissions mapping for one raw CECE species (SPEC-EMISSIONS-002 parity
/// with YAML `emissions_mapping[s]`).
struct GocartEmissionsMappingConfig {
    std::string raw_name; ///< Reserved for diagnostics; unused by the mapping engine (YAML parity).

    /// @brief One split of a raw species into a package target species.
    struct Mapping {
        std::string target_species;             ///< Must match a species name in the config.
        double mass_split_fraction = 0;         ///< Fraction of raw mass routed to the target.
        bool is_modal_mode = false;             ///< True => number emission uses a modal mode.
        double emitted_particle_diameter = 0;   ///< [m] used iff is_modal_mode.
        double lognormal_sigma = 0;             ///< GSD used iff is_modal_mode.
    };
    std::vector<Mapping> mappings; ///< Sequential raw index == position in the parent vector.
};

/// @brief Complete GOCART package configuration — exactly what the YAML
/// document expresses, passed directly in memory. `initialize(const
/// GocartConfig&)` applies it with the same precedence as the YAML path:
/// config curves > runtime file > baked-in tables.
struct GocartConfig {
    std::vector<GocartSpeciesConfig> species; ///< Required, non-empty.
    std::optional<GocartActivationConfig> activation; ///< Absent => defaults (micro+spectral, all species).
    std::vector<GocartEmissionsMappingConfig> emissions_mapping; ///< Optional; may be empty.
};

} // namespace exaero
```

- [X] **Step 4: Add the overload to `IAerosolPackage.hpp`**

After `virtual void initialize(const std::string& config_yaml) = 0;` add:

```cpp
        // Structured (no-YAML) initialization: same semantics and precedence as the
        // YAML string path — config curves > runtime file > baked-in (ADR-003).
        // Default throws so packages without structured support (MAM4xx wrapper,
        // test doubles) compile and behave unchanged (contract §2 additive pattern).
        virtual void initialize(const GocartConfig& config) {
            (void)config;
            throw std::logic_error(
                "EX-aero Error: this package does not support structured configuration");
        }
```

Also add at the top of the file: `#include <exaero/GocartConfig.hpp>` and `#include <stdexcept>`.

- [X] **Step 5: Verify name-hiding safety**

Run: `grep -rn "public IAerosolPackage\|: IAerosolPackage" src tests | grep -v GocartPackage`
For each concrete subclass that calls `initialize(...)` on the *concrete* type anywhere, add `using IAerosolPackage::initialize;` in its public section (calling through an `IAerosolPackage&` always works without it). Expected current hits: test doubles only — inspect and add the using-declaration if any direct concrete call sites exist.

- [X] **Step 6: Build and run the full suite**

Run: `cmake --build build -j8 && ./build/tests/exaero_test_runner 2>&1 | tail -8`
Expected: `Structured config unsupported-package guard: PASS` and all prior tests PASS, exit 0.

- [X] **Step 7: Commit**

```sh
git add src/exaero/GocartConfig.hpp src/exaero/IAerosolPackage.hpp tests/test_main.cpp
git commit -m "feat(api): add GocartConfig structured-configuration surface"
```

---

### Task 2: ConfigReader seam — pure refactor of the YAML path

**Files:**
- Modify: `src/src_impl/gocart/GocartPackage.cpp` (the `initialize(const std::string&)` body, ~lines 108–365)
- Test: no new tests — the ENTIRE existing suite is the regression gate (this task must be behavior-identical).

**Interfaces:**
- Consumes: nothing new.
- Produces: `exaero::PackageConfigReader` (impl-private, anonymous namespace of `GocartPackage.cpp`) with the exact virtuals below; `GocartPackage::initializeImpl(PackageConfigReader&)` (private member, declared in `GocartPackage.hpp`); `exaero::(anon)::YamlReader` implementing the reader over `YAML::Node`. Task 3 consumes the reader interface verbatim.

```cpp
// Semantic view over a package configuration (YAML node or GocartConfig struct).
// initializeImpl() consumes ONLY this interface — it never touches YAML nodes or
// GocartConfig fields directly, so both entry points share one code path and the
// precedence order (reset -> file -> config) cannot drift between them.
struct PackageConfigReader {
    virtual ~PackageConfigReader() = default;
    virtual int numSpecies() const = 0;
    virtual std::string speciesName(int i) const = 0;
    virtual GocartSpeciesParams speciesScalars(int i) const = 0; // 8 doubles + has_optics_lookup; offsets filled by impl
    virtual bool hasOpticsLookup(int i) const = 0;
    virtual GocartLegacyOpticsLookup opticsLookup(int i) const = 0;      // validated: n>=2, equal lengths
    virtual bool hasMieTable(int i) const = 0;
    virtual SpeciesCurveConfig mieTable(int i) const = 0;                // species_name pre-filled
    virtual bool hasActivation() const = 0;
    virtual int activationExtraMask() const = 0;                          // category bits OR-ed on default
    virtual std::vector<std::string> activationSpecies() const = 0;
    virtual std::string activationDataFile() const = 0;
    virtual int numEmissionsMappings() const = 0;
    virtual GocartEmissionsMappingConfig emissionsMapping(int s) const = 0; // raw_name ignored downstream (parity)
};
```

- [X] **Step 1: Declare `initializeImpl` in `GocartPackage.hpp`**

In the private section (after `GocartSolverState* solver_state_ = nullptr;`):

```cpp
        // Shared orchestration for both entry points (design 2026-09-10, Approach A):
        // species params -> store reset -> activation -> curve configs -> pool wiring
        // -> emissions -> solver-state upload. Reads config ONLY through the reader.
        void initializeImpl(class PackageConfigReader& reader);
```

- [X] **Step 2: Move the body**

In `GocartPackage.cpp`: create the anonymous-namespace `PackageConfigReader` interface above (place it next to `attribute_key`), then move the ENTIRE current body of `initialize(const std::string&)` into `GocartPackage::initializeImpl(PackageConfigReader& reader)` replacing every `config[...]`/`s[...]` access with the reader calls:

- `config["species"]` presence check + `species_node.size()` → `reader.numSpecies()` (reader constructor performs the `missing 'species' field` throw — see Step 3).
- species loop: `s["dry_density"].as<double>()` etc. → `reader.speciesScalars(i)` (return a `GocartSpeciesParams` with the 8 scalars + `has_optics_lookup`; leave `curve_offset`/`n_radius`/`n_rh`/`n_lambda` to `initializeImpl` exactly as today).
- legacy lookup block: the four `*_lookup`/`rh_bins` reads + their three fail-fast checks (`omits a required lookup list`, `needs at least 2 points`, `mismatched lengths`) → `reader.hasOpticsLookup(i)` / `reader.opticsLookup(i)`; the pool append order `{rh, ext, ssa, asm}` and `curve_slot` semantics are unchanged in `initializeImpl`.
- `mie_table` loop → `reader.hasMieTable(i)` / `reader.mieTable(i)`; the attribute-name validation (`override names unknown attribute`) stays where it is — it happens inside `apply_curve_config`/`attribute_key`, not in the reader.
- activation block → `reader.hasActivation()`, `reader.activationExtraMask()`, `reader.activationSpecies()`, `reader.activationDataFile()`; default mask `Microphysical | SpectralOptical` stays in `initializeImpl`.
- emissions block → `reader.numEmissionsMappings()` / `reader.emissionsMapping(s)`; the `is_active`/`raw_cece_index = s`/fraction/modal logic is copied unchanged.
- micro/spectral/pmom pool wiring loops, `store.reset_to_baked_in()`, `setAttributeActivation(...)`, `setSpeciesCurveConfig(...)`, and the `create_solver_state` upload: byte-identical, still in the same order.

Then:

```cpp
    void GocartPackage::initialize(const std::string& config_yaml) {
        YamlReader reader{YAML::Load(config_yaml)};
        initializeImpl(reader);
    }
```

- [X] **Step 3: Implement `YamlReader`**

In the same anonymous namespace. Constructor reproduces today's document-level fail-fast:

```cpp
    class YamlReader final : public PackageConfigReader {
    public:
        explicit YamlReader(YAML::Node config) : config_(std::move(config)) {
            if (!config_["species"]) {
                throw std::runtime_error("GOCART YAML config missing 'species' field");
            }
            species_ = config_["species"];
        }
        // ... one method per reader virtual, each a direct port of the deleted
        // inline YAML access, preserving EXACT error strings, e.g.:
        GocartLegacyOpticsLookup opticsLookup(int i) const override {
            auto s = species_[i];
            if (!(s["rh_bins"] && s["ext_lookup"] && s["ssa_lookup"] && s["asm_lookup"])) {
                throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                    "' sets has_optics_lookup but omits a required lookup list");
            }
            // n_rh >= 2 and equal-length checks: same messages as today.
            ...
        }
    private:
        YAML::Node config_, species_;
    };
```

`speciesScalars(i)` must throw on missing `dry_density`/... exactly like `.as<double>()` on an absent node does today (yaml-cpp `TypedKeyNotFound` propagates — do not catch it; parity).

- [X] **Step 4: Build and run the FULL regression suite**

Run: `cmake --build build -j8 && ./build/tests/exaero_test_runner 2>&1 | tail -30 && ctest --test-dir build 2>&1 | tail -8`
Expected: every existing test PASS (including all `test_gocart_yaml_*`, `test_mie_*`, `ExaeroMieDeterminism`), 4/4 ctest. Any failure = the refactor changed behavior — fix by restoring the original semantics, do not adjust the test.

- [X] **Step 5: Commit**

```sh
git add src/src_impl/gocart/GocartPackage.hpp src/src_impl/gocart/GocartPackage.cpp
git commit -m "refactor(gocart): route YAML initialize through PackageConfigReader seam"
```

---

### Task 3: StructReader + `GocartPackage::initialize(const GocartConfig&)` — species parity

**Files:**
- Modify: `src/src_impl/gocart/GocartPackage.cpp` (add `StructReader`, override body)
- Modify: `src/src_impl/gocart/GocartPackage.hpp` (add `void initialize(const GocartConfig& config) override;`)
- Test: `tests/test_main.cpp` — `test_structured_config_species_parity()`

**Interfaces:**
- Consumes: `PackageConfigReader` (Task 2), `GocartConfig` types (Task 1).
- Produces: `GocartPackage::initialize(const GocartConfig&)` — fully functional for species scalars + legacy lookups; activation/curve/emissions blocks wired in Task 4 (the StructReader already delegates them; only validation tests are added there).

- [X] **Step 1: Write the failing parity test**

`kMieSixYaml` (test_main.cpp:312) is the fixture to mirror. Add:

```cpp
// Build the GocartConfig equivalent of kMieSixYaml (six species, no blocks).
static exaero::GocartConfig kMieSixConfig() {
    auto sp = [](const char* n, double dens, double mw, double dpg, double kap,
                 double sig, double dg, double nr, double ni) {
        exaero::GocartSpeciesConfig c;
        c.name = n; c.dry_density = dens; c.molecular_weight = mw;
        c.dry_particle_diameter = dpg; c.hygroscopicity = kap;
        c.lognormal_sigma = sig; c.lognormal_dg = dg;
        c.refractive_index_real = nr; c.refractive_index_imag = ni;
        return c;
    };
    exaero::GocartConfig cfg;
    cfg.species = {
        sp("DU", 2600.0, 100.0, 2.0e-6, 0.1, 1.5, 1.0e-6, 1.55, 0.002),
        sp("SS", 1800.0, 98.0, 0.2e-6, 0.5, 2.0, 0.15e-6, 1.43, 1.0e-8),
        sp("SU", 1800.0, 98.0, 0.2e-6, 0.5, 2.0, 0.15e-6, 1.43, 1.0e-8),
        sp("BC", 1800.0, 12.0, 0.1e-6, 0.0, 1.8, 0.05e-6, 1.85, 0.75),
        sp("OC", 1300.0, 150.0, 0.1e-6, 0.1, 1.8, 0.05e-6, 1.55, 0.0),
        sp("NI", 2150.0, 85.0, 0.3e-6, 0.6, 1.6, 0.1e-6, 1.52, 0.01),
    };
    return cfg;
}

void test_structured_config_species_parity() {
    exaero::GocartPackage pkg_yaml, pkg_cfg;
    pkg_yaml.initialize(kMieSixYaml);
    pkg_cfg.initialize(kMieSixConfig());

    gate(pkg_cfg.get_num_species() == 6, "six species configured");
    for (int i = 0; i < 6; ++i) {
        auto a = pkg_yaml.get_species_params(i);
        auto b = pkg_cfg.get_species_params(i);
        gate(pkg_cfg.getSpeciesName(i) == pkg_yaml.getSpeciesName(i), "species name parity");
        // Bit-identical scalars and pool wiring (offsets/extents included).
        gate(b.dry_density == a.dry_density && b.molecular_weight == a.molecular_weight,
             "density/molecular-weight parity");
        gate(b.dry_particle_diameter == a.dry_particle_diameter, "dpg parity");
        gate(b.hygroscopicity == a.hygroscopicity && b.lognormal_sigma == a.lognormal_sigma,
             "kappa/sigma parity");
        gate(b.lognormal_dg == a.lognormal_dg, "dg parity");
        gate(b.refractive_index_real == a.refractive_index_real, "n parity");
        gate(b.refractive_index_imag == a.refractive_index_imag, "k parity");
        gate(b.has_optics_lookup == a.has_optics_lookup, "lookup flag parity");
        gate(b.curve_offset == a.curve_offset && b.n_rh == a.n_rh, "legacy curve wiring parity");
        gate(b.micro_offset == a.micro_offset && b.spec_offset == a.spec_offset,
             "micro/spectral wiring parity");
        gate(b.pmom_offset == a.pmom_offset, "pmom wiring parity");
    }
    std::cout << "Structured config species parity: PASS" << std::endl;
}
```

Register in `main()` after `test_structured_config_unsupported_package()`.

- [X] **Step 2: Run to verify it fails**

Run: `cmake --build build -j8 2>&1 | tail -5`
Expected: compile error — `GocartPackage` has no `initialize(const GocartConfig&)` overload visible / base default throws at runtime once it compiles.

- [X] **Step 3: Implement `StructReader` + override**

In `GocartPackage.cpp` anonymous namespace:

```cpp
    // Reader over an in-memory GocartConfig (design 2026-09-10 §4). Near-trivial by
    // design: the caller already supplied typed data. Validation duplicates the YAML
    // fail-fast guarantees with the same "EX-aero Error:" style (never silent).
    class StructReader final : public PackageConfigReader {
    public:
        explicit StructReader(const GocartConfig& config) : config_(config) {
            if (config_.species.empty()) {
                throw std::runtime_error("EX-aero Error: structured config has no species");
            }
            for (const auto& s : config_.species) {
                if (s.name.empty()) {
                    throw std::runtime_error("EX-aero Error: structured config species at index has empty name");
                }
            }
        }
        int numSpecies() const override { return static_cast<int>(config_.species.size()); }
        std::string speciesName(int i) const override { return config_.species[i].name; }
        GocartSpeciesParams speciesScalars(int i) const override {
            const auto& s = config_.species[i];
            GocartSpeciesParams p{s.dry_density, s.molecular_weight, s.dry_particle_diameter,
                                  s.hygroscopicity, s.lognormal_sigma, s.lognormal_dg,
                                  s.refractive_index_real, s.refractive_index_imag};
            p.has_optics_lookup = s.optics_lookup.has_value();
            return p;
        }
        bool hasOpticsLookup(int i) const override {
            return config_.species[i].optics_lookup.has_value();
        }
        GocartLegacyOpticsLookup opticsLookup(int i) const override {
            const auto& l = *config_.species[i].optics_lookup; // caller checked engaged
            const int n = static_cast<int>(l.rh.size());
            if (n < 2) {
                throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                    "' rh_bins needs at least 2 points to interpolate");
            }
            if (static_cast<int>(l.ext.size()) != n || static_cast<int>(l.ssa.size()) != n ||
                static_cast<int>(l.asm_.size()) != n) {
                throw std::runtime_error("EX-aero Error: species '" + speciesName(i) +
                    "' lookup lists have mismatched lengths (fail fast, no silent fallback)");
            }
            return l;
        }
        bool hasMieTable(int i) const override { return config_.species[i].mie_table.has_value(); }
        SpeciesCurveConfig mieTable(int i) const override {
            SpeciesCurveConfig cc = *config_.species[i].mie_table; // caller checked engaged
            cc.species_name = config_.species[i].name;            // config name wins (R9)
            return cc;
        }
        bool hasActivation() const override { return config_.activation.has_value(); }
        int activationExtraMask() const override {
            int mask = 0;
            for (auto c : config_.activation->categories) mask |= attribute_category_bit(c);
            return mask;
        }
        std::vector<std::string> activationSpecies() const override {
            return config_.activation ? config_.activation->species : std::vector<std::string>{};
        }
        std::string activationDataFile() const override {
            return config_.activation ? config_.activation->data_file : std::string();
        }
        int numEmissionsMappings() const override {
            return static_cast<int>(config_.emissions_mapping.size());
        }
        GocartEmissionsMappingConfig emissionsMapping(int s) const override {
            return config_.emissions_mapping[s];
        }
    private:
        const GocartConfig& config_;
    };
```

In `GocartPackage.hpp` public section, next to the string overload:

```cpp
        void initialize(const GocartConfig& config) override;
```

In `GocartPackage.cpp`:

```cpp
    void GocartPackage::initialize(const GocartConfig& config) {
        StructReader reader{config};
        initializeImpl(reader);
    }
```

- [X] **Step 4: Run tests**

Run: `cmake --build build -j8 && ./build/tests/exaero_test_runner 2>&1 | tail -6`
Expected: `Structured config species parity: PASS`, all prior tests PASS.

- [X] **Step 5: Commit**

```sh
git add src/src_impl/gocart/GocartPackage.hpp src/src_impl/gocart/GocartPackage.cpp tests/test_main.cpp
git commit -m "feat(gocart): structured initialize via StructReader (species parity)"
```

---

### Task 4: Full-feature equivalence + struct fail-fast tests

**Files:**
- Modify: `tests/test_main.cpp` (two new tests; one helper)
- Possibly modify: `src/src_impl/gocart/GocartPackage.cpp` ONLY if a test exposes a real seam bug (fix the seam, not the test).

**Interfaces:**
- Consumes: everything above; `tools/make_test_fixtures.py`'s `tests/data/mie_dust7.txt` (C12 fixture) for the runtime-file leg; `mie_file_yaml()` helper (exists in test_main.cpp).
- Produces: final proof that the struct path is complete.

- [X] **Step 1: Full-feature equivalence test**

`test_structured_config_full_equivalence()`: build BOTH a YAML string and a `GocartConfig` that exercise every block at once — one species with `optics_lookup` (reuse the 4-point lookup from `test_gocart_arbitrary_lookup_length`), one species with `mie_table` (source `DU`, one inserted `radius_nodes` + one `bext` override, mirroring the C11/C13 blocks in `test_mie_config_curves`), `activation` (categories `[polarized]`, species `["DU"]`... note: YAML `categories` are strings, struct uses `AttributeCategory` — same bits), and an `emissions_mapping` with one modal mapping (mirror test_main.cpp:1503). Initialize two `GocartPackage`s and gate every check with Task 1's `gate()` helper (throw-based, holds under NDEBUG):

1. All `get_species_params(i)` fields bit-identical (as Task 3, plus `emissions_mapping` subfields).
2. `queryAttribute` over a fixed probe matrix (species × category × attribute_index × {rh=0, 0.495, 0.99} × {1e-6, 2e-6, 5e-6} m) returns identical status codes and bit-identical values.
3. `computeAttributes` on a 2×3 grid with `wavelengths {3.0, 5.0}`: identical values and statuses.
4. `momentCounts` identical for every species.

Use the same `store`-ordering care as `test_mie_config_curves` (fresh package per case; the MieTableStore is a singleton reset by `initialize`).

- [X] **Step 2: Runtime-file leg equivalence**

Extend the full-equivalence test (or a sibling) with `activation.data_file = "tests/data/mie_dust7.txt"`-equivalent content: the YAML variant uses `mie_file_yaml(...)` exactly as `test_mie_config_curves` does; the struct variant sets `data_file` to the same path. Assert DUST7-sourced queries (`qext` at rh=0.495, band 2.0 → the known `4.31`, `AvailableFile`) match on both paths.

- [X] **Step 3: Struct fail-fast tests**

`test_structured_config_failfast()` — each case must throw `std::runtime_error` whose message starts `"EX-aero Error:"` and (where applicable) contains the species name:

```cpp
auto expect_throw = [](const exaero::GocartConfig& cfg, const char* what) {
    exaero::GocartPackage pkg;
    bool threw = false;
    try { pkg.initialize(cfg); } catch (const std::runtime_error& e) {
        threw = std::string(e.what()).rfind("EX-aero Error:", 0) == 0;
    }
    gate(threw, what);
};
```

Cases: (a) empty `species`; (b) species with `optics_lookup` of length 1; (c) lookup lists of mismatched lengths; (d) `mie_table` engaged with `overrides` but empty `radius_nodes` (add the check to `StructReader::mieTable` mirroring the YAML guarantee — override values must be resampleable; message: `"...overrides require radius_nodes"`); (e) override `values.size()` != `radius_nodes.size()`.

- [X] **Step 4: Run everything, Debug + Release**

Run: `cmake --build build -j8 && ./build/tests/exaero_test_runner 2>&1 | tail -6 && ctest --test-dir build 2>&1 | tail -5`
Then: `cmake --build build-rel -j8 && ./build-rel/tests/exaero_test_runner 2>&1 | tail -6 && ctest --test-dir build-rel 2>&1 | tail -5`
Expected: all PASS, 4/4 ctest in both configs (Release catches NDEBUG-gated issues — see the intent(inout) lesson).

- [X] **Step 5: Commit**

```sh
git add tests/test_main.cpp src/src_impl/gocart/GocartPackage.cpp
git commit -m "test(gocart): struct/YAML full equivalence + structured fail-fast coverage"
```

---

### Task 5: Documentation + final verification

**Files:**
- Modify: `docs/api-contracts.md` (§4 attribute surface — add §4.x structured config)
- Modify: `docs/codebase-guide.md` (layout + namespaces: `GocartConfig.hpp`)
- Modify: `docs/guides/getting-started.md` ONLY if it shows YAML init (add a short "or pass a `GocartConfig`" snippet)

**Interfaces:**
- Consumes: final API.
- Produces: contract documentation for the new surface.

- [X] **Step 1: Document in `docs/api-contracts.md`**

New subsection under §4: table of the five config structs (fields, units, optionality), the `initialize(const GocartConfig&)` signature, default-throws contract for non-supporting packages, precedence statement (config > runtime file > baked-in, identical on both paths), and the fail-fast list from design §4.

- [X] **Step 2: Update `docs/codebase-guide.md`**

Add `GocartConfig.hpp` to the public-header layout listing and the `exaero::` namespace table (one line: "Structured configuration descriptors — Kokkos-free").

- [X] **Step 3: Final full verification**

Run: `cmake --build build -j8 && ./build/tests/exaero_test_runner && ctest --test-dir build && cmake --build build-rel -j8 && ./build-rel/tests/exaero_test_runner && ctest --test-dir build-rel`
Expected: exit 0 everywhere; 4/4 ctest each config.

- [X] **Step 4: Commit**

```sh
git add docs/api-contracts.md docs/codebase-guide.md
git commit -m "docs: specify structured (no-YAML) GOCART configuration contract"
```

---

## Self-review checklist (run after writing, before executing)

- Spec coverage: design §3.1 structs → Task 1; §3.2 overload → Task 1; §4 seam → Task 2; §4 struct path + validation → Task 3/4; §6 tests 1–4 → Tasks 1/3/4; §8 file list → Tasks 1–5. No gaps.
- No placeholder language ("TBD", "similar to Task N") in code steps — every code block is complete or a direct port instruction with the source line range.
- Type names consistent across tasks: `GocartLegacyOpticsLookup::asm_`, `PackageConfigReader` virtuals identical in Task 2 declaration and Task 3 implementation.
- Gate discipline: every new structured-config test check goes through `gate()` (throw-based), never bare `assert` (user decision 2026-09-10, reviewer plan-mandated finding).
