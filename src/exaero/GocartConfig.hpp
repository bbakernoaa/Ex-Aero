#pragma once
// Structured (no-YAML) configuration surface for the GOCART package.
//
// MUST NOT include <Kokkos_Core.hpp> public, zero-dependency header.
// Every optional block is std::optional; all array lengths are data — no
// fixed-size shapes cross the public boundary. Species curve bindings reuse
// SpeciesCurveConfig from AttributeQuery.hpp ; there is no parallel
// curve type.
#include <exaero/AttributeQuery.hpp>
#include <optional>
#include <string>
#include <vector>

namespace exaero {

/// @brief Legacy band-integrated RH optics lookup (Mode B): four parallel
/// lists indexed by the rh axis. Engaged optics_lookup replaces the ADT
/// analytical solver for that species; rh must hold >= 2 strictly ordered
/// points.
struct GocartLegacyOpticsLookup {
  std::vector<double> rh;  ///< RH axis [fraction 0-1], size >= 2
  std::vector<double> ext; ///< Extinction lookup, same length as rh
  std::vector<double>
      ssa; ///< Single-scattering albedo lookup, same length as rh.
  std::vector<double> asm_; ///< Asymmetry-factor lookup, same length as rh
                            ///< ('asm' is reserved on MSVC).
};

/// @brief One aerosol species: the eight physical scalars plus optional
/// legacy-lookup and MIE-curve blocks.
struct GocartSpeciesConfig {
  std::string name;                 ///< Species label (required; e.g. "DU").
  double dry_density = 0;           ///< [kg m^-3]
  double molecular_weight = 0;      ///< [g mol^-1]
  double dry_particle_diameter = 0; ///< [m]
  double hygroscopicity = 0;        ///< kappa (Kohler), dimensionless
  double lognormal_sigma = 0;       ///< Geometric standard deviation
  double lognormal_dg = 0;          ///< Geometric mean diameter [m]
  double refractive_index_real = 0; ///< n
  double refractive_index_imag = 0; ///< k

  std::optional<GocartLegacyOpticsLookup>
      optics_lookup;                           ///< Engaged => Mode B RH lookup.
  std::optional<SpeciesCurveConfig> mie_table; ///< Engaged => curve binding.
};

/// @brief Top-level attribute activation block (parity with YAML
/// `activation:`). Absent activation engages the defaults: microphysical +
/// spectral-optical categories, all species, no runtime file.
struct GocartActivationConfig {
  std::vector<AttributeCategory>
      categories; ///< Extra category bits OR-ed onto the default mask.
  std::vector<std::string> species; ///< Activated species; empty => all.
  std::string data_file; ///< Optional runtime MIE table file; empty => none.
};

/// @brief Emissions mapping for one raw CECE species (parity
/// with YAML `emissions_mapping[s]`).
struct GocartEmissionsMappingConfig {
  std::string raw_name; ///< Reserved for diagnostics; unused by the mapping
                        ///< engine (YAML parity).

  /// @brief One split of a raw species into a package target species.
  struct Mapping {
    std::string target_species; ///< Must match a species name in the config.
    double mass_split_fraction =
        0;                      ///< Fraction of raw mass routed to the target.
    bool is_modal_mode = false; ///< True => number emission uses a modal mode.
    double emitted_particle_diameter = 0; ///< [m] used iff is_modal_mode.
    double lognormal_sigma = 0;           ///< GSD used iff is_modal_mode.
  };
  std::vector<Mapping>
      mappings; ///< Sequential raw index == position in the parent vector.
};

/// @brief Complete GOCART package configuration — exactly what the YAML
/// document expresses, passed directly in memory. `initialize(const
/// GocartConfig&)` applies it with the same precedence as the YAML path:
/// config curves > runtime file > baked-in tables.
struct GocartConfig {
  std::vector<GocartSpeciesConfig> species; ///< Required, non-empty.
  std::optional<GocartActivationConfig>
      activation; ///< Absent => defaults (micro+spectral, all species).
  std::vector<GocartEmissionsMappingConfig>
      emissions_mapping; ///< Optional; may be empty.
};

} // namespace exaero
