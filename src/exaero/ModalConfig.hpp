#pragma once
/// @file ModalConfig.hpp
/// @brief Generic modal configuration; no Kokkos dependency.
/// @author EX-aero contributors
/// @details History: 2026-10 generic modal backend. Input/output: in-memory
/// configuration only. Usage: initializeModal(package, config).
#include <exaero/AttributeQuery.hpp>
#include <optional>
#include <string>
#include <vector>

namespace exaero {
/// @brief Interpretation of a table's size axis. Unknown is never guessed.
enum class ModalTableKind { Unknown, ParticleResolved, LognormalAverage };
/// @brief Select constituent addition or a supplied internally mixed provider.
enum class ModalOpticalTreatment { ComponentTables, InternalMixture };

/// @brief Optical semantics supplied alongside an existing GEOSmie curve.
/// Radius must represent dry particle radius or dry geometric mean radius,
/// respectively. Coefficients are per kg of dry material and already include
/// the source's RH response. No additional water-growth factor is applied.
struct ModalTableBinding {
  SpeciesCurveConfig curve;
  ModalTableKind kind = ModalTableKind::Unknown;
  double geometric_std_dev = 0;         ///< Required for LognormalAverage.
  std::vector<double> wavelength_nodes; ///< Physical wavelengths [m], one per
                                        ///< source spectral node. Required:
                                        ///< baked RRTMG axes contain band IDs.
  std::string data_file; ///< Optional native EX-aero MIE override file.
};

/// @brief One internally mixed material component; molecular weight in kg/mol.
struct ModalSpeciesConfig {
  std::string name;
  double density = 0;          ///< kg/m3 of dry material.
  double molecular_weight = 0; ///< kg/mol (explicitly SI).
  double hygroscopicity = 0;   ///< Kappa, dimensionless and nonnegative.
  std::optional<ModalTableBinding> optics;
};

/// @brief One fixed-width lognormal mode; radius evolves with mass and number.
struct ModalModeConfig {
  std::string name;
  double geometric_radius = 0;      ///< Reference dry geometric radius [m].
  double geometric_std_dev = 0;     ///< Fixed geometric standard deviation > 1.
  std::vector<std::string> species; ///< Permitted species; empty means all.
  ModalOpticalTreatment optical_treatment =
      ModalOpticalTreatment::ComponentTables;
};

/// @brief Destination for an intermode Brownian collision (unordered pair).
struct ModalCoagulationRoute {
  std::string first, second, destination;
};
/// @brief Growth-tail transfer; routes must form an acyclic graph.
struct ModalTransferRoute {
  std::string source, destination;
  double radius_threshold = 0; ///< Dry particle radius [m].
};
/// @brief MKPP-owned empirical sulfuric-acid nucleation, J=A*C^p.
/// C is H2SO4 molecules/m3, J is particles/m3/s; A consequently has units
/// (particles/m3/s)/(molecules/m3)^p. No chemistry is integrated by EX-aero.
struct ModalNucleationConfig {
  std::string destination;
  double prefactor = 0;
  double exponent = 1;
  double seed_radius =
      0; ///< Dry geometric radius [m], using destination width.
  std::vector<double> mass_fractions; ///< Species order; sum must equal one.
  double precursor_molecular_weight = 0.098079; ///< H2SO4 kg/mol.
};
/// @brief Complete configuration, identical for YAML and structured callers.
struct ModalConfig {
  std::vector<ModalSpeciesConfig> species;
  std::vector<ModalModeConfig> modes;
  bool coagulation = true;
  std::vector<ModalCoagulationRoute> coagulation_routes;
  std::vector<ModalTransferRoute> transfer_routes;
  std::optional<ModalNucleationConfig> nucleation;
};
} // namespace exaero
