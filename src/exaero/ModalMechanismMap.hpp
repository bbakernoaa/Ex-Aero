#pragma once
/// @file ModalMechanismMap.hpp
/// @brief Name-based, explicit EX-aero/MKPP unit and index mapping.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: build once per compiled mechanism; no I/O.
#include <algorithm>
#include <cmath>
#include <exaero/ModalCoupling.hpp>
#include <set>

namespace exaero {
/// @brief One requested aerosol-field to mechanism-variable association.
struct ModalMechanismBinding {
  std::string aerosol_field, mechanism_variable;
};
/// @brief Resolved mapping for mechanisms using molecules/cm3 and
/// particles/cm3. Multiply SI state/tendencies by to_mechanism on entry; divide
/// on return. Number remains an independent particle-number variable, never a
/// molecule.
struct ModalMechanismSlot {
  int aerosol_index = -1, mechanism_index = -1;
  StateQuantity quantity = StateQuantity::Mass;
  double to_mechanism = 0;
};
/// @brief Resolve every configured aerosol field exactly once; reject missing,
/// duplicate, or unknown names. Chemistry mechanisms may contain other species.
inline std::vector<ModalMechanismSlot>
modalMechanismMap(IAerosolPackage &package,
                  std::span<const std::string> mechanism_names,
                  std::span<const ModalMechanismBinding> bindings) {
  const auto &config = modalCoupling(package).modalConfig();
  const auto fields = package.stateFields();
  if (fields.empty() || bindings.size() != fields.size())
    throw std::invalid_argument(
        "FATAL ERROR: bind every modal state field exactly once");
  std::set<std::string> names(mechanism_names.begin(), mechanism_names.end());
  if (names.size() != mechanism_names.size() || names.count(""))
    throw std::invalid_argument(
        "FATAL ERROR: duplicate or empty mechanism name");
  std::set<int> seen_aerosol, seen_mechanism;
  std::vector<ModalMechanismSlot> result;
  for (const auto &binding : bindings) {
    const int f = package.getSpeciesIndex(binding.aerosol_field);
    const auto it = std::find(mechanism_names.begin(), mechanism_names.end(),
                              binding.mechanism_variable);
    const int g = static_cast<int>(it - mechanism_names.begin());
    if (f < 0 || it == mechanism_names.end() ||
        !seen_aerosol.insert(f).second || !seen_mechanism.insert(g).second)
      throw std::invalid_argument(
          "FATAL ERROR: invalid or duplicate modal mechanism binding");
    const bool number = fields[f].quantity == StateQuantity::Number;
    const auto s = static_cast<std::size_t>(f) % (config.species.size() + 1);
    const double scale =
        number ? 1e-6 : 6.02214076e17 / config.species.at(s).molecular_weight;
    if (!std::isfinite(scale) || scale <= 0)
      throw std::invalid_argument(
          "FATAL ERROR: mechanism unit conversion overflows");
    result.push_back({f, g, fields[f].quantity, scale});
  }
  return result;
}
/// @brief Vapor binding in the equivalent-material basis used by phase
/// exchange.
struct ModalVaporMechanismBinding {
  std::string aerosol_species, mechanism_variable;
  double vapor_molecules_per_material_unit =
      1; ///< Stoichiometric molecule/formula-unit ratio.
};
/// @brief One combined local aerosol-then-vapor to mechanism unit/index
/// mapping.
struct ModalChemicalMechanismSlot {
  int local_index = -1, mechanism_index = -1;
  StateQuantity quantity = StateQuantity::Mass;
  double to_mechanism = 0;
};
/// @brief Extend the existing aerosol mapping with explicitly ordered vapors.
/// Concentrations use molecules/cm3 (number fields particles/cm3). Gas
/// conversion uses target-material molecular weight and explicit stoichiometry,
/// not the physical vapor molecular weight used by the transport kernel. Other
/// reaction partners/coproducts belong to the host mechanism.
inline std::vector<ModalChemicalMechanismSlot> modalChemicalMechanismMap(
    IAerosolPackage &package, std::span<const std::string> mechanism_names,
    std::span<const ModalMechanismBinding> aerosol_bindings,
    std::span<const ModalVaporMechanismBinding> vapor_bindings) {
  const auto aerosol =
      modalMechanismMap(package, mechanism_names, aerosol_bindings);
  const auto &config = modalCoupling(package).modalConfig();
  std::vector<ModalChemicalMechanismSlot> result;
  std::set<int> seen_variables, seen_species;
  for (const auto &slot : aerosol) {
    result.push_back({slot.aerosol_index, slot.mechanism_index, slot.quantity,
                      slot.to_mechanism});
    seen_variables.insert(slot.mechanism_index);
  }
  for (std::size_t v = 0; v < vapor_bindings.size(); ++v) {
    const auto &binding = vapor_bindings[v];
    const auto material = std::find_if(
        config.species.begin(), config.species.end(),
        [&](const auto &s) { return s.name == binding.aerosol_species; });
    const auto variable =
        std::find(mechanism_names.begin(), mechanism_names.end(),
                  binding.mechanism_variable);
    if (material == config.species.end() || variable == mechanism_names.end() ||
        !seen_species.insert(material - config.species.begin()).second ||
        !seen_variables.insert(variable - mechanism_names.begin()).second ||
        !std::isfinite(binding.vapor_molecules_per_material_unit) ||
        binding.vapor_molecules_per_material_unit <= 0)
      throw std::invalid_argument(
          "FATAL ERROR: invalid or duplicate vapor mechanism binding");
    const double scale = 6.02214076e17 / material->molecular_weight *
                         binding.vapor_molecules_per_material_unit;
    if (!std::isfinite(scale) || scale <= 0)
      throw std::invalid_argument(
          "FATAL ERROR: vapor mechanism scale overflows");
    result.push_back({static_cast<int>(aerosol.size() + v),
                      static_cast<int>(variable - mechanism_names.begin()),
                      StateQuantity::Mass, scale});
  }
  return result;
}
/// @brief Chain-rule multiplier for J_mechanism(i,j)=scale_i/scale_j*J_SI(i,j).
inline double
modalChemicalJacobianScale(const ModalChemicalMechanismSlot &row,
                           const ModalChemicalMechanismSlot &column) {
  if (!std::isfinite(row.to_mechanism) || row.to_mechanism <= 0 ||
      !std::isfinite(column.to_mechanism) || column.to_mechanism <= 0)
    throw std::invalid_argument(
        "FATAL ERROR: invalid chemical Jacobian unit scale");
  const double ratio = row.to_mechanism / column.to_mechanism;
  if (!std::isfinite(ratio) || ratio <= 0)
    throw std::invalid_argument(
        "FATAL ERROR: invalid chemical Jacobian unit scale");
  return ratio;
}
} // namespace exaero
