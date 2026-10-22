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
} // namespace exaero
