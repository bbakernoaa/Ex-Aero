#pragma once
/// @file StateField.hpp
/// @brief State-axis metadata shared by bulk and modal packages.
#include <string>

namespace exaero {
/// @brief Aerosol population representation; numeric values also form the C
/// ABI.
enum class AerosolRepresentation { Bulk = 0, Modal = 1 };
/// @brief Quantity in a packed state slot (concentrations, not mixing ratios).
enum class StateQuantity { Mass = 0, Number = 1 };
/// @brief Metadata for one (cell, level, field) column-major state slot.
struct StateField {
  std::string name;                             ///< Unique lookup key.
  StateQuantity quantity = StateQuantity::Mass; ///< kg/m3 or particles/m3.
  int mode = -1; ///< -1 for bulk; MAM4 order: accumulation, Aitken, coarse,
                 ///< primary carbon.
  bool cloudborne = false; ///< Cloudborne rather than interstitial population.
};
/// @brief Supported operations, queryable before calling an optional surface.
enum PackageCapability : unsigned {
  Microphysics = 1,
  Diagnostics = 2,
  Optics = 4,
  CCN = 8,
  Emissions = 16
};
} // namespace exaero
