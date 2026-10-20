#pragma once
/// @file Mam4xxConfig.hpp
/// @brief Stable modal packed-state layout; no MAM4xx or Kokkos dependency.
#include <array>
#include <string_view>

namespace exaero::mam4xx_layout {
inline constexpr int modes = 4; ///< MAM4 mode order.
inline constexpr int species =
    7; ///< Global mam4::AeroId order, not mode-local order.
inline constexpr int fields_per_mode =
    species + 1; ///< Seven masses and one number.
inline constexpr int fields_per_phase = modes * fields_per_mode;
inline constexpr int fields =
    2 * fields_per_phase; ///< Interstitial, then cloudborne.
inline constexpr std::array<std::string_view, modes> mode_names = {
    "accumulation", "aitken", "coarse", "primary_carbon"};
inline constexpr std::array<std::string_view, species> species_names = {
    "SOA", "SO4", "POM", "BC", "NaCl", "DST", "MOM"};
/// @brief Mass field [kg/m3]; aerosol is a global AeroId (0..6).
constexpr int mass(int mode, int aerosol, bool cloudborne = false) {
  return (cloudborne ? fields_per_phase : 0) + mode * fields_per_mode + aerosol;
}
/// @brief Number field [particles/m3], independently prognosed for every mode.
constexpr int number(int mode, bool cloudborne = false) {
  return mass(mode, species, cloudborne);
}
} // namespace exaero::mam4xx_layout
