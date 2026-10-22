/// @file modal_mam4_reference.cpp
/// @brief Native MAM4xx oracle and Kokkos lifecycle, isolated from mdspan.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: comparison driver. Host scalar kernels;
/// no file I/O or device-view dereferencing. Native arrays retain upstream ABI.
#include "modal_mam4_reference.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <mam4xx/aging.hpp>
#include <mam4xx/coagulation.hpp>
#include <mam4xx/gasaerexch.hpp>
#include <mam4xx/nucleation.hpp>
#include <numbers>
#include <stdexcept>
#ifdef EXAERO_COMPARISON_MKPP
#include "modal_mam4_comparison/host_declarations.hpp"
#endif

namespace modal_comparison {
double diffusivity(double temperature, double pressure) {
  return mam4::gasaerexch::gas_diffusivity(temperature, pressure / 101325.,
                                           98.079, 42.88);
}
std::array<double, 4> uptake(double temperature, double diffusion,
                             std::span<const double, 4> radius,
                             std::span<const double, 4> width) {
  double diameter[4], lnwidth[4];
  std::array<double, 4> result{};
  for (int m = 0; m < 4; ++m) {
    diameter[m] = 2 * radius[m];
    lnwidth[m] = std::log(width[m]);
  }
  const double speed =
      mam4::gasaerexch::mean_molecular_speed(temperature, 98.079, 0, 0);
  mam4::gasaerexch::gas_aer_uptkrates_1box1gas(
      .65, diffusion, 3 * diffusion / speed, 0, diameter, lnwidth,
      result.data());
  return result;
}
double nucleation(double acid_cm, int choice) {
  double j = 0, logj = -1e100, ntot = 0, nacid = 0, nh3 = 0, radius = 0;
  int actual = 0;
  mam4::nucleation::pbl_nuc_wang2008(acid_cm, std::numbers::pi, choice, 1,
                                     actual, j, logj, ntot, nacid, nh3, radius);
  if (actual != choice || !std::isfinite(j) || j < 0)
    throw std::runtime_error("native Wang cluster-rate evaluation failed");
  return j * 1e6;
}
void aging(double temperature, double pressure, std::span<double, 32> state,
           std::span<const double, 4> radius,
           std::span<const double, 32> condensation,
           std::span<const double, 32> coagulation_gain, double dt) {
  // Native aging's effective molecular weights are explicit, including its
  // 115 g/mol sulfate convention. Convert every species back on output.
  constexpr double mw_kg_kmol[7] = {150, 115, 150, 12, 58.44, 135, 250093};
  const double air = pressure / (1000 * mam4::Constants::r_gas * temperature);
  double mass[7][4]{}, cond[7][4]{}, coag[7][4]{}, incoming[7][1]{};
  double number[4]{}, ncond[4]{}, ncoag[4]{}, diameter[4]{};
  for (int m = 0; m < 4; ++m) {
    number[m] = state[m * 8 + 7] / air;
    diameter[m] = 2 * radius[m];
    for (int s = 0; s < 7; ++s) {
      mass[s][m] = state[m * 8 + s] / (air * mw_kg_kmol[s]);
      cond[s][m] = condensation[m * 8 + s] * dt / (air * mw_kg_kmol[s]);
    }
  }
  for (int s = 0; s < 7; ++s)
    incoming[s][0] = coagulation_gain[3 * 8 + s] / (air * mw_kg_kmol[s]);
  const auto species = mam4::default_aero_species();
  const auto host_species = [&](int s) { return species[s]; };
  mam4::aging::mam_pcarbon_aging_1subarea(host_species, 8, diameter, number,
                                          ncond, ncoag, mass, cond, coag,
                                          incoming);
  for (int m = 0; m < 4; ++m) {
    state[m * 8 + 7] = number[m] * air;
    for (int s = 0; s < 7; ++s)
      state[m * 8 + s] = mass[s][m] * air * mw_kg_kmol[s];
  }
}
#ifdef EXAERO_COMPARISON_MKPP
std::array<double, 39> generatedRhs(std::span<const double, 39> state,
                                    double temperature) {
  namespace generated = mkpp::generated::modal_mam4_comparison::detail;
  std::array<double, 39> rates{};
  generated::compute_rates_chunk_0(
      generated::StrideView<const double, 1>(state.data(), 1),
      generated::StrideView<double, 1>(rates.data(), 1), nullptr, temperature,
      0, 0);
  generated::compute_rates_chunk_1(
      generated::StrideView<const double, 1>(state.data(), 1),
      generated::StrideView<double, 1>(rates.data(), 1), nullptr, temperature,
      0, 0);
  return rates;
}
#endif
std::array<double, 32> coagulation(double temperature, double pressure,
                                   std::span<double, 32> state,
                                   std::span<const double, 4> radius,
                                   std::span<const double, 4> density,
                                   double dt) {
  const double air = pressure / (mam4::Constants::r_gas * temperature);
  double number[4]{}, mass[7][4]{}, transfer[7][1]{}, diameter[4]{};
  std::array<double, 32> gains{};
  for (int m = 0; m < 4; ++m) {
    diameter[m] = 2 * radius[m];
    number[m] = state[m * 8 + 7] / air;
    for (int s = 0; s < 7; ++s)
      mass[s][m] = state[m * 8 + s];
  }
  mam4::coagulation::mam_coag_1subarea(dt, temperature, pressure, air, diameter,
                                       density.data(), number, mass, transfer);
  for (int s = 0; s < 7; ++s) {
    gains[s] = std::max(0., mass[s][0] - state[s]);
    gains[24 + s] = transfer[s][0];
  }
  for (int m = 0; m < 4; ++m) {
    state[m * 8 + 7] = number[m] * air;
    for (int s = 0; s < 7; ++s)
      state[m * 8 + s] = mass[s][m];
  }
  return gains;
}
} // namespace modal_comparison

/// @brief Initialize Kokkos before public-package calls and release it on exit.
int main(int argc, char **argv) {
  Kokkos::ScopeGuard guard(argc, argv);
  return runComparison(argc, argv);
}
