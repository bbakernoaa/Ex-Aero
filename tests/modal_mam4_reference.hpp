#pragma once
/// @file modal_mam4_reference.hpp
/// @brief Host-only scalar boundary isolating Kokkos from public mdspan
/// headers.
/// @author EX-aero contributors
/// @details History: 2026-10. No I/O; spans are contiguous host memory.
#include <array>
#include <span>
namespace modal_comparison {
/// @brief Native dry uptake per particle, m3/s, using native molecular speed.
std::array<double, 4> uptake(double temperature, double diffusivity,
                             std::span<const double, 4> radius,
                             std::span<const double, 4> width);
/// @brief Native vapor diffusivity in m2/s for H2SO4.
double diffusivity(double temperature, double pressure);
/// @brief Native Wang cluster formation, particles/m3/s.
double nucleation(double acid_cm, int choice);
/// @brief Direct native coagulation oracle; dry RH=0, input geometric radii
/// and mixture densities. Returns accepted gross coating gains in kg/m3.
std::array<double, 32> coagulation(double temperature, double pressure,
                                   std::span<double, 32> state,
                                   std::span<const double, 4> radius,
                                   std::span<const double, 4> density,
                                   double dt);
/// @brief Native full finite-step aging; state is mode-major, seven masses and
/// number per mode. Coating gains are integrated kg/m3; cond is kg/m3/s.
void aging(double temperature, double pressure, std::span<double, 32> state,
           std::span<const double, 4> radius,
           std::span<const double, 32> condensation,
           std::span<const double, 32> coagulation_gain, double dt);
#ifdef EXAERO_COMPARISON_MKPP
/// @brief Execute the compiled comparison YAML in molecular/cm3 coordinates.
std::array<double, 39> generatedRhs(std::span<const double, 39> state,
                                    double temperature);
#endif
} // namespace modal_comparison
/// @brief Public mdspan driver called within a separately owned Kokkos scope.
int runComparison(int argc, char **argv);
