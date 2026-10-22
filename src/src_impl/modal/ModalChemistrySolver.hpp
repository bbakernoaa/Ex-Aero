#pragma once
/// @file ModalChemistrySolver.hpp
/// @brief Kokkos-isolated chemistry pointer bridge; no public mdspan headers.
/// @author EX-aero contributors
/// @details History: 2026-10. Host column-major arrays, synchronous returns; no
/// I/O.
#include <cstddef>
#include <vector>
namespace exaero::modal {
/// @brief Resolved, validated aging/coating configuration.
struct AgingParameters {
  std::vector<double> density, sigma, thickness;
  std::vector<int> source, destination, coating; // route-major species mask.
};
/// @brief Balanced exchange RHS/Jacobian with precomputed dry uptake geometry.
void exchange(std::size_t nc, std::size_t nl, int ns, int nm,
              const std::vector<int> &species, const std::vector<int> &allowed,
              int ni, const double *state, const double *gas,
              const double *equilibrium, const double *uptake,
              const double *uptake_jac, const double *equilibrium_jac,
              double *aerosol_rates, double *gas_rates, double *jac);
/// @brief Shell-coating aging rate/Jacobian; all outputs are signed.
void aging(std::size_t nc, std::size_t nl, const AgingParameters &, int ni,
           const double *state, const double *condensation,
           const double *coagulation, const double *condensation_jac,
           int condensation_rows, const double *coagulation_jac,
           int coagulation_rows, double *rates, double *jac);
} // namespace exaero::modal
