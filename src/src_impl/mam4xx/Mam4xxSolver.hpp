#pragma once
/// @file Mam4xxSolver.hpp
/// @brief Private host-pointer bridge, isolated from both mdspan and MAM4xx
/// headers.
#include <cstddef>

namespace exaero {
/// @brief Read-only contiguous host buffers, each (cell, level), layout_left.
struct Mam4xxEnvironment {
  const double *temperature; ///< Air temperature [K].
  const double *pressure;    ///< Air pressure [Pa]; unused for CCN.
  const double
      *density;     ///< Air density [kg/m3]; validated on the public boundary.
  const double *rh; ///< Relative humidity [fraction].
  const double *thickness; ///< Layer depth [m]; unused for CCN.
};
/// @brief Run upstream clear-air coagulation on a host (cell, level, field)
/// state. Stages to the default device, fences, and copies back before
/// returning.
void mam4xxMicrophysics(std::size_t cells, std::size_t levels,
                        const Mam4xxEnvironment &env, double *state, double dt);
/// @brief Compute modal integral diagnostics into a host layout_left slab.
void mam4xxDiagnostics(std::size_t cells, std::size_t levels,
                       const Mam4xxEnvironment &env, const double *state,
                       double *output);
/// @brief Modal lognormal kappa-Kohler CCN spectrum [particles/m3].
void mam4xxCCN(std::size_t cells, std::size_t levels,
               const Mam4xxEnvironment &env, const double *state,
               const double *supersaturations, std::size_t count,
               double *output);
} // namespace exaero
