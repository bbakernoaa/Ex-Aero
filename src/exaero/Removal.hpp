#pragma once
/// @file Removal.hpp
/// @brief Common synchronous removal adapters; the host owns rates and updates.
#include <exaero/IAerosolPackage.hpp>
namespace exaero {
/// @brief Apply host nonnegative /s coefficients; no integration or state
/// update.
/// @details State, mass rates and mass output use (cell,level,field). Modal
/// mass-rate number slots must be zero. Number rates are required for modal,
/// absent for bulk; number output uses (cell,level,population). Bulk number is
/// diagnostic, modal number independent. Outputs kg/m3/s and particles/m3/s.
/// Outputs must be disjoint from all inputs and each other; failure is atomic.
void removalRateTendencies(IAerosolPackage &, const View3D<const double> &state,
                           const View3D<const double> &mass_rate,
                           const View3D<const double> *number_rate,
                           View3D<double> &mass_tendency,
                           View3D<double> &number_tendency);
/// @brief Map population velocities [m/s] and positive layer thickness [m].
/// @details Velocity (cell,level,population); output layouts as above. The same
/// fractional rate v/dz applies to a population's mass constituents and number.
void removalVelocityTendencies(IAerosolPackage &,
                               const View2D<const double> &thickness,
                               const View3D<const double> &state,
                               const View3D<const double> &velocity,
                               View3D<double> &mass_tendency,
                               View3D<double> &number_tendency);
} // namespace exaero
