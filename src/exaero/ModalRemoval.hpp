#pragma once
/// @file ModalRemoval.hpp
/// @brief Explicit host utilities for generic modal removal and number bounds.
/// @details All arrays are column-major; operations publish outputs atomically.
#include <cstddef>
#include <exaero/IAerosolPackage.hpp>

namespace exaero {
enum class ModalBoundStatus : int {
  Disabled = 0,
  Empty = 1,
  Within = 2,
  BelowRadius = 3,
  AboveRadius = 4
};

/// @brief Report dry radius and allowed number interval without changing state.
/// @param[in] package Initialized generic modal package.
/// @param[in] state Packed mass and number, (cell,level,field).
/// @param[out] report Radius, Nmin, Nmax, (cell,level,mode,3).
/// @param[out] status ModalBoundStatus, (cell,level,mode).
void inspectModalBounds(IAerosolPackage &package,
                        const View3D<const double> &state,
                        View4D<double> &report, View3D<int> &status);

/// @brief Clamp number to configured dry radius bounds, preserving mass.
/// @param[in] package Initialized generic modal package.
/// @param[in,out] state Packed aerosol state; only number fields may change.
/// @param[out] delta_number Signed number corrections, (cell,level,mode).
/// @param[out] status Status of original state, (cell,level,mode).
void enforceModalBounds(IAerosolPackage &package, View3D<double> &state,
                        View3D<double> &delta_number, View3D<int> &status);

/// @brief Return dry/wet number-median diameters for host removal models.
/// @param[in] package Initialized generic modal package.
/// @param[in] env Temperature and RH supplied by host.
/// @param[in] state Packed aerosol state.
/// @param[out] diameters (cell,level,mode,2), metres.
/// @param[out] populated (cell,level,mode), 0 or 1.
void modalScavengingSizes(IAerosolPackage &package,
                          const EnvironmentalStateView &env,
                          const View3D<const double> &state,
                          View4D<double> &diameters, View3D<int> &populated);

/// @brief Map host first-order mode/moment coefficients to removal tendencies.
/// @details Synchronous host column-major views. State is packed mass/number
/// (cell,level,field). Both coefficient inputs are (cell,level,mode), in [1/s].
/// A mode's mass coefficient removes all its mixed constituents proportionally;
/// its independent number coefficient may represent a different distribution
/// moment. The host derives physically consistent coefficients from current
/// composition, geometry, activation and precipitation; this computes no rates.
/// Mass output uses state layout [kg/m3/s], with number slots zero. Number
/// output is (cell,level,mode) [particles/m3/s]. Inputs remain unchanged.
/// Outputs must be disjoint from inputs and each other; failure is atomic.
void modalRemovalRateTendencies(
    IAerosolPackage &package, const View3D<const double> &state,
    const View3D<const double> &mass_scavenging_rate,
    const View3D<const double> &number_scavenging_rate,
    View3D<double> &mass_tendency, View3D<double> &number_tendency);

/// @brief Return mode-uniform tendencies -velocity/layer_thickness times state.
/// @details Velocity [m/s] is (cell,level,mode); thickness [m] is (cell,level),
/// finite and positive. The same velocity applies to mass and number. Layout,
/// output units and atomicity match modalRemovalRateTendencies. Hosts needing
/// size-selective moment losses supply separate coefficients to that API.
/// No bounds correction, timestep integration, or state update is performed.
void modalRemovalVelocityTendencies(IAerosolPackage &package,
                                    const View2D<const double> &layer_thickness,
                                    const View3D<const double> &state,
                                    const View3D<const double> &velocity,
                                    View3D<double> &mass_tendency,
                                    View3D<double> &number_tendency);

/// @brief Route nonnegative returned material rates once from a snapshot.
/// @param[in] package Initialized generic modal package.
/// @param[in] mass_rates (cell,level,field), kg/m3/s, zero number slots.
/// @param[out] routed_rates Same shape and units; no particle number inferred.
void redistributeModalResuspension(IAerosolPackage &package,
                                   const View3D<const double> &mass_rates,
                                   View3D<double> &routed_rates);

/// @brief Rebin deposited material budgets by dry geometric diameter.
/// @param[in] package Initialized generic modal package.
/// @param[in] deposited_budget (cell,level,field), kg/m2 or particles/m2.
/// @param[in] species_index Zero-based configured dry material index.
/// @param[in] diameter_edges nbin+1 increasing edges in metres.
/// @param[out] bin_mass (cell,level,nbin), kg/m2.
/// @param[out] tail_mass (cell,level,2), below/above bin range, kg/m2.
void rebinModalDeposition(IAerosolPackage &package,
                          const View3D<const double> &deposited_budget,
                          std::size_t species_index,
                          const View1D<const double> &diameter_edges,
                          View3D<double> &bin_mass, View3D<double> &tail_mass);
} // namespace exaero
