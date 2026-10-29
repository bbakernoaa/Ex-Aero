#pragma once
/// @file DepositionCoupling.hpp
/// @brief Common host-return routing and deposited-budget binning; no removal
/// scheme.
#include <exaero/HostCoupling.hpp>
#include <span>
namespace exaero {
struct ResuspensionFieldRoute {
  std::size_t source_field, destination_field;
};
/// @brief Configured modal routes as mass fields; bulk has no configured
/// routes.
std::vector<ResuspensionFieldRoute>
configuredResuspensionRoutes(IAerosolPackage &);
/// @brief Route original nonnegative kg/m3/s sources once, preserving material
/// mass.
/// @details State-shaped input/output; modal number slots zero. Routes cannot
/// duplicate sources, cycle, target prohibited material or cross modal species.
/// Bulk field-to-field chemistry identity is the host's explicit mapping
/// choice. No particle number inferred; disjoint outputs publish atomically.
void redistributeResuspension(IAerosolPackage &, const View3D<const double> &,
                              std::span<const ResuspensionFieldRoute>,
                              View3D<double> &);
/// @brief Rebin deposited budgets by dry diameter edges [m], with explicit
/// tails.
/// @details Species index is package-local: one bulk bin or one modal material
/// across modes. Budgets use kg/m2 and modal particles/m2 in number fields.
/// Output bin mass (cell,level,nbin) and tails (cell,level,2), kg/m2.
void rebinDeposition(IAerosolPackage &, const View3D<const double> &,
                     std::size_t species_index,
                     const View1D<const double> &diameter_edges,
                     View3D<double> &bin_mass, View3D<double> &tail_mass);
} // namespace exaero
