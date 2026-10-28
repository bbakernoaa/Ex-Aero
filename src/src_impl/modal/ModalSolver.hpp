#pragma once
/// @file ModalSolver.hpp
/// @brief Private pointer bridge; neither mdspan nor Kokkos headers cross it.
/// @author EX-aero contributors
/// @details History: 2026-10. Host buffers use (cell,level,field), layout_left.
/// Usage: package adapter. No file I/O.
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace exaero::modal {
/// @brief Device-independent configuration prepared by the validated adapter.
struct Parameters {
  std::vector<double> density, kappa, radius, sigma;
  std::vector<int> allowed, destinations, transfer_source, transfer_destination;
  std::vector<double> transfer_radius;
  bool coagulation = true;
  bool mam4xx_water = false;
  double gravity = 9.80665, water_molecular_weight = .01801528;
  double gas_constant = 8.314462618, surface_tension = .076;
  std::vector<double> crystallization_rh, deliquescence_rh;
};
/// @brief Immutable host environmental buffers, (cell,level), layout_left.
struct Environment {
  const double *temperature, *pressure, *density, *rh, *thickness;
};
/// @brief Geometry-only vapor properties in SI units.
struct Vapor {
  double diffusivity, molecular_weight, accommodation;
};
/// @brief Immutable optical table snapshot in radius-major, RH, wavelength
/// order.
struct Table {
  int kind = 0; ///< 1 particle-resolved, 2 lognormal-average.
  double sigma = 0;
  std::string source_label, source_version, citation;
  std::vector<double> radius, rh, wavelength, extinction, scattering, asymmetry,
      backscatter;
  /// Independent source-band axis and coefficients. Empty means this table
  /// contains no native band data; wavelength nodes are never substituted.
  std::vector<double> band_source_ids, band_extinction, band_scattering,
      band_asymmetry, band_backscatter;
};
/// @brief Opaque Kokkos configuration owner; release before Kokkos
/// finalization.
class Solver;
/// @brief Upload a validated configuration to the default execution space.
std::shared_ptr<Solver> create(const Parameters &, const std::vector<Table> &);
/// @brief Advance coagulation and directed growth transfer atomically.
void advance(const Solver &, std::size_t, std::size_t, const Environment &,
             double *, double, double *coagulation_gain = nullptr,
             bool do_coagulation = true, bool do_transfer = true);
/// @brief Modal properties and optional analytic derivatives on host arrays.
void properties(const Solver &, std::size_t, std::size_t, const Environment &,
                const double *, double *, double *);
/// @brief Fuchs-Sutugin dry uptake coefficients and optional derivatives.
void uptake(const Solver &, std::size_t, std::size_t, const Environment &,
            const double *, const std::vector<Vapor> &, double *, double *);
/// @brief Fill standard integral diagnostics, including level-zero columns.
void diagnostics(const Solver &, std::size_t, std::size_t, const Environment &,
                 const double *, double *);
/// @brief Lognormal kappa-Kohler CCN spectra.
void ccn(const Solver &, std::size_t, std::size_t, const Environment &,
         const double *, const double *, std::size_t, double *);
/// @brief Component-table mode integrals: (cell,level,mode,band,9),
/// layout_left.
void optics(const Solver &, std::size_t, std::size_t, const Environment &,
            const double *, const double *, std::size_t, double *);
} // namespace exaero::modal
