#pragma once
/// @file GocartSpeciesParams.hpp
/// @brief POD species parameters + curve-pool offsets for the GOCART solver.
///
/// This is the transfer record between configuration (host) and the
/// Kokkos kernels (device). It is a plain-old-data type carrying only
/// scalars and *offsets/extents* into the package's single flat device
/// curve pool — no fixed-size arrays, so every table axis length is data
/// supplied by config or the pinned GEOSmie tables, never code.
#include <exaero/AerosolIndices.hpp>

namespace exaero {

/// @brief Raw-emissions -> package-species mapping parameters.
///
/// One entry per configured CECE raw species: how much of the raw mass
/// is routed to this package species, and (optionally) the modal
/// distribution used to convert emitted mass into particle number.
struct EmissionsMappingParams {
  bool is_active = false;  ///< False => this species receives no emissions.
  int raw_cece_index = -1; ///< Index of the raw species in the input flux
                           ///< axis (cell, level, raw_species).
  double mass_split_fraction = 0.0; ///< Fraction [0,1] of raw mass routed
                                    ///< to this species.

  // Modal Parameters for Mass-to-Number conversion
  bool is_modal_mode = false; ///< True => number emission derived from a
                              ///< lognormal mode rather than a fixed ratio.
  double emitted_particle_diameter = 1.0; ///< Mode mean diameter [m].
  double lognormal_sigma = 1.0; ///< Geometric standard deviation (GSD).
};

/// @brief Per-species physical parameters and curve-pool bindings.
///
/// The eight analytical scalars drive the Mode A (ADT-style) fallback
/// path; the offset/extents pairs bind the species to table-backed
/// curves (RH optics lookup, microphysical growth, spectral MIE, and
/// polarized moments) inside the device curve pool. When an offset is
/// -1 the corresponding table is absent and the solver falls back to
/// the analytical path — it never silently clamps to an unrelated curve.
struct GocartSpeciesParams {
  double dry_density;           ///< [kg/m³] particle material density (rho_0).
  double molecular_weight;      ///< [g/mol] for Köhler solute terms.
  double dry_particle_diameter; ///< [m] dry mode diameter (D_0).
  double hygroscopicity;        ///< kappa (k) hygroscopicity parameter, [-].
  double lognormal_sigma;       ///< Geometric standard deviation (GSD), [-].
  double lognormal_dg;          ///< Geometric mean diameter (GMD) [m].
  double refractive_index_real; ///< Real part n of the complex RI, [-].
  double refractive_index_imag; ///< Imaginary part k of the complex RI, [-].

  // RH optical lookup (Mode B). The tabulated curves live in ONE flat device
  // double pool ; this POD carries only the offset + extents,
  // never a fixed array. Block layout at pool[curve_offset...]: four
  // contiguous n_rh-long slots in the documented order {rh_axis, ext, ssa,
  // asm}. Lengths are data, not capacity.
  bool has_optics_lookup; ///< True => Mode B RH lookup; false => Mode A
                          ///< analytical ADT path.
  int curve_offset =
      -1;           ///< Index into the device curve pool; -1 => no curve bound.
  int n_radius = 0; ///< Radius/bin extent (1 for band-integrated lookup).
  int n_rh = 0;     ///< RH-axis length (arbitrary, from config or table).
  int n_lambda =
      0; ///< Wavelength/band extent (0 for band-integrated legacy lookup).

  // Table-backed microphysical curves for the device hot path (
  //). When micro_offset >= 0 the diagnostics kernel reads hygroscopic
  // growth and wet particle density from the GEOSmie curve instead of the
  // kappa-Kohler approximation. Block layout at pool[micro_offset...]
  // (documented order, extents are data):
  // [0, nH) rh axis
  // [nH, nH + nR*nH) growth_factor (radius-major: r*nH + h)
  // [nH + nR*nH, nH + 2*nR*nH) wet_particle_density (same layout)
  int micro_offset =
      -1; ///< Index into the device curve pool; -1 => no micro curve.
  int n_micro_radius = 0; ///< Radius-bin extent of the micro block.
  int n_micro_rh = 0;     ///< RH-axis extent of the micro block.
  int solver_radius_node =
      0; ///< Size-bin the hot path queries (clamped to n_micro_radius-1).

  // Table-backed spectral curves for the optics hot path.
  // When spec_offset >= 0 the optics kernels read mass extinction,
  // single-scattering albedo and asymmetry from the GEOSmie curve at the solver
  // radius node, interpolating linearly in RH and linear-in-log over the band
  // coordinate (the RRTMG lambda axis), replacing the ADT analytical solver.
  // Block layout at pool[spec_offset...]; the extents are data (never
  // literals):
  // [0, nH) rh axis
  // [nH, nH + nL) lambda (band coordinate) axis
  // [nH + nL,... + nR*nH*nL) bext (mass extinction) (r-major:
  // (r*nH+h)*nL+l)
  // [...,... + nR*nH*nL) ssa (single-scattering albedo)
  // [...,... + nR*nH*nL) g (asymmetry factor)
  int spec_offset =
      -1; ///< Index into the device curve pool; -1 => no spectral curve.
  int n_spec_radius = 0; ///< Radius-bin extent of the spectral block.
  int n_spec_rh = 0;     ///< RH-axis extent of the spectral block.
  int n_spec_lambda = 0; ///< Band-axis extent of the spectral block.

  // Table-backed polarized phase-function moments for the device hot path.
  // When pmom_offset >= 0 the rank-5 moment
  // array is device-resident in the same single-upload pool; a polarized
  // radiative-transfer consumer reads it with zero H2D in the timestep loop.
  // Block layout at pool[pmom_offset...]; extents are data (never literals):
  // [0, nH) rh axis
  // [nH, nH + nL) lambda (band coordinate) axis
  // [nH + nL,...) pmom, C-order over
  // (radius,rh,lambda,pol,moment)
  // i.e. (((r*nH+h)*nL+l)*nP+e)*nM+m
  int pmom_offset = -1; ///< Index into the device curve pool; -1 => no moments.
  int n_pmom_radius = 0; ///< Radius-bin extent of the moment block.
  int n_pmom_rh = 0;     ///< RH-axis extent of the moment block.
  int n_pmom_lambda = 0; ///< Band-axis extent of the moment block.
  int n_pmom_pol = 0;    ///< Polarized element count (P11..P44 ordering).
  int n_pmom_moment = 0; ///< Moment count (M), per-species data.

  EmissionsMappingParams emissions_mapping; ///< Emissions routing for this
                                            ///< species (inactive by default).
};

/// @brief Named slots within a species' Mode-B curve-pool block
/// (documented order; four contiguous n_rh-long segments).
namespace curve_slot {
constexpr int RH_AXIS = 0; ///< RH axis [fraction 0-1], strictly ordered.
constexpr int EXT = 1;     ///< Extinction lookup, parallel to rh axis.
constexpr int SSA = 2;     ///< Single-scattering albedo lookup in [0, 1].
constexpr int ASM = 3;     ///< Asymmetry-factor lookup in [-1, 1].
constexpr int NUM = 4;     ///< Slot count (block stride = NUM * n_rh).
} // namespace curve_slot

} // namespace exaero
