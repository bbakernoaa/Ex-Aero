#pragma once
#include <exaero/AerosolIndices.hpp>

namespace exaero {

struct EmissionsMappingParams {
  bool is_active = false;
  int raw_cece_index = -1;
  double mass_split_fraction = 0.0;

  // Modal Parameters for Mass-to-Number conversion
  bool is_modal_mode = false;
  double emitted_particle_diameter = 1.0;
  double lognormal_sigma = 1.0;
};

struct GocartSpeciesParams {
  double dry_density;           // [kg/m³]
  double molecular_weight;      // [g/mol]
  double dry_particle_diameter; // [m]
  double hygroscopicity;        // kappa value
  double lognormal_sigma;       // GSD
  double lognormal_dg;          // GMD [m]
  double refractive_index_real; // n
  double refractive_index_imag; // k

  // RH optical lookup (Mode B). The tabulated curves live in ONE flat device
  // double pool (ADR-003 R10); this POD carries only the offset + extents,
  // never a fixed array. Block layout at pool[curve_offset ...]: four
  // contiguous n_rh-long slots in the documented order {rh_axis, ext, ssa,
  // asm}. Lengths are data, not capacity.
  bool has_optics_lookup; // Flag to enable RH lookup tables (Mode A ADT when
                          // false)
  int curve_offset =
      -1;           // index into the device curve pool; -1 => no curve bound
  int n_radius = 0; // radius/bin extent (1 for band-integrated legacy lookup)
  int n_rh = 0;     // RH-axis length (arbitrary, from config or table)
  int n_lambda =
      0; // wavelength/band extent (0 for band-integrated legacy lookup)

  // Table-backed microphysical curves for the device hot path (T021, ADR-003
  // R10). When micro_offset >= 0 the diagnostics kernel reads hygroscopic
  // growth and wet particle density from the GEOSmie curve instead of the
  // kappa-Kohler approximation. Block layout at pool[micro_offset ...]
  // (documented order, extents are data):
  //   [0, nH)                     rh axis
  //   [nH, nH + nR*nH)            growth_factor   (radius-major: r*nH + h)
  //   [nH + nR*nH, nH + 2*nR*nH)  wet_particle_density (same layout)
  int micro_offset =
      -1; // index into the device curve pool; -1 => no micro curve
  int n_micro_radius = 0; // radius-bin extent of the micro block
  int n_micro_rh = 0;     // RH-axis extent of the micro block
  int solver_radius_node =
      0; // size-bin the hot path queries (clamped to n_micro_radius-1)

  // Table-backed spectral curves for the optics hot path (T028, ADR-003 R10).
  // When spec_offset >= 0 the optics kernels read mass extinction,
  // single-scattering albedo and asymmetry from the GEOSmie curve at the solver
  // radius node, interpolating linearly in RH and linear-in-log over the band
  // coordinate (the RRTMG lambda axis), replacing the ADT analytical solver.
  // Block layout at pool[spec_offset ...]; the extents are data (never
  // literals):
  //   [0, nH)                          rh axis
  //   [nH, nH + nL)                    lambda (band coordinate) axis
  //   [nH + nL, ... + nR*nH*nL)        bext (mass extinction)  (r-major:
  //   (r*nH+h)*nL+l)
  //   [... , ... + nR*nH*nL)           ssa  (single-scattering albedo)
  //   [... , ... + nR*nH*nL)           g    (asymmetry factor)
  int spec_offset =
      -1; // index into the device curve pool; -1 => no spectral curve
  int n_spec_radius = 0; // radius-bin extent of the spectral block
  int n_spec_rh = 0;     // RH-axis extent of the spectral block
  int n_spec_lambda = 0; // band-axis extent of the spectral block

  // Table-backed polarized phase-function moments for the device hot path
  // (T033, FR-003/FR-015, ADR-003 R10). When pmom_offset >= 0 the rank-5 moment
  // array is device-resident in the same single-upload pool; a polarized
  // radiative-transfer consumer reads it with zero H2D in the timestep loop.
  // Block layout at pool[pmom_offset ...]; extents are data (never literals):
  //   [0, nH)                              rh axis
  //   [nH, nH + nL)                        lambda (band coordinate) axis
  //   [nH + nL, ... )                      pmom, C-order over
  //   (radius,rh,lambda,pol,moment)
  //                                        i.e. (((r*nH+h)*nL+l)*nP+e)*nM+m
  int pmom_offset = -1;  // index into the device curve pool; -1 => no moments
  int n_pmom_radius = 0; // radius-bin extent of the moment block
  int n_pmom_rh = 0;     // RH-axis extent of the moment block
  int n_pmom_lambda = 0; // band-axis extent of the moment block
  int n_pmom_pol = 0;    // polarized element count (P11..P44 ordering)
  int n_pmom_moment = 0; // moment count (M)

  EmissionsMappingParams emissions_mapping;
};

// Named slots within a species' curve-pool block (documented order; R10).
namespace curve_slot {
constexpr int RH_AXIS = 0;
constexpr int EXT = 1;
constexpr int SSA = 2;
constexpr int ASM = 3;
constexpr int NUM = 4;
} // namespace curve_slot

} // namespace exaero
