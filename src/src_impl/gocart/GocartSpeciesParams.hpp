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

  // RH optical lookup (Mode B). The tabulated curves live in ONE flat device double
  // pool (ADR-003 R10); this POD carries only the offset + extents, never a fixed array.
  // Block layout at pool[curve_offset ...]: four contiguous n_rh-long slots in the
  // documented order {rh_axis, ext, ssa, asm}. Lengths are data, not capacity.
  bool has_optics_lookup; // Flag to enable RH lookup tables (Mode A ADT when false)
  int  curve_offset = -1; // index into the device curve pool; -1 => no curve bound
  int  n_radius     = 0;  // radius/bin extent (1 for band-integrated legacy lookup)
  int  n_rh         = 0;  // RH-axis length (arbitrary, from config or table)
  int  n_lambda     = 0;  // wavelength/band extent (0 for band-integrated legacy lookup)

  EmissionsMappingParams emissions_mapping;
};

// Named slots within a species' curve-pool block (documented order; R10).
namespace curve_slot {
constexpr int RH_AXIS = 0;
constexpr int EXT     = 1;
constexpr int SSA     = 2;
constexpr int ASM     = 3;
constexpr int NUM     = 4;
} // namespace curve_slot

} // namespace exaero
