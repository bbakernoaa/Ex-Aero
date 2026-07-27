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

  bool has_optics_lookup; // Flag to enable RH lookup tables
  double rh_bins[8];      // RH values
  double ext_lookup[8];   // Pre-tabulated mass extinction [m²/g]
  double ssa_lookup[8];   // Pre-tabulated SSA [fraction]
  double asm_lookup[8];   // Pre-tabulated asymmetry [g]

  EmissionsMappingParams emissions_mapping;
};

} // namespace exaero
