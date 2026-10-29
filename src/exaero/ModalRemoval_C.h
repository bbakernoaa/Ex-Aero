#pragma once
/// @file ModalRemoval_C.h
/// @brief C ABI for explicit generic modal removal utilities.
#include <exaero/IAerosolPackage_C.h>
/// @brief Integer codes returned in modal bound status buffers.
enum {
  EXAERO_MODAL_BOUNDS_DISABLED = 0,
  EXAERO_MODAL_BOUNDS_EMPTY = 1,
  EXAERO_MODAL_BOUNDS_WITHIN = 2,
  EXAERO_MODAL_BOUNDS_BELOW_RADIUS = 3,
  EXAERO_MODAL_BOUNDS_ABOVE_RADIUS = 4
};
#ifdef __cplusplus
extern "C" {
#endif
/// @brief All error buffers below require at least 256 characters.
void exaero_modal_inspect_bounds(exaero_package_t pkg, int nc, int nl, int nf,
                                 int nm, const double *state, double *report,
                                 int *status, char *errmsg, int *errflg);
void exaero_modal_enforce_bounds(exaero_package_t pkg, int nc, int nl, int nf,
                                 int nm, double *state, double *delta,
                                 int *status, char *errmsg, int *errflg);
void exaero_modal_scavenging_sizes(exaero_package_t pkg, int nc, int nl, int nf,
                                   int nm, const double *temp, const double *rh,
                                   const double *state, double *diameters,
                                   int *populated, char *errmsg, int *errflg);
/// @brief Host first-order mass/number coefficients (nc,nl,nm) [1/s] to sinks.
/// @details Host column-major state/mass output are (nc,nl,nf); number output
/// is (nc,nl,nm). Mass coefficient is common to a mode's constituents; number
/// is independent. Output units are kg/m3/s and particles/m3/s; packed number
/// slots in mass output are zero. Outputs cannot overlap inputs or each other.
/// No state update or bounds correction; errors preserve both output buffers.
void exaero_modal_removal_rate_tendencies(
    exaero_package_t pkg, int nc, int nl, int nf, int nm, const double *state,
    const double *mass_scavenging_rate, const double *number_scavenging_rate,
    double *mass_tendency, double *number_tendency, char *errmsg, int *errflg);
/// @brief Mode velocities (nc,nl,nm) [m/s] and positive thickness (nc,nl) [m].
/// @details Return -velocity/thickness times each mode's mass and number;
/// output shape, units, non-overlap and atomicity match the rate adapter.
void exaero_modal_removal_velocity_tendencies(
    exaero_package_t pkg, int nc, int nl, int nf, int nm, const double *state,
    const double *velocity, const double *layer_thickness,
    double *mass_tendency, double *number_tendency, char *errmsg, int *errflg);
void exaero_modal_resuspension(exaero_package_t pkg, int nc, int nl, int nf,
                               const double *rates, double *routed,
                               char *errmsg, int *errflg);
void exaero_modal_rebin_deposition(exaero_package_t pkg, int nc, int nl, int nf,
                                   int species_index, int nbin,
                                   const double *budget, const double *edges,
                                   double *bin_mass, double *tail_mass,
                                   char *errmsg, int *errflg);
#ifdef __cplusplus
}
#endif
