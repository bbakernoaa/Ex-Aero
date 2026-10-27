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
