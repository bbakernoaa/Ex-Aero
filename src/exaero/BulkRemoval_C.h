#pragma once
/// @file BulkRemoval_C.h
/// @brief C ABI for explicit bulk removal utilities; errmsg needs 256 bytes.
#include <exaero/IAerosolPackage_C.h>
enum { EXAERO_BULK_SIZE_EMPTY = -1 };
enum { EXAERO_BULK_SIZE_CONFIGURED = 0 };
enum {
  EXAERO_BULK_RADIUS_UNAVAILABLE = 0,
  EXAERO_BULK_RADIUS_CONFIGURED = 1,
  EXAERO_BULK_RADIUS_DERIVED = 2
};
#ifdef __cplusplus
extern "C" {
#endif
/// @brief Diagnose number (cell,level,species) [particles/m3] from mass.
void exaero_bulk_fixed_number(exaero_package_t pkg, int nc, int nl, int nf,
                              const double *state, double *number, char *errmsg,
                              int *errflg);
/// @brief Physical dry lower/upper/effective radius and code, both (nf,3).
void exaero_bulk_reference_radii(exaero_package_t pkg, int nf, double *radii,
                                 int *status, char *errmsg, int *errflg);
/// @brief Dry/wet diameters (cell,level,species,2) [m] and status.
void exaero_bulk_scavenging_sizes(exaero_package_t pkg, int nc, int nl, int nf,
                                  const double *rh, const double *state,
                                  double *diameters, int *status, char *errmsg,
                                  int *errflg);
/// @brief Route returned kg/m3/s with zero-based source/destination indices.
void exaero_bulk_resuspension(exaero_package_t pkg, int nc, int nl, int nf,
                              int nroute, const int *sources,
                              const int *destinations, const double *rates,
                              double *routed, char *errmsg, int *errflg);
/// @brief Rebin selected-species kg/m2 into bins and below/above tails.
void exaero_bulk_rebin_deposition(exaero_package_t pkg, int nc, int nl, int nf,
                                  int species_index, int nbin,
                                  const double *budget, const double *edges,
                                  double *bins, double *tails, char *errmsg,
                                  int *errflg);
#ifdef __cplusplus
}
#endif
