#pragma once
/// @file DepositionCoupling_C.h
/// @brief Common host budget utilities; zero-based fields/species, no state
/// update.
#include <exaero/IAerosolPackage_C.h>
#ifdef __cplusplus
extern "C" {
#endif
/// @brief capacity=0 queries count with NULL arrays; insufficient capacity is
/// atomic.
void exaero_resuspension_routes(exaero_package_t, int capacity, int *nroute,
                                int *sources, int *destinations, char *errmsg,
                                int *errflg);
/// @brief Route state-shaped material kg/m3/s from original source fields once.
void exaero_redistribute_resuspension(exaero_package_t, int nc, int nl, int nf,
                                      int nroute, const int *sources,
                                      const int *destinations,
                                      const double *rates, double *routed,
                                      char *errmsg, int *errflg);
/// @brief Species/bin bulk selection or modal material across modes; budget
/// kg/m2.
/// @details Modal number fields use particles/m2. Outputs bins (nc,nl,nbin),
/// tails (nc,nl,2); dry diameter edges (nbin+1) in m. errmsg requires 256
/// bytes.
void exaero_rebin_deposition(exaero_package_t, int nc, int nl, int nf,
                             int species_index, int nbin, const double *budget,
                             const double *edges, double *bin_mass,
                             double *tail_mass, char *errmsg, int *errflg);
#ifdef __cplusplus
}
#endif
