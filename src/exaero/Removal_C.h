#pragma once
/// @file Removal_C.h
/// @brief Common host sinks; column-major arrays, no accepted-state mutation.
#include <exaero/IAerosolPackage_C.h>
#ifdef __cplusplus
extern "C" {
#endif
/// @brief Mass inputs/output (nc,nl,nf); number output (nc,nl,np).
/// @details Number rate (nc,nl,np) required for modal, NULL for bulk. Modal
/// number slots in mass rate must be zero. Outputs kg/m3/s and particles/m3/s.
/// Failure preserves disjoint output buffers; errmsg requires 256 bytes.
void exaero_removal_rate_tendencies(
    exaero_package_t, int nc, int nl, int nf, int np, const double *state,
    const double *mass_rate, const double *number_rate, double *mass_tendency,
    double *number_tendency, char *errmsg, int *errflg);
/// @brief Population velocity (nc,nl,np) [m/s], thickness (nc,nl) [m].
void exaero_removal_velocity_tendencies(
    exaero_package_t, int nc, int nl, int nf, int np, const double *state,
    const double *velocity, const double *thickness, double *mass_tendency,
    double *number_tendency, char *errmsg, int *errflg);
#ifdef __cplusplus
}
#endif
