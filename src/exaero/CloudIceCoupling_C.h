#pragma once
/// @file CloudIceCoupling_C.h
/// @brief C/Fortran cloud and ice diagnostics; errmsg requires 256 bytes.
#include <exaero/IAerosolPackage_C.h>
enum {
  EXAERO_ICE_DUST = 0,
  EXAERO_ICE_SULFATE = 1,
  EXAERO_ICE_BLACK_CARBON = 2,
  EXAERO_ICE_STRATOSPHERIC_SULFATE = 3
};
enum {
  EXAERO_CLOUD_ICE_ABSENT = -2,
  EXAERO_CLOUD_ICE_EMPTY = -1,
  EXAERO_CLOUD_ICE_NOT_CONFIGURED = 0,
  EXAERO_CLOUD_ICE_CONFIGURED = 1
};
enum {
  EXAERO_CLOUD_ICE_MASS = 0,
  EXAERO_CLOUD_ICE_NUMBER = 1,
  EXAERO_CLOUD_ICE_DRY_NUMBER_MEDIAN_RADIUS = 2,
  EXAERO_CLOUD_ICE_DRY_SCAVENGING_RADIUS = 3,
  EXAERO_CLOUD_ICE_WET_SCAVENGING_RADIUS = 4,
  EXAERO_CLOUD_ICE_LOWER_RADIUS = 5,
  EXAERO_CLOUD_ICE_UPPER_RADIUS = 6,
  EXAERO_CLOUD_ICE_EFFECTIVE_RADIUS = 7,
  EXAERO_CLOUD_ICE_TYPE_WEIGHT = 8,
  EXAERO_CLOUD_ICE_SIZE_WEIGHT = 9,
  EXAERO_CLOUD_ICE_ELIGIBLE_NUMBER = 10,
  EXAERO_CLOUD_ICE_HETEROGENEOUS_FREEZING = 11,
  EXAERO_CLOUD_ICE_PROPERTY_COUNT = 12
};
#ifdef __cplusplus
extern "C" {
#endif
/// @brief Common population cloud/ice slabs (nc,nl,np,12); cloud_state nullable.
void exaero_cloud_ice(exaero_package_t,int nc,int nl,int nf,int np,
    const double *temperature,const double *rh,const double *state,const double *cloud_state,
    int category,int preexisting_ice,double *ambient_out,double *cloud_out,
    int *ambient_status,int *cloud_status,char *errmsg,int *errflg);
/// @brief Common fractions (nc,nl,np,2*nss), number then mass per supersaturation.
void exaero_activation_fractions(exaero_package_t,int nc,int nl,int nf,int np,int nss,
    const double *temperature,const double *rh,const double *state,const double *supersaturation,
    double *fractions,char *errmsg,int *errflg);
/// @brief Bulk outputs (nc,nl,nf,12); cloud_state may be NULL.
void exaero_bulk_cloud_ice(
    exaero_package_t pkg, int nc, int nl, int nf, const double *temperature,
    const double *rh, const double *state, const double *cloud_state,
    int category, int preexisting_ice, double *ambient_out, double *cloud_out,
    int *ambient_status, int *cloud_status, char *errmsg, int *errflg);
/// @brief Modal outputs (nc,nl,nm,12); nf is packed mass+number fields.
void exaero_modal_cloud_ice(
    exaero_package_t pkg, int nc, int nl, int nf, const double *temperature,
    const double *rh, const double *state, const double *cloud_state,
    int category, int preexisting_ice, double *ambient_out, double *cloud_out,
    int *ambient_status, int *cloud_status, char *errmsg, int *errflg);
/// @brief Interstitial fractions (nc,nl,nf,2*nss), number then mass.
void exaero_bulk_activation_fractions(
    exaero_package_t pkg, int nc, int nl, int nf, int nss,
    const double *temperature, const double *rh, const double *state,
    const double *supersaturation, double *fractions, char *errmsg,
    int *errflg);
/// @brief Interstitial fractions (nc,nl,nm,2*nss), number then mass.
void exaero_modal_activation_fractions(
    exaero_package_t pkg, int nc, int nl, int nf, int nss,
    const double *temperature, const double *rh, const double *state,
    const double *supersaturation, double *fractions, char *errmsg,
    int *errflg);
#ifdef __cplusplus
}
#endif
