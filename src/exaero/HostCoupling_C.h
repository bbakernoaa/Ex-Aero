#pragma once
/// @file HostCoupling_C.h
/// @brief Zero-based host metadata snapshots; caller buffers, no state updates.
#include <exaero/IAerosolPackage_C.h>
enum { EXAERO_NUMBER_DIAGNOSED=0, EXAERO_NUMBER_PROGNOSTIC=1 };
enum { EXAERO_HOST_POPULATION_PROPERTIES=0, EXAERO_HOST_CLOUD_ICE=1,
  EXAERO_HOST_ACTIVATION=2, EXAERO_HOST_REMOVAL_RATES=3,
  EXAERO_HOST_REMOVAL_VELOCITY=4, EXAERO_HOST_RESUSPENSION=5,
  EXAERO_HOST_DEPOSITION_REBIN=6, EXAERO_HOST_INSPECT_NUMBER_BOUNDS=7,
  EXAERO_HOST_ENFORCE_NUMBER_BOUNDS=8 };
#ifdef __cplusplus
extern "C" {
#endif
/// @brief Errors preserve descriptor outputs; errmsg requires 256 bytes.
void exaero_host_coupling_counts(exaero_package_t,int *representation,int *nf,int *np,int *nconstituent,char *errmsg,int *errflg);
void exaero_host_population_info(exaero_package_t,int population,int *index,int *number_kind,int *number_field,char *name,int name_capacity,char *errmsg,int *errflg);
void exaero_host_constituent_info(exaero_package_t,int constituent,int *population,int *species,int *mass_field,int *permitted,char *name,int name_capacity,char *errmsg,int *errflg);
void exaero_supports_host_operation(exaero_package_t,int operation,int *supported,char *errmsg,int *errflg);
#ifdef __cplusplus
}
#endif
