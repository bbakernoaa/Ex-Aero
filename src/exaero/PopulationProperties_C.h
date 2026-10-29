#pragma once
/// @file PopulationProperties_C.h
/// @brief Physical property slabs (cell,level,population,property), SI units.
#include <exaero/IAerosolPackage_C.h>
enum {
  EXAERO_PROPERTY_UNAVAILABLE = 0,
  EXAERO_PROPERTY_CONFIGURED = 1,
  EXAERO_PROPERTY_DERIVED = 2
};
enum { EXAERO_POPULATION_EMPTY = 0, EXAERO_POPULATION_POPULATED = 1 };
enum {
  EXAERO_POP_MASS = 0,
  EXAERO_POP_NUMBER = 1,
  EXAERO_POP_DRY_NUMBER_MEDIAN_RADIUS = 2,
  EXAERO_POP_WET_NUMBER_MEDIAN_RADIUS = 3,
  EXAERO_POP_DRY_SCAVENGING_DIAMETER = 4,
  EXAERO_POP_WET_SCAVENGING_DIAMETER = 5,
  EXAERO_POP_DRY_VOLUME = 6,
  EXAERO_POP_DRY_SURFACE_AREA = 7,
  EXAERO_POP_KAPPA = 8,
  EXAERO_POP_DRY_DENSITY = 9,
  EXAERO_POP_WET_DENSITY = 10,
  EXAERO_POP_GEOMETRIC_STD_DEV = 11,
  EXAERO_POP_LOWER_DRY_RADIUS = 12,
  EXAERO_POP_UPPER_DRY_RADIUS = 13,
  EXAERO_POP_EFFECTIVE_DRY_RADIUS = 14,
  EXAERO_POP_PROPERTY_COUNT = 15
};
#ifdef __cplusplus
extern "C" {
#endif
/// @brief Values and per-value availability (nc,nl,np,15); status (nc,nl,np).
/// @details All outputs disjoint; atomic failure. errmsg requires 256 bytes.
void exaero_population_properties(exaero_package_t, int nc, int nl, int nf,
                                  int np, const double *temperature,
                                  const double *rh, const double *state,
                                  double *values, int *availability,
                                  int *population_status, char *errmsg,
                                  int *errflg);
#ifdef __cplusplus
}
#endif
