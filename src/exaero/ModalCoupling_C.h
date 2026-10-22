#pragma once
/** @file ModalCoupling_C.h
 * @brief Additive generic modal coupling ABI, host column-major arrays.
 * @author EX-aero contributors
 * @details History: 2026-10. Usage: opaque modal package handle; no file I/O.
 */
#include <exaero/IAerosolPackage_C.h>
#ifdef __cplusplus
extern "C" {
#endif
/** @brief Mode count after initialization; -1 for unsupported/uninitialized. */
int exaero_get_modal_mode_count(exaero_package_t pkg);
/** @brief Query (nc,nl,nm,9) properties; optional jac is (nc,nl,nm,9*nf).
 * Properties: dry radius [m], wet radius [m], number [m^-3], mass [kg/m3],
 * dry volume [m3/m3], kappa, wet area [m2/m3], wet density [kg/m3], water
 * [kg/m3].
 * @param errmsg At least 256 bytes. @param errflg Zero on success.
 * Every environmental buffer is (nc,nl); all buffers are host memory.
 */
void exaero_modal_properties(exaero_package_t pkg, int nc, int nl, int nf,
                             int nm, const double *temp, const double *rh,
                             const double *state, double *properties,
                             double *jac, char *errmsg, int *errflg);
/** @brief Dry Fuchs-Sutugin uptake [s^-1], (nc,nl,nm,nv).
 * Gas arrays have nv elements: diffusivity [m2/s], molecular weight [kg/mol],
 * accommodation (0,1]. Optional jac: (nc,nl,nm,nv*nf).
 */
void exaero_modal_uptake(exaero_package_t pkg, int nc, int nl, int nf, int nm,
                         int nv, const double *temp, const double *rh,
                         const double *state, const double *diffusivity,
                         const double *molecular_weight,
                         const double *accommodation, double *uptake,
                         double *jac, char *errmsg, int *errflg);
/** @brief Evaluate MKPP-owned seed source terms without advancing any state.
 * acid is kg H2SO4/m3. mass_rates has ns configured species, kg/m3/s.
 * destination is zero-based; number_rate is particles/m3/s; loss is kg/m3/s.
 */
void exaero_modal_nucleation(exaero_package_t pkg, double acid, int ns,
                             int *destination, double *number_rate,
                             double *precursor_loss, double *mass_rates,
                             char *errmsg, int *errflg);
/** @brief Infer emitted number from ns species masses [kg/m3] and mode sizing.
 * The same linear conversion applies to consistent mass-rate or flux units.
 */
void exaero_modal_emitted_number(exaero_package_t pkg, int mode, int ns,
                                 const double *masses, double *number,
                                 char *errmsg, int *errflg);
#ifdef __cplusplus
}
#endif
