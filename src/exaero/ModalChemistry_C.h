#pragma once
/** @file ModalChemistry_C.h
 * @brief Additive host column-major phase-exchange/aging coupling ABI.
 * @author EX-aero contributors
 * @details History: 2026-10. No integration or file I/O; error buffers >=256
 * bytes.
 */
#include <exaero/IAerosolPackage_C.h>
#ifdef __cplusplus
extern "C" {
#endif
/** @brief Balanced phase exchange. Species indices are zero-based, length nv.
 * gas(nc,nl,nv) and eq(nc,nl,nm,nv): kg target-material-equivalent/m3.
 * ar(nc,nl,nf), gr(nc,nl,nv): tendencies per second.
 * Optional jac(nc,nl,nf+nv,ni), eq_jac(nc,nl,nm*nv,ni), ni>=nf+nv.
 * Inputs/outputs cannot overlap. See ModalChemistry.hpp for active branches.
 */
void exaero_modal_phase_exchange(
    exaero_package_t pkg, int nc, int nl, int nf, int nm, int nv, int ni,
    const int *species, const double *diff, const double *mw,
    const double *alpha, const double *temp, const double *rh,
    const double *state, const double *gas, const double *eq, double *ar,
    double *gr, double *jac, const double *eq_jac, char *errmsg, int *errflg);
/** @brief Aging RHS from signed condensation and gross coagulation mass rates.
 * rates(nc,nl,nf); optional jac(nc,nl,nf,ni), ni>=nf, aerosol columns first.
 * Optional budget derivatives have cond_rows/coag_rows rows (nf..ni), ni
 * columns. Number entries of both budgets must be zero. No state update is
 * performed.
 */
void exaero_modal_aging_terms(exaero_package_t pkg, int nc, int nl, int nf,
                              int ni, const double *temp, const double *rh,
                              const double *state, const double *cond,
                              const double *coag, double *rates, double *jac,
                              const double *cond_jac, int cond_rows,
                              const double *coag_jac, int coag_rows,
                              char *errmsg, int *errflg);
/** @brief Advance microphysics and overwrite gross incoming coagulation mass
 * gain(nc,nl,nf), kg/m3 over dt; its number entries are zero. Both outputs
 * atomic.
 */
void exaero_modal_microphysics_with_coating(
    exaero_package_t pkg, int nc, int nl, int nf, const double *temp,
    const double *press, const double *rho, const double *rh, const double *dz,
    double *state, double dt, double *gain, char *errmsg, int *errflg);
/** @brief Host-owned growth-tail renaming, once after an accepted chemistry
 * interval. */
void exaero_modal_rename_modes(exaero_package_t pkg, int nc, int nl, int nf,
                               const double *temp, const double *press,
                               const double *rho, const double *rh,
                               const double *dz, double *state, char *errmsg,
                               int *errflg);
#ifdef __cplusplus
}
#endif
