/// @file ModalChemistry_C.cpp
/// @brief Exception-safe additive C/Fortran chemical coupling adapters.
/// @author EX-aero contributors
/// @details History: 2026-10. Host layout_left views; no file I/O/integration.
#include <cstring>
#include <exaero/ModalChemistry.hpp>
#include <exaero/ModalChemistry_C.h>
#include <limits>
#include <vector>

namespace {
/// @brief Translate exceptions without allowing a C ABI escape.
template <class F> void guarded(char *error, int *status, F action) {
  try {
    action();
    if (error)
      error[0] = '\0';
    if (status)
      *status = 0;
    return;
  } catch (const std::exception &e) {
    if (error) {
      std::strncpy(error, e.what(), 255);
      error[255] = '\0';
    }
  } catch (...) {
    if (error)
      std::strcpy(error, "FATAL ERROR: modal chemistry failed");
  }
  if (status)
    *status = 1;
}
/// @brief Host-boundary assertion before constructing nonowning views.
void require(bool value) {
  if (!value)
    throw std::invalid_argument(
        "FATAL ERROR: invalid modal chemistry C buffers/dimensions");
}
/// @brief Validate the opaque handle and generic field layout.
exaero::IAerosolPackage &package(exaero_package_t handle, int nc, int nl,
                                 int nf) {
  require(handle && nc > 0 && nl > 0 && nf > 0);
  auto &p = *static_cast<exaero::IAerosolPackage *>(handle);
  const auto &c = exaero::modalCoupling(p).modalConfig();
  require(!c.modes.empty() && static_cast<std::size_t>(nf) ==
                                  c.modes.size() * (c.species.size() + 1));
  require(static_cast<std::size_t>(nc) <=
          static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl / nf);
  return p;
}
/// @brief Bound dense Jacobian products before constructing mdspan extents.
void matrix(int nc, int nl, int rows, int columns) {
  require(rows > 0 && columns > 0 &&
          static_cast<std::size_t>(nc) <=
              static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl /
                  rows / columns);
}
/// @brief Environment view; pressure/density/thickness optional for RHS
/// queries.
exaero::EnvironmentalStateView environment(int nc, int nl, const double *t,
                                           const double *rh,
                                           const double *press = nullptr,
                                           const double *rho = nullptr,
                                           const double *dz = nullptr) {
  return {exaero::View2D<const double>(t, nc, nl),
          exaero::View2D<const double>(press, nc, nl),
          exaero::View2D<const double>(rho, nc, nl),
          exaero::View2D<const double>(rh, nc, nl),
          exaero::View2D<const double>(dz, nc, nl)};
}
} // namespace
extern "C" {
void exaero_modal_phase_exchange(
    exaero_package_t handle, int nc, int nl, int nf, int nm, int nv, int ni,
    const int *species, const double *diff, const double *mw,
    const double *alpha, const double *temp, const double *rh,
    const double *state, const double *gas, const double *eq, double *ar,
    double *gr, double *jac, const double *eq_jac, char *error, int *status) {
  guarded(error, status, [&] {
    auto &p = package(handle, nc, nl, nf);
    const auto &c = exaero::modalCoupling(p).modalConfig();
    require(nm > 0 && static_cast<std::size_t>(nm) == c.modes.size() &&
            nv > 0 && static_cast<std::size_t>(nv) <= c.species.size());
    require(species && diff && mw && alpha);
    const int nq = nf + nv;
    require(ni >= nq);
    if (jac)
      matrix(nc, nl, nq, ni);
    if (eq_jac)
      matrix(nc, nl, nm * nv, ni);
    std::vector<exaero::ModalCondensable> vapors;
    for (int v = 0; v < nv; ++v) {
      require(species[v] >= 0 &&
              static_cast<std::size_t>(species[v]) < c.species.size());
      vapors.push_back(
          {c.species[species[v]].name, {diff[v], mw[v], alpha[v]}});
    }
    exaero::View3D<double> aout(ar, nc, nl, nf), gout(gr, nc, nl, nv);
    exaero::View4D<double> jout(jac, nc, nl, nq, ni);
    const exaero::View4D<const double> ej(eq_jac, nc, nl, nm * nv, ni);
    exaero::modalPhaseExchange(
        p, environment(nc, nl, temp, rh),
        exaero::View3D<const double>(state, nc, nl, nf), vapors,
        exaero::View3D<const double>(gas, nc, nl, nv),
        exaero::View4D<const double>(eq, nc, nl, nm, nv), aout, gout,
        jac ? &jout : nullptr, eq_jac ? &ej : nullptr);
  });
}
void exaero_modal_aging_terms(exaero_package_t handle, int nc, int nl, int nf,
                              int ni, const double *temp, const double *rh,
                              const double *state, const double *cond,
                              const double *coag, double *rates, double *jac,
                              const double *cond_jac, int cond_rows,
                              const double *coag_jac, int coag_rows,
                              char *error, int *status) {
  guarded(error, status, [&] {
    auto &p = package(handle, nc, nl, nf);
    require(ni >= nf);
    if (jac)
      matrix(nc, nl, nf, ni);
    if (cond_jac) {
      require(cond_rows >= nf && cond_rows <= ni);
      matrix(nc, nl, cond_rows, ni);
    } else
      cond_rows = 0;
    if (coag_jac) {
      require(coag_rows >= nf && coag_rows <= ni);
      matrix(nc, nl, coag_rows, ni);
    } else
      coag_rows = 0;
    exaero::View3D<double> out(rates, nc, nl, nf);
    exaero::View4D<double> jout(jac, nc, nl, nf, ni);
    const exaero::View4D<const double> dc(cond_jac, nc, nl, cond_rows, ni),
        dg(coag_jac, nc, nl, coag_rows, ni);
    exaero::modalAgingTerms(p, environment(nc, nl, temp, rh),
                            exaero::View3D<const double>(state, nc, nl, nf),
                            exaero::View3D<const double>(cond, nc, nl, nf),
                            exaero::View3D<const double>(coag, nc, nl, nf), out,
                            jac ? &jout : nullptr, cond_jac ? &dc : nullptr,
                            coag_jac ? &dg : nullptr);
  });
}
void exaero_modal_microphysics_with_coating(
    exaero_package_t handle, int nc, int nl, int nf, const double *temp,
    const double *press, const double *rho, const double *rh, const double *dz,
    double *state, double dt, double *gain, char *error, int *status) {
  guarded(error, status, [&] {
    auto &p = package(handle, nc, nl, nf);
    exaero::View3D<double> a(state, nc, nl, nf), g(gain, nc, nl, nf);
    exaero::advanceModalWithCoating(
        p, environment(nc, nl, temp, rh, press, rho, dz), a, dt, g);
  });
}
void exaero_modal_rename_modes(exaero_package_t handle, int nc, int nl, int nf,
                               const double *temp, const double *press,
                               const double *rho, const double *rh,
                               const double *dz, double *state, char *error,
                               int *status) {
  guarded(error, status, [&] {
    auto &p = package(handle, nc, nl, nf);
    exaero::View3D<double> a(state, nc, nl, nf);
    exaero::renameModalModes(p, environment(nc, nl, temp, rh, press, rho, dz),
                             a);
  });
}
} // extern C
