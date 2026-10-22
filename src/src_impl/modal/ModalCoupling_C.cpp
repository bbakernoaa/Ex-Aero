/// @file ModalCoupling_C.cpp
/// @brief Exception-safe C/Fortran modal coupling adapters.
/// @author EX-aero contributors
/// @details History: 2026-10. Host layout_left inputs/outputs, no file I/O.
#include <algorithm>
#include <cstring>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalCoupling_C.h>
#include <limits>

namespace {
/// @brief Translate exceptions to the existing CCPP error protocol.
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
      std::strcpy(error, "FATAL ERROR: modal coupling failed");
  }
  if (status)
    *status = 1;
}
/// @brief Validate the opaque handle without changing the base interface ABI.
exaero::IModalCoupling &coupling(exaero_package_t pkg) {
  if (!pkg)
    throw std::invalid_argument("FATAL ERROR: null modal handle");
  return exaero::modalCoupling(*static_cast<exaero::IAerosolPackage *>(pkg));
}
/// @brief Validate dimensions and products before constructing mdspan extents.
void dimensions(exaero_package_t pkg, int nc, int nl, int nf, int nm, int nq) {
  auto &p = coupling(pkg);
  if (nc <= 0 || nl <= 0 || nf <= 0 || nm <= 0 || nq <= 0 ||
      static_cast<std::size_t>(nm) != p.modalConfig().modes.size() ||
      static_cast<std::size_t>(nf) !=
          nm * (p.modalConfig().species.size() + 1) ||
      static_cast<std::size_t>(nc) >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl / nm /
              nq / nf)
    throw std::invalid_argument(
        "FATAL ERROR: invalid modal coupling dimensions");
}
/// @brief Properties and uptake need only temperature and RH.
exaero::EnvironmentalStateView environment(int nc, int nl, const double *t,
                                           const double *rh) {
  return {exaero::View2D<const double>(t, nc, nl),
          {},
          {},
          exaero::View2D<const double>(rh, nc, nl),
          {}};
}
} // namespace
extern "C" {
int exaero_get_modal_mode_count(exaero_package_t pkg) {
  try {
    const auto &c = coupling(pkg).modalConfig();
    return c.modes.empty() ? -1 : static_cast<int>(c.modes.size());
  } catch (...) {
    return -1;
  }
}
void exaero_modal_properties(exaero_package_t pkg, int nc, int nl, int nf,
                             int nm, const double *temp, const double *rh,
                             const double *state, double *out, double *jac,
                             char *error, int *status) {
  guarded(error, status, [&] {
    dimensions(pkg, nc, nl, nf, nm, exaero::modal_properties::COUNT);
    exaero::View4D<double> result(out, nc, nl, nm,
                                  exaero::modal_properties::COUNT);
    exaero::View4D<double> derivatives(jac, nc, nl, nm,
                                       exaero::modal_properties::COUNT * nf);
    coupling(pkg).modalProperties(
        environment(nc, nl, temp, rh),
        exaero::View3D<const double>(state, nc, nl, nf), result,
        jac ? &derivatives : nullptr);
  });
}
void exaero_modal_uptake(exaero_package_t pkg, int nc, int nl, int nf, int nm,
                         int nv, const double *temp, const double *rh,
                         const double *state, const double *diff,
                         const double *mw, const double *alpha, double *out,
                         double *jac, char *error, int *status) {
  guarded(error, status, [&] {
    dimensions(pkg, nc, nl, nf, nm, nv);
    if (!diff || !mw || !alpha)
      throw std::invalid_argument("FATAL ERROR: null vapor properties");
    std::vector<exaero::ModalVapor> vapors;
    for (int v = 0; v < nv; ++v)
      vapors.push_back({diff[v], mw[v], alpha[v]});
    exaero::View4D<double> result(out, nc, nl, nm, nv),
        derivatives(jac, nc, nl, nm, nv * nf);
    coupling(pkg).modalUptake(environment(nc, nl, temp, rh),
                              exaero::View3D<const double>(state, nc, nl, nf),
                              vapors, result, jac ? &derivatives : nullptr);
  });
}
void exaero_modal_nucleation(exaero_package_t pkg, double acid, int ns,
                             int *dest, double *number, double *loss,
                             double *mass, char *error, int *status) {
  guarded(error, status, [&] {
    auto &p = coupling(pkg);
    if (ns <= 0 ||
        static_cast<std::size_t>(ns) != p.modalConfig().species.size() ||
        !dest || !number || !loss || !mass)
      throw std::invalid_argument(
          "FATAL ERROR: invalid nucleation output buffers");
    const auto terms = p.nucleationTerms(acid);
    std::copy(terms.mass_rates.begin(), terms.mass_rates.end(), mass);
    *dest = terms.destination;
    *number = terms.number_rate;
    *loss = terms.precursor_loss;
  });
}
void exaero_modal_emitted_number(exaero_package_t pkg, int mode, int ns,
                                 const double *mass, double *number,
                                 char *error, int *status) {
  guarded(error, status, [&] {
    if (mode < 0 || ns <= 0 || !mass || !number)
      throw std::invalid_argument(
          "FATAL ERROR: invalid emitted sizing buffers");
    *number =
        coupling(pkg).emittedNumber(mode, std::span<const double>(mass, ns));
  });
}
}
