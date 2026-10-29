/// @file ModalRemoval_C.cpp
/// @brief Checked C adapters for generic modal removal operations.
#include <cstring>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemoval.hpp>
#include <exaero/ModalRemoval_C.h>
#include <exception>
#include <limits>
#include <stdexcept>

namespace {
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
      std::strcpy(error, "FATAL ERROR: modal removal failed");
  }
  if (status)
    *status = 1;
}
void require(bool ok) {
  if (!ok)
    throw std::invalid_argument(
        "FATAL ERROR: invalid modal removal C arguments");
}
exaero::IAerosolPackage &package(exaero_package_t handle, int nc, int nl,
                                 int nf, int *nm = nullptr) {
  require(handle && nc > 0 && nl > 0 && nf > 0 &&
          static_cast<std::size_t>(nc) <=
              static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl /
                  nf);
  auto &p = *static_cast<exaero::IAerosolPackage *>(handle);
  const auto &cfg = exaero::modalCoupling(p).modalConfig();
  require(cfg.modes.size() * (cfg.species.size() + 1) ==
          static_cast<std::size_t>(nf));
  if (nm)
    require(*nm > 0 && static_cast<std::size_t>(*nm) == cfg.modes.size());
  return p;
}
exaero::EnvironmentalStateView environment(int nc, int nl, const double *t,
                                           const double *rh) {
  require(t && rh);
  return {exaero::View2D<const double>(t, nc, nl),
          exaero::View2D<const double>(nullptr, nc, nl),
          exaero::View2D<const double>(nullptr, nc, nl),
          exaero::View2D<const double>(rh, nc, nl),
          exaero::View2D<const double>(nullptr, nc, nl)};
}
} // namespace
extern "C" {
void exaero_modal_removal_rate_tendencies(exaero_package_t h, int nc, int nl,
                                          int nf, int nm, const double *state,
                                          const double *mass_rate,
                                          const double *number_rate,
                                          double *mass, double *number,
                                          char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf, &nm);
    require(state && mass_rate && number_rate && mass && number);
    exaero::View3D<double> mv(mass, nc, nl, nf), nv(number, nc, nl, nm);
    exaero::modalRemovalRateTendencies(
        p, exaero::View3D<const double>(state, nc, nl, nf),
        exaero::View3D<const double>(mass_rate, nc, nl, nm),
        exaero::View3D<const double>(number_rate, nc, nl, nm), mv, nv);
  });
}
void exaero_modal_removal_velocity_tendencies(
    exaero_package_t h, int nc, int nl, int nf, int nm, const double *state,
    const double *velocity, const double *thickness, double *mass,
    double *number, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf, &nm);
    require(state && velocity && thickness && mass && number);
    exaero::View3D<double> mv(mass, nc, nl, nf), nv(number, nc, nl, nm);
    exaero::modalRemovalVelocityTendencies(
        p, exaero::View2D<const double>(thickness, nc, nl),
        exaero::View3D<const double>(state, nc, nl, nf),
        exaero::View3D<const double>(velocity, nc, nl, nm), mv, nv);
  });
}
void exaero_modal_inspect_bounds(exaero_package_t h, int nc, int nl, int nf,
                                 int nm, const double *state, double *report,
                                 int *status, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf, &nm);
    require(state && report && status);
    exaero::View4D<double> r(report, nc, nl, nm, 3);
    exaero::View3D<int> st(status, nc, nl, nm);
    exaero::inspectModalBounds(
        p, exaero::View3D<const double>(state, nc, nl, nf), r, st);
  });
}
void exaero_modal_enforce_bounds(exaero_package_t h, int nc, int nl, int nf,
                                 int nm, double *state, double *delta,
                                 int *status, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf, &nm);
    require(state && delta && status);
    exaero::View3D<double> a(state, nc, nl, nf), d(delta, nc, nl, nm);
    exaero::View3D<int> st(status, nc, nl, nm);
    exaero::enforceModalBounds(p, a, d, st);
  });
}
void exaero_modal_scavenging_sizes(exaero_package_t h, int nc, int nl, int nf,
                                   int nm, const double *t, const double *rh,
                                   const double *state, double *diameters,
                                   int *pop, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf, &nm);
    require(state && diameters && pop);
    exaero::View4D<double> d(diameters, nc, nl, nm, 2);
    exaero::View3D<int> q(pop, nc, nl, nm);
    exaero::modalScavengingSizes(
        p, environment(nc, nl, t, rh),
        exaero::View3D<const double>(state, nc, nl, nf), d, q);
  });
}
void exaero_modal_resuspension(exaero_package_t h, int nc, int nl, int nf,
                               const double *rates, double *routed, char *err,
                               int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf);
    require(rates && routed);
    exaero::View3D<double> out(routed, nc, nl, nf);
    exaero::redistributeModalResuspension(
        p, exaero::View3D<const double>(rates, nc, nl, nf), out);
  });
}
void exaero_modal_rebin_deposition(exaero_package_t h, int nc, int nl, int nf,
                                   int species, int nbin, const double *budget,
                                   const double *edges, double *bins,
                                   double *tails, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf);
    require(species >= 0 && nbin > 0 &&
            nbin < std::numeric_limits<int>::max() && budget && edges && bins &&
            tails &&
            static_cast<std::size_t>(nc) <=
                static_cast<std::size_t>(std::numeric_limits<int>::max()) / nl /
                    nbin);
    exaero::View3D<double> b(bins, nc, nl, nbin), t(tails, nc, nl, 2);
    exaero::rebinModalDeposition(
        p, exaero::View3D<const double>(budget, nc, nl, nf),
        static_cast<std::size_t>(species),
        exaero::View1D<const double>(edges, nbin + 1), b, t);
  });
}
}
