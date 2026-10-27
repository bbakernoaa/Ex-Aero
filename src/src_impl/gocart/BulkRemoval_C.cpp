/// @file BulkRemoval_C.cpp
/// @brief Exception-safe C adapters for explicit bulk removal utilities.
#include <cstring>
#include <exaero/BulkRemoval.hpp>
#include <exaero/BulkRemoval_C.h>
#include <exception>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
void require(bool ok) {
  if (!ok)
    throw std::invalid_argument(
        "FATAL ERROR: invalid bulk removal C arguments");
}
template <class F> void guarded(char *error, int *flag, F action) {
  try {
    action();
    if (error)
      error[0] = '\0';
    if (flag)
      *flag = 0;
    return;
  } catch (const std::exception &e) {
    if (error) {
      std::strncpy(error, e.what(), 255);
      error[255] = '\0';
    }
  } catch (...) {
    if (error)
      std::strcpy(error, "FATAL ERROR: bulk removal failed");
  }
  if (flag)
    *flag = 1;
}
exaero::IAerosolPackage &package(exaero_package_t handle, int nc, int nl,
                                 int nf) {
  require(handle && nc > 0 && nl > 0 && nf > 0 &&
          nc <= std::numeric_limits<int>::max() / nl / nf);
  auto &p = *static_cast<exaero::IAerosolPackage *>(handle);
  require(p.stateFields().size() == static_cast<std::size_t>(nf));
  return p;
}
} // namespace

extern "C" {
void exaero_bulk_fixed_number(exaero_package_t h, int nc, int nl, int nf,
                              const double *state, double *number, char *err,
                              int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf);
    require(state && number);
    exaero::View3D<double> out(number, nc, nl, nf);
    exaero::bulkFixedNumber(p, exaero::View3D<const double>(state, nc, nl, nf),
                            out);
  });
}
void exaero_bulk_reference_radii(exaero_package_t h, int nf, double *radii,
                                 int *status, char *err, int *flag) {
  guarded(err, flag, [&] {
    require(h && nf > 0 && radii && status &&
            nf <= std::numeric_limits<int>::max() / 3);
    auto &p = *static_cast<exaero::IAerosolPackage *>(h);
    require(p.stateFields().size() == static_cast<std::size_t>(nf));
    exaero::View2D<double> r(radii, nf, 3);
    exaero::View2D<int> st(status, nf, 3);
    exaero::bulkReferenceRadii(p, r, st);
  });
}
void exaero_bulk_scavenging_sizes(exaero_package_t h, int nc, int nl, int nf,
                                  const double *rh, const double *state,
                                  double *diameters, int *status, char *err,
                                  int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf);
    require(rh && state && diameters && status);
    exaero::View4D<double> d(diameters, nc, nl, nf, 2);
    exaero::View3D<int> st(status, nc, nl, nf);
    exaero::bulkScavengingSizes(p, exaero::View2D<const double>(rh, nc, nl),
                                exaero::View3D<const double>(state, nc, nl, nf),
                                d, st);
  });
}
void exaero_bulk_resuspension(exaero_package_t h, int nc, int nl, int nf,
                              int nroute, const int *sources,
                              const int *destinations, const double *rates,
                              double *routed, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf);
    require(nroute >= 0 && rates && routed &&
            (nroute == 0 || (sources && destinations)));
    std::vector<exaero::BulkResuspensionRoute> routes;
    routes.reserve(static_cast<std::size_t>(nroute));
    for (int i = 0; i < nroute; ++i) {
      require(sources[i] >= 0 && destinations[i] >= 0);
      routes.push_back({static_cast<std::size_t>(sources[i]),
                        static_cast<std::size_t>(destinations[i])});
    }
    exaero::View3D<double> out(routed, nc, nl, nf);
    exaero::redistributeBulkResuspension(
        p, exaero::View3D<const double>(rates, nc, nl, nf), routes, out);
  });
}
void exaero_bulk_rebin_deposition(exaero_package_t h, int nc, int nl, int nf,
                                  int species, int nbin, const double *budget,
                                  const double *edges, double *bins,
                                  double *tails, char *err, int *flag) {
  guarded(err, flag, [&] {
    auto &p = package(h, nc, nl, nf);
    require(species >= 0 && nbin > 0 &&
            nbin < std::numeric_limits<int>::max() &&
            nc <= std::numeric_limits<int>::max() / nl / nbin && budget &&
            edges && bins && tails);
    exaero::View3D<double> b(bins, nc, nl, nbin), t(tails, nc, nl, 2);
    exaero::rebinBulkDeposition(
        p, exaero::View3D<const double>(budget, nc, nl, nf),
        static_cast<std::size_t>(species),
        exaero::View1D<const double>(edges, nbin + 1), b, t);
  });
}
}
