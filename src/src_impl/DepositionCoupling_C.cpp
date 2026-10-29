/// @file DepositionCoupling_C.cpp
/// @brief Checked budget routing and rebinning for C/Fortran hosts.
#include "HostCouplingInternal.hpp"
#include <exaero/DepositionCoupling.hpp>
#include <exaero/DepositionCoupling_C.h>
using namespace exaero;
using namespace exaero::host_detail;
extern "C" {
void exaero_resuspension_routes(exaero_package_t h, int capacity, int *nroute,
                                int *sources, int *destinations, char *error,
                                int *flag) {
  guarded(error, flag, [&] {
    require(nroute && capacity >= 0, "invalid route-query buffers");
    const auto routes = configuredResuspensionRoutes(package(h));
    if (capacity == 0) {
      *nroute = int(routes.size());
      return;
    }
    require(std::size_t(capacity) >= routes.size(),
            "route-query capacity too small");
    require(routes.empty() || (sources && destinations),
            "null route-query buffer");
    if (!routes.empty()) {
      separate(sources, routes.size() * sizeof(int), destinations,
               routes.size() * sizeof(int));
      separate(sources, routes.size() * sizeof(int), nroute, sizeof(int));
      separate(destinations, routes.size() * sizeof(int), nroute, sizeof(int));
    }
    for (std::size_t i = 0; i < routes.size(); ++i) {
      sources[i] = int(routes[i].source_field);
      destinations[i] = int(routes[i].destination_field);
    }
    *nroute = int(routes.size());
  });
}
void exaero_redistribute_resuspension(exaero_package_t h, int nc, int nl,
                                      int nf, int nroute, const int *sources,
                                      const int *destinations,
                                      const double *rates, double *out,
                                      char *error, int *flag) {
  guarded(error, flag, [&] {
    auto &p = package(h);
    const auto info = hostCouplingInfo(p);
    dimensions(p, nc, nl, nf, int(info.populations.size()));
    require(rates && out && nroute >= 0 && nroute <= nf &&
                (nroute == 0 || (sources && destinations)),
            "invalid resuspension C arrays");
    if (nroute > 0) {
      separate(sources, std::size_t(nroute) * sizeof(int), out,
               std::size_t(nc) * nl * nf * sizeof(double));
      separate(destinations, std::size_t(nroute) * sizeof(int), out,
               std::size_t(nc) * nl * nf * sizeof(double));
    }
    std::vector<ResuspensionFieldRoute> routes;
    for (int i = 0; i < nroute; ++i) {
      require(sources[i] >= 0 && destinations[i] >= 0, "negative route field");
      routes.push_back({std::size_t(sources[i]), std::size_t(destinations[i])});
    }
    View3D<double> result(out, nc, nl, nf);
    redistributeResuspension(p, View3D<const double>(rates, nc, nl, nf), routes,
                             result);
  });
}
void exaero_rebin_deposition(exaero_package_t h, int nc, int nl, int nf,
                             int species, int nbin, const double *budget,
                             const double *edges, double *bins, double *tails,
                             char *error, int *flag) {
  guarded(error, flag, [&] {
    auto &p = package(h);
    const auto info = hostCouplingInfo(p);
    require(species >= 0 && nbin > 0 &&
                nbin < std::numeric_limits<int>::max() && budget && edges &&
                bins && tails,
            "invalid rebin C arguments");
    dimensions(p, nc, nl, nf, int(info.populations.size()), std::max(nbin, 2));
    View3D<double> b(bins, nc, nl, nbin), t(tails, nc, nl, 2);
    rebinDeposition(p, View3D<const double>(budget, nc, nl, nf),
                    std::size_t(species), View1D<const double>(edges, nbin + 1),
                    b, t);
  });
}
}
