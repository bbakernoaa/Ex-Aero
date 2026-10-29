/// @file CommonCloudIce_C.cpp
/// @brief Validated common cloud/ice and activation C bindings.
#include "HostCouplingInternal.hpp"
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceCoupling_C.h>
#include <optional>
using namespace exaero;
using namespace exaero::host_detail;
extern "C" {
void exaero_cloud_ice(exaero_package_t h, int nc, int nl, int nf, int np,
                      const double *t, const double *rh, const double *state,
                      const double *cloud, int category, int ice, double *a,
                      double *b, int *as, int *bs, char *error, int *flag) {
  guarded(error, flag, [&] {
    auto &p = package(h);
    dimensions(p, nc, nl, nf, np, cloud_ice_properties::COUNT);
    require(state && a && b && as && bs &&
                validIceCategory(static_cast<IceCategory>(category)) &&
                (ice == 0 || ice == 1),
            "invalid common cloud/ice arguments");
    std::optional<View3D<const double>> cv;
    if (cloud)
      cv.emplace(cloud, nc, nl, nf);
    View4D<double> av(a, nc, nl, np, cloud_ice_properties::COUNT),
        bv(b, nc, nl, np, cloud_ice_properties::COUNT);
    View3D<int> asi(as, nc, nl, np), bsi(bs, nc, nl, np);
    cloudIce(p, environment(nc, nl, t, rh),
             View3D<const double>(state, nc, nl, nf), cv ? &*cv : nullptr,
             static_cast<IceCategory>(category), ice == 1, av, bv, asi, bsi);
  });
}
void exaero_activation_fractions(exaero_package_t h, int nc, int nl, int nf,
                                 int np, int nss, const double *t,
                                 const double *rh, const double *state,
                                 const double *ss, double *out, char *error,
                                 int *flag) {
  guarded(error, flag, [&] {
    auto &p = package(h);
    require(nss > 0 && nss <= std::numeric_limits<int>::max() / 2,
            "invalid supersaturation count");
    dimensions(p, nc, nl, nf, np, 2 * nss);
    require(state && ss && out, "null activation array");
    View4D<double> result(out, nc, nl, np, 2 * nss);
    activationFractions(p, environment(nc, nl, t, rh),
                        View3D<const double>(state, nc, nl, nf),
                        View1D<const double>(ss, nss), result);
  });
}
}
