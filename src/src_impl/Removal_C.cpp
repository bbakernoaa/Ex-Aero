/// @file Removal_C.cpp
/// @brief Shape-checked C/Fortran forwarding for atomic host removal adapters.
#include "HostCouplingInternal.hpp"
#include <exaero/Removal.hpp>
#include <exaero/Removal_C.h>
#include <optional>
using namespace exaero;
using namespace exaero::host_detail;
extern "C" {
void exaero_removal_rate_tendencies(exaero_package_t h, int nc, int nl, int nf,
                                    int np, const double *state,
                                    const double *mass_rate,
                                    const double *number_rate, double *mass,
                                    double *number, char *error, int *flag) {
  guarded(error, flag, [&] {
    auto &p = package(h);
    dimensions(p, nc, nl, nf, np);
    require(state && mass_rate && mass && number, "null removal array");
    std::optional<View3D<const double>> nr;
    if (number_rate)
      nr.emplace(number_rate, nc, nl, np);
    View3D<double> mv(mass, nc, nl, nf), nv(number, nc, nl, np);
    removalRateTendencies(p, View3D<const double>(state, nc, nl, nf),
                          View3D<const double>(mass_rate, nc, nl, nf),
                          nr ? &*nr : nullptr, mv, nv);
  });
}
void exaero_removal_velocity_tendencies(exaero_package_t h, int nc, int nl,
                                        int nf, int np, const double *state,
                                        const double *velocity,
                                        const double *dz, double *mass,
                                        double *number, char *error,
                                        int *flag) {
  guarded(error, flag, [&] {
    auto &p = package(h);
    dimensions(p, nc, nl, nf, np);
    require(state && velocity && dz && mass && number,
            "null removal velocity array");
    View3D<double> mv(mass, nc, nl, nf), nv(number, nc, nl, np);
    removalVelocityTendencies(p, View2D<const double>(dz, nc, nl),
                              View3D<const double>(state, nc, nl, nf),
                              View3D<const double>(velocity, nc, nl, np), mv,
                              nv);
  });
}
}
