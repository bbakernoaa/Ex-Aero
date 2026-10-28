/// @file CloudIceCoupling_C.cpp
/// @brief Exception-safe zero-copy adapters for cloud/ice diagnostics.
#include <cstring>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceCoupling_C.h>
#include <exaero/ModalCoupling.hpp>
#include <limits>
#include <optional>
#include <stdexcept>

namespace {
void require(bool okay) {
  if (!okay)
    throw std::invalid_argument("FATAL ERROR: invalid cloud/ice C arguments");
}
template <class F> void guarded(char *error, int *flag, F action) {
  try {
    action();
    if (error)
      error[0] = '\0';
    if (flag)
      *flag = 0;
  } catch (const std::exception &e) {
    if (error) {
      std::strncpy(error, e.what(), 255);
      error[255] = '\0';
    }
    if (flag)
      *flag = 1;
  } catch (...) {
    if (error)
      std::strcpy(error, "FATAL ERROR: cloud/ice coupling failed");
    if (flag)
      *flag = 1;
  }
}
exaero::IAerosolPackage &package(exaero_package_t handle, int nc, int nl,
                                 int nf) {
  require(handle && nc > 0 && nl > 0 && nf > 0 &&
          nc <= std::numeric_limits<int>::max() / nl / nf);
  auto &p = *static_cast<exaero::IAerosolPackage *>(handle);
  require(p.stateFields().size() == static_cast<std::size_t>(nf));
  return p;
}
std::size_t populations(exaero::IAerosolPackage &p, int nf, bool modal) {
  if (!modal)
    return static_cast<std::size_t>(nf);
  auto *coupling = dynamic_cast<exaero::IModalCoupling *>(&p);
  require(coupling != nullptr);
  const auto &config = coupling->modalConfig();
  require(!config.modes.empty() && !config.species.empty() &&
          static_cast<std::size_t>(nf) ==
              config.modes.size() * (config.species.size() + 1));
  return config.modes.size();
}
exaero::EnvironmentalStateView environment(int nc, int nl, const double *temp,
                                           const double *rh) {
  require(temp && rh);
  return {exaero::View2D<const double>(temp, nc, nl),
          {},
          {},
          exaero::View2D<const double>(rh, nc, nl),
          {}};
}
void cloudIce(exaero_package_t handle, int nc, int nl, int nf,
              const double *temperature, const double *rh, const double *state,
              const double *cloud_state, int category, int preexisting_ice,
              double *ambient_out, double *cloud_out, int *ambient_status,
              int *cloud_status, bool modal) {
  auto &p = package(handle, nc, nl, nf);
  const auto np = populations(p, nf, modal);
  require(state && ambient_out && cloud_out && ambient_status && cloud_status &&
          category >= EXAERO_ICE_DUST &&
          category <= EXAERO_ICE_STRATOSPHERIC_SULFATE &&
          (preexisting_ice == 0 || preexisting_ice == 1) &&
          np <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                    EXAERO_CLOUD_ICE_PROPERTY_COUNT / nc / nl);
  const auto env = environment(nc, nl, temperature, rh);
  const exaero::View3D<const double> ambient(state, nc, nl, nf);
  std::optional<exaero::View3D<const double>> cloudborne;
  if (cloud_state)
    cloudborne.emplace(cloud_state, nc, nl, nf);
  exaero::View4D<double> a(ambient_out, nc, nl, np,
                           exaero::cloud_ice_properties::COUNT);
  exaero::View4D<double> b(cloud_out, nc, nl, np,
                           exaero::cloud_ice_properties::COUNT);
  exaero::View3D<int> as(ambient_status, nc, nl, np);
  exaero::View3D<int> bs(cloud_status, nc, nl, np);
  if (modal)
    exaero::modalCloudIce(
        p, env, ambient, cloudborne ? &*cloudborne : nullptr,
        static_cast<exaero::IceCategory>(category), preexisting_ice == 1, a, b,
        as, bs);
  else
    exaero::bulkCloudIce(
        p, env, ambient, cloudborne ? &*cloudborne : nullptr,
        static_cast<exaero::IceCategory>(category), preexisting_ice == 1, a, b,
        as, bs);
}
void fractions(exaero_package_t handle, int nc, int nl, int nf, int nss,
               const double *temperature, const double *rh,
               const double *state, const double *supersaturation,
               double *output, bool modal) {
  auto &p = package(handle, nc, nl, nf);
  const auto np = populations(p, nf, modal);
  require(nss > 0 && state && supersaturation && output &&
          np <= static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                    2 / nss / nc / nl);
  const auto env = environment(nc, nl, temperature, rh);
  const exaero::View3D<const double> ambient(state, nc, nl, nf);
  const exaero::View1D<const double> ss(supersaturation, nss);
  exaero::View4D<double> out(output, nc, nl, np, 2 * nss);
  if (modal)
    exaero::modalActivationFractions(p, env, ambient, ss, out);
  else
    exaero::bulkActivationFractions(p, env, ambient, ss, out);
}
} // namespace

extern "C" {
void exaero_bulk_cloud_ice(
    exaero_package_t pkg, int nc, int nl, int nf, const double *temperature,
    const double *rh, const double *state, const double *cloud_state,
    int category, int preexisting_ice, double *ambient_out, double *cloud_out,
    int *ambient_status, int *cloud_status, char *errmsg, int *errflg) {
  guarded(errmsg, errflg, [&] {
    cloudIce(pkg, nc, nl, nf, temperature, rh, state, cloud_state, category,
             preexisting_ice, ambient_out, cloud_out, ambient_status,
             cloud_status, false);
  });
}
void exaero_modal_cloud_ice(
    exaero_package_t pkg, int nc, int nl, int nf, const double *temperature,
    const double *rh, const double *state, const double *cloud_state,
    int category, int preexisting_ice, double *ambient_out, double *cloud_out,
    int *ambient_status, int *cloud_status, char *errmsg, int *errflg) {
  guarded(errmsg, errflg, [&] {
    cloudIce(pkg, nc, nl, nf, temperature, rh, state, cloud_state, category,
             preexisting_ice, ambient_out, cloud_out, ambient_status,
             cloud_status, true);
  });
}
void exaero_bulk_activation_fractions(
    exaero_package_t pkg, int nc, int nl, int nf, int nss,
    const double *temperature, const double *rh, const double *state,
    const double *supersaturation, double *output, char *errmsg,
    int *errflg) {
  guarded(errmsg, errflg, [&] {
    fractions(pkg, nc, nl, nf, nss, temperature, rh, state, supersaturation,
              output, false);
  });
}
void exaero_modal_activation_fractions(
    exaero_package_t pkg, int nc, int nl, int nf, int nss,
    const double *temperature, const double *rh, const double *state,
    const double *supersaturation, double *output, char *errmsg,
    int *errflg) {
  guarded(errmsg, errflg, [&] {
    fractions(pkg, nc, nl, nf, nss, temperature, rh, state, supersaturation,
              output, true);
  });
}
}
