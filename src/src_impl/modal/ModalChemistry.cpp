/// @file ModalChemistry.cpp
/// @brief Validated, additive host chemistry interfaces; no integration or I/O.
/// @author EX-aero contributors
/// @details History: 2026-10. Host layout_left; Kokkos kernels are isolated.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exaero/ModalChemistry.hpp>
#include <limits>
#include <modal/ModalChemistrySolver.hpp>
#include <modal/ModalPackage.hpp>
#include <set>

namespace exaero {
namespace {
/// @brief Fail before mutating any output.
void require(bool value, const char *message) {
  if (!value)
    throw std::invalid_argument(std::string("FATAL ERROR: modal chemistry: ") +
                                message);
}
/// @brief Checked concrete boundary without changing public package vtables.
ModalPackage &checkedPackage(IAerosolPackage &package) {
  auto *p = dynamic_cast<ModalPackage *>(&package);
  require(p, "generic modal package required");
  return *p;
}
/// @brief Validate extents, nonnull data and bounded product before allocation.
template <class View, class... Sizes>
void shape(const View &view, Sizes... sizes) {
  const std::size_t expected[] = {static_cast<std::size_t>(sizes)...};
  require(view.data_handle(), "null buffer");
  std::size_t product = 1;
  for (std::size_t d = 0; d < sizeof...(sizes); ++d) {
    require(expected[d] > 0 && view.extent(d) == expected[d],
            "buffer shape mismatch");
    require(product <=
                static_cast<std::size_t>(std::numeric_limits<int>::max()) /
                    expected[d],
            "workspace exceeds index range");
    product *= expected[d];
  }
}
/// @brief Host array interval for overlap checks.
struct Buffer {
  const double *data;
  std::size_t size;
};
/// @brief Construct a validated-view interval.
template <class View> Buffer buffer(const View &view) {
  return {view.data_handle(), view.size()};
}
/// @brief Reject outputs that alias input state or another output.
void distinct(std::vector<Buffer> outputs, const std::vector<Buffer> &inputs) {
  const auto overlaps = [](const Buffer &a, const Buffer &b) {
    const auto x = reinterpret_cast<std::uintptr_t>(a.data),
               y = reinterpret_cast<std::uintptr_t>(b.data);
    return x <= y ? (y - x) / sizeof(double) < a.size
                  : (x - y) / sizeof(double) < b.size;
  };
  for (std::size_t i = 0; i < outputs.size(); ++i) {
    for (const auto &input : inputs)
      require(!overlaps(outputs[i], input), "outputs must not alias inputs");
    for (std::size_t j = 0; j < i; ++j)
      require(!overlaps(outputs[i], outputs[j]), "outputs must not overlap");
  }
}
/// @brief Validate signed/nonnegative scalar arrays before launching kernels.
template <class View> void finite(const View &view, bool nonnegative = false) {
  for (std::size_t i = 0; i < view.size(); ++i)
    require(std::isfinite(view.data_handle()[i]) &&
                (!nonnegative || view.data_handle()[i] >= 0),
            "nonfinite or negative input");
}
/// @brief Resolve a species name in a validated configuration.
int species(const ModalConfig &c, const std::string &name) {
  const auto it = std::find_if(c.species.begin(), c.species.end(),
                               [&](const auto &s) { return s.name == name; });
  require(it != c.species.end(), "unknown material binding");
  return it - c.species.begin();
}
/// @brief Resolve a mode name in a validated configuration.
int mode(const ModalConfig &c, const std::string &name) {
  const auto it = std::find_if(c.modes.begin(), c.modes.end(),
                               [&](const auto &m) { return m.name == name; });
  require(it != c.modes.end(), "unknown aging mode");
  return it - c.modes.begin();
}
} // namespace

void modalPhaseExchange(IAerosolPackage &package,
                        const EnvironmentalStateView &env,
                        const View3D<const double> &state,
                        std::span<const ModalCondensable> bindings,
                        const View3D<const double> &gas,
                        const View4D<const double> &equilibrium,
                        View3D<double> &aerosol_rates,
                        View3D<double> &gas_rates, View4D<double> *jac,
                        const View4D<const double> *equilibrium_jac) {
  auto &p = checkedPackage(package);
  p.validateChemistry(env, state);
  const auto &config = p.modalConfig();
  const auto nc = state.extent(0), nl = state.extent(1), nf = state.extent(2),
             nm = config.modes.size(), ns = config.species.size(),
             nv = bindings.size(), nq = nf + nv;
  require(nv > 0 && nv <= ns, "bind one distinct vapor per dry material");
  shape(gas, nc, nl, nv);
  shape(equilibrium, nc, nl, nm, nv);
  shape(aerosol_rates, nc, nl, nf);
  shape(gas_rates, nc, nl, nv);
  std::vector<Buffer> inputs{buffer(state), buffer(gas), buffer(equilibrium),
                             buffer(env.temperature),
                             buffer(env.relative_humidity)};
  std::vector<Buffer> outputs{buffer(aerosol_rates), buffer(gas_rates)};
  const auto ni = jac ? jac->extent(3) : nq;
  require(ni >= nq,
          "exchange Jacobian needs aerosol and bound gas columns first");
  if (jac) {
    shape(*jac, nc, nl, nq, ni);
    outputs.push_back(buffer(*jac));
  }
  require(!equilibrium_jac || jac,
          "equilibrium derivatives require a Jacobian output");
  if (equilibrium_jac) {
    shape(*equilibrium_jac, nc, nl, nm * nv, ni);
    finite(*equilibrium_jac);
    inputs.push_back(buffer(*equilibrium_jac));
  }
  distinct(outputs, inputs);
  finite(gas, true);
  finite(equilibrium, true);
  std::vector<int> resolved, allowed(nm * ns);
  std::set<int> seen;
  std::vector<ModalVapor> vapors;
  for (const auto &v : bindings) {
    const int s = species(config, v.species);
    require(seen.insert(s).second, "duplicate vapor/material binding");
    resolved.push_back(s);
    vapors.push_back(v.transport);
  }
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t s = 0; s < ns; ++s)
      allowed[m * ns + s] =
          config.modes[m].species.empty() ||
          std::find(config.modes[m].species.begin(),
                    config.modes[m].species.end(),
                    config.species[s].name) != config.modes[m].species.end();
  std::vector<double> uptake(nc * nl * nm * nv),
      derivative(jac ? uptake.size() * nf : 0);
  View4D<double> ku(uptake.data(), nc, nl, nm, nv),
      dk(derivative.data(), nc, nl, nm, jac ? nv * nf : 0);
  p.modalUptake(env, state, vapors, ku, jac ? &dk : nullptr);
  modal::exchange(nc, nl, ns, nm, resolved, allowed, ni, state.data_handle(),
                  gas.data_handle(), equilibrium.data_handle(), uptake.data(),
                  jac ? derivative.data() : nullptr,
                  equilibrium_jac ? equilibrium_jac->data_handle() : nullptr,
                  aerosol_rates.data_handle(), gas_rates.data_handle(),
                  jac ? jac->data_handle() : nullptr);
}

void modalAgingTerms(IAerosolPackage &package,
                     const EnvironmentalStateView &env,
                     const View3D<const double> &state,
                     const View3D<const double> &condensation,
                     const View3D<const double> &coagulation,
                     View3D<double> &rates, View4D<double> *jac,
                     const View4D<const double> *condensation_jac,
                     const View4D<const double> *coagulation_jac) {
  auto &p = checkedPackage(package);
  p.validateChemistry(env, state);
  const auto &config = p.modalConfig();
  const auto nc = state.extent(0), nl = state.extent(1), nf = state.extent(2),
             ns = config.species.size(), nm = config.modes.size();
  const auto ni = jac ? jac->extent(3) : nf;
  require(ni >= nf, "aging Jacobian needs all aerosol columns first");
  shape(condensation, nc, nl, nf);
  shape(coagulation, nc, nl, nf);
  shape(rates, nc, nl, nf);
  finite(condensation);
  finite(coagulation, true);
  std::vector<Buffer> inputs{buffer(state), buffer(condensation),
                             buffer(coagulation), buffer(env.temperature),
                             buffer(env.relative_humidity)};
  std::vector<Buffer> outputs{buffer(rates)};
  if (jac) {
    shape(*jac, nc, nl, nf, ni);
    outputs.push_back(buffer(*jac));
  }
  require(jac || (!condensation_jac && !coagulation_jac),
          "budget derivatives require a Jacobian output");
  for (const auto *d : {condensation_jac, coagulation_jac})
    if (d) {
      require(d->extent(2) >= nf && d->extent(2) <= ni,
              "budget Jacobian rows must be between nf and ni");
      shape(*d, nc, nl, d->extent(2), ni);
      finite(*d);
      inputs.push_back(buffer(*d));
    }
  distinct(outputs, inputs);
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c)
        require(condensation(c, k, m * (ns + 1) + ns) == 0 &&
                    coagulation(c, k, m * (ns + 1) + ns) == 0,
                "coating budgets contain mass only, not particle number");
  modal::AgingParameters params;
  for (const auto &s : config.species)
    params.density.push_back(s.density);
  for (const auto &m : config.modes)
    params.sigma.push_back(m.geometric_std_dev);
  for (const auto &route : config.aging_routes) {
    params.source.push_back(mode(config, route.source));
    params.destination.push_back(mode(config, route.destination));
    params.thickness.push_back(route.coating_thickness);
    for (const auto &s : config.species)
      params.coating.push_back(std::find(route.coating_species.begin(),
                                         route.coating_species.end(), s.name) !=
                               route.coating_species.end());
  }
  modal::aging(nc, nl, params, ni, state.data_handle(),
               condensation.data_handle(), coagulation.data_handle(),
               condensation_jac ? condensation_jac->data_handle() : nullptr,
               condensation_jac ? condensation_jac->extent(2) : 0,
               coagulation_jac ? coagulation_jac->data_handle() : nullptr,
               coagulation_jac ? coagulation_jac->extent(2) : 0,
               rates.data_handle(), jac ? jac->data_handle() : nullptr);
}

void advanceModalWithCoating(IAerosolPackage &package,
                             const EnvironmentalStateView &env,
                             View3D<double> &state, double dt,
                             View3D<double> &gain) {
  const auto &config = modalCoupling(package).modalConfig();
  require(config.transfer_routes.empty() ||
              config.transfer_owner == ModalTransferOwner::Host,
          "coating-budget steps require host-owned transfer routes");
  shape(gain, state.extent(0), state.extent(1), state.extent(2));
  distinct({buffer(state), buffer(gain)},
           {buffer(env.temperature), buffer(env.pressure),
            buffer(env.air_density), buffer(env.relative_humidity),
            buffer(env.layer_thickness)});
  checkedPackage(package).advanceWithCoating(env, state, dt, &gain);
}
void renameModalModes(IAerosolPackage &package,
                      const EnvironmentalStateView &env,
                      View3D<double> &state) {
  distinct({buffer(state)},
           {buffer(env.temperature), buffer(env.pressure),
            buffer(env.air_density), buffer(env.relative_humidity),
            buffer(env.layer_thickness)});
  checkedPackage(package).renameModes(env, state);
}
} // namespace exaero
