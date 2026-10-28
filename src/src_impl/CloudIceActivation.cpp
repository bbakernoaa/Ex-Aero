/// @file CloudIceActivation.cpp
/// @brief Per-population liquid activation fractions, no state transfer.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/ModalCoupling.hpp>
#include <limits>
#include <stdexcept>
#include <vector>

namespace exaero {
namespace {
void require(bool okay, const char *why) {
  if (!okay)
    throw std::invalid_argument(std::string("FATAL ERROR: ") + why);
}
template <class V>
void shape(const V &view, std::initializer_list<std::size_t> extents) {
  require(view.data_handle() != nullptr, "null activation array");
  std::size_t axis = 0;
  for (auto extent : extents)
    require(view.extent(axis++) == extent, "activation shape mismatch");
}
template <class A, class B> void disjoint(const A &a, const B &b) {
  const auto ap = reinterpret_cast<std::uintptr_t>(a.data_handle());
  const auto bp = reinterpret_cast<std::uintptr_t>(b.data_handle());
  const auto an = a.size() * sizeof(*a.data_handle());
  const auto bn = b.size() * sizeof(*b.data_handle());
  require(ap <= std::numeric_limits<std::uintptr_t>::max() - an &&
              bp <= std::numeric_limits<std::uintptr_t>::max() - bn &&
              (ap + an <= bp || bp + bn <= ap),
          "activation arrays overlap or address overflow");
}
void environment(const EnvironmentalStateView &env, std::size_t nc,
                 std::size_t nl) {
  shape(env.temperature, {nc, nl});
  for (std::size_t i = 0; i < nc * nl; ++i)
    require(std::isfinite(env.temperature.data_handle()[i]) &&
                env.temperature.data_handle()[i] > 0,
            "invalid activation temperature");
}
void axis(const View1D<const double> &ss) {
  require(ss.data_handle() != nullptr && ss.extent(0) > 0,
          "empty activation supersaturation axis");
  for (std::size_t i = 0; i < ss.size(); ++i)
    require(std::isfinite(ss(i)) && ss(i) >= 0 && ss(i) <= 1,
            "invalid activation supersaturation");
}
void finiteOutput(const std::vector<double> &values) {
  for (double value : values)
    require(std::isfinite(value) && value >= 0 && value <= 1,
            "invalid activation fraction");
}
} // namespace

void bulkActivationFractions(IAerosolPackage &package,
                             const EnvironmentalStateView &env,
                             const View3D<const double> &ambient,
                             const View1D<const double> &supersaturation,
                             View4D<double> &fractions) {
  require(package.representation() == AerosolRepresentation::Bulk,
          "bulk package required for activation");
  auto *size = dynamic_cast<IBulkRemovalMetadata *>(&package);
  require(size != nullptr, "bulk size metadata unavailable");
  const auto nc = ambient.extent(0), nl = ambient.extent(1);
  const auto ns = package.stateFields().size();
  require(nc > 0 && nl > 0 && ns > 0, "empty bulk activation grid");
  shape(ambient, {nc, nl, ns});
  environment(env, nc, nl);
  axis(supersaturation);
  shape(fractions, {nc, nl, ns, 2 * supersaturation.extent(0)});
  disjoint(ambient, fractions);
  disjoint(env.temperature, fractions);
  disjoint(supersaturation, fractions);
  std::vector<double> values(fractions.size(), 0);
  View4D<double> result(values.data(), nc, nl, ns,
                        2 * supersaturation.extent(0));
  for (std::size_t s = 0; s < ns; ++s) {
    const auto physical = size->bulkSpeciesSize(s);
    require(std::isfinite(physical.dry_particle_diameter) &&
                physical.dry_particle_diameter > 0 &&
                std::isfinite(physical.hygroscopicity) &&
                physical.hygroscopicity >= 0,
            "invalid bulk activation material");
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const double mass = ambient(c, k, s);
        require(std::isfinite(mass) && mass >= 0,
                "invalid bulk activation mass");
        if (mass == 0 || physical.hygroscopicity == 0)
          continue;
        const double kelvin = 1.2e-9 / env.temperature(c, k);
        const double critical =
            std::sqrt(4 * std::pow(kelvin, 3) /
                      (27 * physical.hygroscopicity *
                           std::pow(physical.dry_particle_diameter, 3) +
                       1e-30));
        require(std::isfinite(critical), "invalid bulk activation threshold");
        for (std::size_t j = 0; j < supersaturation.extent(0); ++j)
          if (supersaturation(j) >= critical) {
            result(c, k, s, 2 * j) = 1;
            result(c, k, s, 2 * j + 1) = 1;
          }
      }
  }
  finiteOutput(values);
  std::copy(values.begin(), values.end(), fractions.data_handle());
}

void modalActivationFractions(IAerosolPackage &package,
                              const EnvironmentalStateView &env,
                              const View3D<const double> &ambient,
                              const View1D<const double> &supersaturation,
                              View4D<double> &fractions) {
  require(package.representation() == AerosolRepresentation::Modal,
          "modal package required for activation");
  auto *modal = dynamic_cast<IModalCoupling *>(&package);
  require(modal != nullptr, "generic modal coupling unavailable");
  const auto &config = modal->modalConfig();
  const auto nc = ambient.extent(0), nl = ambient.extent(1);
  const auto nm = config.modes.size(), ns = config.species.size();
  require(nc > 0 && nl > 0 && nm > 0 && ns > 0,
          "empty modal activation grid");
  shape(ambient, {nc, nl, nm * (ns + 1)});
  environment(env, nc, nl);
  axis(supersaturation);
  shape(fractions, {nc, nl, nm, 2 * supersaturation.extent(0)});
  disjoint(ambient, fractions);
  disjoint(env.temperature, fractions);
  disjoint(supersaturation, fractions);
  std::vector<double> properties(nc * nl * nm * modal_properties::COUNT);
  View4D<double> property_view(properties.data(), nc, nl, nm,
                               modal_properties::COUNT);
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        double mass = 0;
        for (std::size_t s = 0; s < ns; ++s) {
          const double part = ambient(c, k, m * (ns + 1) + s);
          require(std::isfinite(part) && part >= 0,
                  "invalid modal activation mass");
          if (!config.modes[m].species.empty())
            require(
                part == 0 ||
                    std::find(config.modes[m].species.begin(),
                              config.modes[m].species.end(),
                              config.species[s].name) !=
                        config.modes[m].species.end(),
                "disallowed modal activation material populated");
          mass += part;
        }
        const double number = ambient(c, k, m * (ns + 1) + ns);
        require(std::isfinite(number) && number >= 0 &&
                    std::isfinite(mass) && (mass > 0) == (number > 0),
                "invalid modal activation mass and number pair");
      }
  modal->modalProperties(env, ambient, property_view);
  std::vector<double> values(fractions.size(), 0);
  View4D<double> result(values.data(), nc, nl, nm,
                        2 * supersaturation.extent(0));
  const auto &physics = config.constants;
  require(std::isfinite(physics.surface_tension) &&
              physics.surface_tension > 0 &&
              std::isfinite(physics.water_molecular_weight) &&
              physics.water_molecular_weight > 0 &&
              std::isfinite(physics.gas_constant) && physics.gas_constant > 0,
          "invalid modal activation constants");
  for (std::size_t m = 0; m < nm; ++m) {
    const double log_sigma = std::log(config.modes[m].geometric_std_dev);
    require(std::isfinite(log_sigma) && log_sigma > 0,
            "invalid modal activation distribution");
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const double number = ambient(c, k, m * (ns + 1) + ns);
        const double radius =
            property_view(c, k, m, modal_properties::DRY_RADIUS);
        const double kappa =
            property_view(c, k, m, modal_properties::KAPPA);
        require(std::isfinite(radius) && radius >= 0 &&
                    std::isfinite(kappa) && kappa >= 0,
                "invalid modal activation properties");
        if (number == 0 || kappa == 0)
          continue;
        require(radius > 0, "nonempty modal activation radius is zero");
        const double kelvin =
            4 * physics.surface_tension * physics.water_molecular_weight /
            (physics.gas_constant * env.temperature(c, k) * 1000);
        const double critical =
            std::sqrt(4 * std::pow(kelvin, 3) /
                      (27 * kappa * std::pow(2 * radius, 3)));
        require(std::isfinite(critical) && critical > 0,
                "invalid modal activation threshold");
        for (std::size_t j = 0; j < supersaturation.extent(0); ++j) {
          const double ss = supersaturation(j);
          if (ss == 0)
            continue;
          const double x = (2.0 / 3.0) * std::log(critical / ss) /
                           (std::sqrt(2.0) * log_sigma);
          const double number_fraction = .5 * std::erfc(x);
          const double mass_fraction =
              .5 * std::erfc(x - 1.5 * std::sqrt(2.0) * log_sigma);
          result(c, k, m, 2 * j) = std::clamp(number_fraction, 0., 1.);
          result(c, k, m, 2 * j + 1) =
              std::clamp(mass_fraction, 0., 1.);
        }
      }
  }
  finiteOutput(values);
  std::copy(values.begin(), values.end(), fractions.data_handle());
}
} // namespace exaero
