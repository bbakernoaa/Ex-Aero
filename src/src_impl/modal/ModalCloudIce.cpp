/// @file ModalCloudIce.cpp
/// @brief Transactional cloud/ice diagnostics for generic modal aerosol.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceMath.hpp>
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
  require(view.data_handle() != nullptr, "null modal cloud/ice array");
  std::size_t axis = 0;
  for (auto extent : extents)
    require(view.extent(axis++) == extent, "modal cloud/ice shape mismatch");
}
template <class A, class B> void disjoint(const A &a, const B &b) {
  const auto ap = reinterpret_cast<std::uintptr_t>(a.data_handle());
  const auto bp = reinterpret_cast<std::uintptr_t>(b.data_handle());
  const auto an = a.size() * sizeof(*a.data_handle());
  const auto bn = b.size() * sizeof(*b.data_handle());
  require(ap <= std::numeric_limits<std::uintptr_t>::max() - an &&
              bp <= std::numeric_limits<std::uintptr_t>::max() - bn &&
              (ap + an <= bp || bp + bn <= ap),
          "modal cloud/ice arrays overlap or address overflow");
}
void validateState(const ModalConfig &config,
                   const View3D<const double> &state) {
  const auto ns = config.species.size();
  for (std::size_t m = 0; m < config.modes.size(); ++m)
    for (std::size_t k = 0; k < state.extent(1); ++k)
      for (std::size_t c = 0; c < state.extent(0); ++c) {
        double mass = 0;
        for (std::size_t s = 0; s < ns; ++s) {
          const double value = state(c, k, m * (ns + 1) + s);
          require(std::isfinite(value) && value >= 0,
                  "invalid modal material mass");
          if (!config.modes[m].species.empty())
            require(
                value == 0 ||
                    std::find(config.modes[m].species.begin(),
                              config.modes[m].species.end(),
                              config.species[s].name) !=
                        config.modes[m].species.end(),
                "disallowed modal material populated");
          mass += value;
        }
        const double number = state(c, k, m * (ns + 1) + ns);
        require(std::isfinite(number) && number >= 0 &&
                    std::isfinite(mass) && (mass > 0) == (number > 0),
                "invalid modal mass and number pair");
      }
}
void fill(IModalCoupling &modal, const EnvironmentalStateView &env,
          const View3D<const double> &state, IceCategory category,
          bool preexisting_ice, View4D<double> &out, View3D<int> &status) {
  const auto &config = modal.modalConfig();
  const auto nc = state.extent(0), nl = state.extent(1);
  const auto nm = config.modes.size(), ns = config.species.size();
  validateState(config, state);
  std::vector<double> properties(nc * nl * nm * modal_properties::COUNT);
  View4D<double> property_view(properties.data(), nc, nl, nm,
                               modal_properties::COUNT);
  modal.modalProperties(env, state, property_view);
  for (std::size_t m = 0; m < nm; ++m) {
    const auto &mode = config.modes[m];
    std::vector<IceCategory> available;
    for (std::size_t s = 0; s < ns; ++s)
      if (mode.species.empty() ||
          std::find(mode.species.begin(), mode.species.end(),
                    config.species[s].name) != mode.species.end())
        for (auto tag : config.species[s].ice_categories)
          if (std::find(available.begin(), available.end(), tag) ==
              available.end())
            available.push_back(tag);
    validateIceConfiguration(available, mode.ice_rules);
    const auto rule = std::find_if(mode.ice_rules.begin(), mode.ice_rules.end(),
                                   [&](const auto &r) {
                                     return r.category == category;
                                   });
    const double log_sigma = std::log(mode.geometric_std_dev);
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        double mass = 0, tagged_mass = 0;
        for (std::size_t s = 0; s < ns; ++s) {
          const double part = state(c, k, m * (ns + 1) + s);
          mass += part;
          if (std::find(config.species[s].ice_categories.begin(),
                        config.species[s].ice_categories.end(),
                        category) != config.species[s].ice_categories.end())
            tagged_mass += part;
        }
        const double number = state(c, k, m * (ns + 1) + ns);
        const double dry =
            property_view(c, k, m, modal_properties::DRY_RADIUS);
        const double wet =
            property_view(c, k, m, modal_properties::WET_RADIUS);
        require(std::isfinite(dry) && dry >= 0 && std::isfinite(wet) &&
                    wet >= 0,
                "invalid modal property radius");
        out(c, k, m, cloud_ice_properties::MASS) = mass;
        out(c, k, m, cloud_ice_properties::NUMBER) = number;
        out(c, k, m, cloud_ice_properties::DRY_NUMBER_MEDIAN_RADIUS) = dry;
        out(c, k, m, cloud_ice_properties::DRY_SCAVENGING_RADIUS) = dry;
        out(c, k, m, cloud_ice_properties::WET_SCAVENGING_RADIUS) = wet;
        out(c, k, m, cloud_ice_properties::LOWER_RADIUS) =
            mode.dry_geometric_radius_min.value_or(0);
        out(c, k, m, cloud_ice_properties::UPPER_RADIUS) =
            mode.dry_geometric_radius_max.value_or(0);
        out(c, k, m, cloud_ice_properties::EFFECTIVE_RADIUS) =
            dry * std::exp(2.5 * log_sigma * log_sigma);
        out(c, k, m, cloud_ice_properties::TYPE_WEIGHT) =
            mass > 0 ? tagged_mass / mass : 0;
        if (mass == 0) {
          status(c, k, m) = static_cast<int>(CloudIceStatus::Empty);
          continue;
        }
        if (rule == mode.ice_rules.end()) {
          status(c, k, m) = static_cast<int>(CloudIceStatus::NotConfigured);
          continue;
        }
        const auto &size_rule = preexisting_ice && rule->preexisting_ice
                                    ? *rule->preexisting_ice
                                    : rule->normal;
        const double weight = iceSizeWeight(size_rule, 2 * dry,
                                            mode.geometric_std_dev);
        out(c, k, m, cloud_ice_properties::SIZE_WEIGHT) = weight;
        out(c, k, m, cloud_ice_properties::ELIGIBLE_NUMBER) =
            number * out(c, k, m, cloud_ice_properties::TYPE_WEIGHT) * weight;
        out(c, k, m, cloud_ice_properties::HETEROGENEOUS_FREEZING) =
            rule->heterogeneous_freezing ? 1 : 0;
        status(c, k, m) = static_cast<int>(CloudIceStatus::Configured);
      }
  }
  for (std::size_t i = 0; i < out.size(); ++i)
    require(std::isfinite(out.data_handle()[i]),
            "nonfinite modal cloud/ice output");
}
} // namespace

void modalCloudIce(IAerosolPackage &package, const EnvironmentalStateView &env,
                   const View3D<const double> &ambient,
                   const View3D<const double> *cloudborne,
                   IceCategory category, bool preexisting_ice,
                   View4D<double> &ambient_out, View4D<double> &cloudborne_out,
                   View3D<int> &ambient_status,
                   View3D<int> &cloudborne_status) {
  require(package.representation() == AerosolRepresentation::Modal,
          "modal package required");
  auto *modal = dynamic_cast<IModalCoupling *>(&package);
  require(modal != nullptr && validIceCategory(category),
          "generic modal coupling or category invalid");
  const auto &config = modal->modalConfig();
  const auto nc = ambient.extent(0), nl = ambient.extent(1);
  const auto nm = config.modes.size(), ns = config.species.size();
  require(nc > 0 && nl > 0 && nm > 0 && ns > 0,
          "empty modal cloud/ice grid");
  shape(ambient, {nc, nl, nm * (ns + 1)});
  if (cloudborne)
    shape(*cloudborne, {nc, nl, nm * (ns + 1)});
  shape(env.temperature, {nc, nl});
  shape(env.relative_humidity, {nc, nl});
  shape(ambient_out, {nc, nl, nm, cloud_ice_properties::COUNT});
  shape(cloudborne_out, {nc, nl, nm, cloud_ice_properties::COUNT});
  shape(ambient_status, {nc, nl, nm});
  shape(cloudborne_status, {nc, nl, nm});
  disjoint(ambient, ambient_out);
  disjoint(ambient, cloudborne_out);
  disjoint(ambient, ambient_status);
  disjoint(ambient, cloudborne_status);
  if (cloudborne) {
    disjoint(*cloudborne, ambient_out);
    disjoint(*cloudborne, cloudborne_out);
    disjoint(*cloudborne, ambient_status);
    disjoint(*cloudborne, cloudborne_status);
  }
  disjoint(ambient_out, cloudborne_out);
  disjoint(ambient_out, ambient_status);
  disjoint(ambient_out, cloudborne_status);
  disjoint(cloudborne_out, ambient_status);
  disjoint(cloudborne_out, cloudborne_status);
  disjoint(ambient_status, cloudborne_status);
  for (std::size_t i = 0; i < env.temperature.size(); ++i)
    require(std::isfinite(env.temperature.data_handle()[i]) &&
                env.temperature.data_handle()[i] > 0 &&
                std::isfinite(env.relative_humidity.data_handle()[i]) &&
                env.relative_humidity.data_handle()[i] >= 0 &&
                env.relative_humidity.data_handle()[i] <= 1,
            "invalid modal cloud environment");
  std::vector<double> ambient_values(ambient_out.size(), 0);
  std::vector<double> cloud_values(cloudborne_out.size(), 0);
  std::vector<int> ambient_codes(ambient_status.size());
  std::vector<int> cloud_codes(cloudborne_status.size(),
                               static_cast<int>(CloudIceStatus::Absent));
  View4D<double> a(ambient_values.data(), nc, nl, nm,
                   cloud_ice_properties::COUNT);
  View4D<double> b(cloud_values.data(), nc, nl, nm,
                   cloud_ice_properties::COUNT);
  View3D<int> a_status(ambient_codes.data(), nc, nl, nm);
  View3D<int> b_status(cloud_codes.data(), nc, nl, nm);
  fill(*modal, env, ambient, category, preexisting_ice, a, a_status);
  if (cloudborne)
    fill(*modal, env, *cloudborne, category, preexisting_ice, b, b_status);
  std::copy(ambient_values.begin(), ambient_values.end(),
            ambient_out.data_handle());
  std::copy(cloud_values.begin(), cloud_values.end(),
            cloudborne_out.data_handle());
  std::copy(ambient_codes.begin(), ambient_codes.end(),
            ambient_status.data_handle());
  std::copy(cloud_codes.begin(), cloud_codes.end(),
            cloudborne_status.data_handle());
}
} // namespace exaero
