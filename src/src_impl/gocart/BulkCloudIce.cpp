/// @file BulkCloudIce.cpp
/// @brief Transactional cloud/ice candidate diagnostics for mass-only bulk.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceMath.hpp>
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
  require(view.data_handle() != nullptr, "null cloud/ice array");
  std::size_t axis = 0;
  for (auto extent : extents)
    require(view.extent(axis++) == extent, "cloud/ice shape mismatch");
}
template <class A, class B> void disjoint(const A &a, const B &b) {
  const auto ap = reinterpret_cast<std::uintptr_t>(a.data_handle());
  const auto bp = reinterpret_cast<std::uintptr_t>(b.data_handle());
  const auto an = a.size() * sizeof(*a.data_handle());
  const auto bn = b.size() * sizeof(*b.data_handle());
  require(ap <= std::numeric_limits<std::uintptr_t>::max() - an &&
              bp <= std::numeric_limits<std::uintptr_t>::max() - bn &&
              (ap + an <= bp || bp + bn <= ap),
          "cloud/ice arrays overlap or address overflow");
}
void fill(IAerosolPackage &package, const EnvironmentalStateView &env,
          const View3D<const double> &state, IceCategory category,
          bool preexisting_ice, IBulkCloudIceMetadata &ice,
          const View2D<double> &radii, View4D<double> &out,
          View3D<int> &status) {
  const auto nc = state.extent(0), nl = state.extent(1), ns = state.extent(2);
  std::vector<double> numbers(nc * nl * ns);
  View3D<double> number(numbers.data(), nc, nl, ns);
  bulkFixedNumber(package, state, number);
  std::vector<double> sizes(nc * nl * ns * 2);
  std::vector<int> size_codes(nc * nl * ns);
  View4D<double> diameters(sizes.data(), nc, nl, ns, 2);
  View3D<int> size_status(size_codes.data(), nc, nl, ns);
  bulkScavengingSizes(package, env.relative_humidity, state, diameters,
                      size_status);
  for (std::size_t s = 0; s < ns; ++s) {
    const auto physical = dynamic_cast<IBulkRemovalMetadata &>(package)
                              .bulkSpeciesSize(s);
    const auto tags = ice.bulkIceCategories(s);
    const auto rules = ice.bulkIceRules(s);
    validateIceConfiguration(tags, rules);
    const bool matching =
        std::find(tags.begin(), tags.end(), category) != tags.end();
    const auto rule = std::find_if(rules.begin(), rules.end(),
                                   [&](const auto &v) {
                                     return v.category == category;
                                   });
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const double mass = state(c, k, s);
        out(c, k, s, cloud_ice_properties::MASS) = mass;
        out(c, k, s, cloud_ice_properties::NUMBER) = number(c, k, s);
        out(c, k, s, cloud_ice_properties::DRY_NUMBER_MEDIAN_RADIUS) =
            physical.lognormal_dg / 2;
        out(c, k, s, cloud_ice_properties::DRY_SCAVENGING_RADIUS) =
            diameters(c, k, s, 0) / 2;
        out(c, k, s, cloud_ice_properties::WET_SCAVENGING_RADIUS) =
            diameters(c, k, s, 1) / 2;
        out(c, k, s, cloud_ice_properties::LOWER_RADIUS) = radii(s, 0);
        out(c, k, s, cloud_ice_properties::UPPER_RADIUS) = radii(s, 1);
        out(c, k, s, cloud_ice_properties::EFFECTIVE_RADIUS) = radii(s, 2);
        out(c, k, s, cloud_ice_properties::TYPE_WEIGHT) = matching ? 1 : 0;
        if (mass == 0) {
          status(c, k, s) = static_cast<int>(CloudIceStatus::Empty);
          continue;
        }
        if (rule == rules.end()) {
          status(c, k, s) = static_cast<int>(CloudIceStatus::NotConfigured);
          continue;
        }
        const auto &size_rule = preexisting_ice && rule->preexisting_ice
                                    ? *rule->preexisting_ice
                                    : rule->normal;
        const double weight =
            iceSizeWeight(size_rule, physical.lognormal_dg,
                          physical.lognormal_sigma);
        out(c, k, s, cloud_ice_properties::SIZE_WEIGHT) = weight;
        out(c, k, s, cloud_ice_properties::ELIGIBLE_NUMBER) =
            number(c, k, s) * weight;
        out(c, k, s, cloud_ice_properties::HETEROGENEOUS_FREEZING) =
            rule->heterogeneous_freezing ? 1 : 0;
        status(c, k, s) = static_cast<int>(CloudIceStatus::Configured);
      }
  }
  for (std::size_t i = 0; i < out.size(); ++i)
    require(std::isfinite(out.data_handle()[i]),
            "nonfinite bulk cloud/ice output");
}
} // namespace

void bulkCloudIce(IAerosolPackage &package, const EnvironmentalStateView &env,
                  const View3D<const double> &ambient,
                  const View3D<const double> *cloudborne,
                  IceCategory category, bool preexisting_ice,
                  View4D<double> &ambient_out, View4D<double> &cloudborne_out,
                  View3D<int> &ambient_status,
                  View3D<int> &cloudborne_status) {
  require(package.representation() == AerosolRepresentation::Bulk,
          "bulk package required");
  auto *ice = dynamic_cast<IBulkCloudIceMetadata *>(&package);
  require(ice != nullptr && validIceCategory(category),
          "bulk cloud/ice metadata or category invalid");
  const auto fields = package.stateFields();
  const auto nc = ambient.extent(0), nl = ambient.extent(1), ns = fields.size();
  require(nc > 0 && nl > 0 && ns > 0, "empty bulk cloud/ice grid");
  shape(ambient, {nc, nl, ns});
  if (cloudborne)
    shape(*cloudborne, {nc, nl, ns});
  shape(env.relative_humidity, {nc, nl});
  shape(ambient_out, {nc, nl, ns, cloud_ice_properties::COUNT});
  shape(cloudborne_out, {nc, nl, ns, cloud_ice_properties::COUNT});
  shape(ambient_status, {nc, nl, ns});
  shape(cloudborne_status, {nc, nl, ns});
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
  std::vector<double> radius_values(ns * 3);
  std::vector<int> radius_codes(ns * 3);
  View2D<double> radii(radius_values.data(), ns, 3);
  View2D<int> radius_status(radius_codes.data(), ns, 3);
  bulkReferenceRadii(package, radii, radius_status);
  std::vector<double> ambient_values(ambient_out.size(), 0);
  std::vector<double> cloud_values(cloudborne_out.size(), 0);
  std::vector<int> ambient_codes(ambient_status.size());
  std::vector<int> cloud_codes(cloudborne_status.size(),
                               static_cast<int>(CloudIceStatus::Absent));
  View4D<double> a(ambient_values.data(), nc, nl, ns,
                   cloud_ice_properties::COUNT);
  View4D<double> b(cloud_values.data(), nc, nl, ns,
                   cloud_ice_properties::COUNT);
  View3D<int> a_status(ambient_codes.data(), nc, nl, ns);
  View3D<int> b_status(cloud_codes.data(), nc, nl, ns);
  fill(package, env, ambient, category, preexisting_ice, *ice, radii, a,
       a_status);
  if (cloudborne)
    fill(package, env, *cloudborne, category, preexisting_ice, *ice, radii, b,
         b_status);
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
