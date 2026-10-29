/// @file PopulationProperties.cpp
/// @brief Staged common dry/wet properties using existing physical
/// metadata/kernels.
#include "HostCouplingInternal.hpp"
#include <exaero/BulkRemoval.hpp>
#include <exaero/ModalCoupling.hpp>
#include <exaero/PopulationProperties.hpp>
#include <numbers>
namespace exaero {
void populationProperties(IAerosolPackage &p, const EnvironmentalStateView &env,
                          const View3D<const double> &state,
                          View4D<double> &values, View4D<int> &availability,
                          View3D<int> &status) {
  using namespace host_detail;
  const auto info = hostCouplingInfo(p);
  const auto nc = state.extent(0), nl = state.extent(1),
             np = info.populations.size();
  constexpr auto count = std::size_t(PopulationProperty::Count);
  require(nc > 0 && nl > 0, "empty host property grid");
  std::size_t elements = 1;
  for (auto d : {nc, nl, np, std::max(count, std::size_t(info.field_count))}) {
    require(d > 0 && elements <= std::numeric_limits<std::size_t>::max() / d,
            "property extent overflow");
    elements *= d;
  }
  require(elements <= std::numeric_limits<std::size_t>::max() / sizeof(double),
          "property byte extent overflow");
  shape(state, {nc, nl, std::size_t(info.field_count)});
  shape(values, {nc, nl, np, count});
  shape(availability, {nc, nl, np, count});
  shape(status, {nc, nl, np});
  shape(env.temperature, {nc, nl});
  shape(env.relative_humidity, {nc, nl});
  for (const auto &input : {state}) {
    disjoint(input, values);
    disjoint(input, availability);
    disjoint(input, status);
  }
  disjoint(values, availability);
  disjoint(values, status);
  disjoint(availability, status);
  for (const auto &input : {env.temperature, env.relative_humidity}) {
    disjoint(input, values);
    disjoint(input, availability);
    disjoint(input, status);
  }
  for (std::size_t k = 0; k < nl; ++k)
    for (std::size_t c = 0; c < nc; ++c)
      require(std::isfinite(env.temperature(c, k)) &&
                  env.temperature(c, k) > 0 &&
                  std::isfinite(env.relative_humidity(c, k)) &&
                  env.relative_humidity(c, k) >= 0 &&
                  env.relative_humidity(c, k) <= 1,
              "invalid property temperature/RH");
  std::vector<double> scratch(values.size(), 0);
  std::vector<int> codes(availability.size(),
                         int(PropertyAvailability::Unavailable)),
      states(status.size());
  View4D<double> out(scratch.data(), nc, nl, np, count);
  View4D<int> available(codes.data(), nc, nl, np, count);
  View3D<int> state_status(states.data(), nc, nl, np);
  const bool modal = info.representation == AerosolRepresentation::Modal;
  const ModalConfig *config = modal ? &modalCoupling(p).modalConfig() : nullptr;
  std::vector<double> mp(nc * nl * np * modal_properties::COUNT),
      numbers(nc * nl * np), radii(np * 3), diameters(nc * nl * np * 2);
  std::vector<int> radius_codes(np * 3), populated(nc * nl * np);
  View4D<double> props(mp.data(), nc, nl, np, modal_properties::COUNT),
      diam(diameters.data(), nc, nl, np, 2);
  View3D<double> number(numbers.data(), nc, nl, np);
  View2D<double> reference(radii.data(), np, 3);
  View2D<int> reference_codes(radius_codes.data(), np, 3);
  View3D<int> populated_view(populated.data(), nc, nl, np);
  if (modal)
    modalCoupling(p).modalProperties(env, state, props);
  else {
    bulkFixedNumber(p, state, number);
    bulkReferenceRadii(p, reference, reference_codes);
    bulkScavengingSizes(p, env.relative_humidity, state, diam, populated_view);
  }
  for (std::size_t m = 0; m < np; ++m) {
    BulkSpeciesSize bulk;
    if (!modal)
      bulk = dynamic_cast<IBulkRemovalMetadata &>(p).bulkSpeciesSize(m);
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        auto set = [&](PopulationProperty q, double value,
                       PropertyAvailability code =
                           PropertyAvailability::Derived) {
          require(std::isfinite(value) && value >= 0,
                  "nonfinite or negative population property");
          out(c, k, m, int(q)) = value;
          available(c, k, m, int(q)) = int(code);
        };
        using Q = PopulationProperty;
        using A = PropertyAvailability;
        const double mass =
            modal ? props(c, k, m, modal_properties::MASS) : state(c, k, m);
        const double n =
            modal ? props(c, k, m, modal_properties::NUMBER) : number(c, k, m);
        const bool empty = mass == 0;
        state_status(c, k, m) =
            int(empty ? PopulationStatus::Empty : PopulationStatus::Populated);
        const double volume = modal ? props(c, k, m, modal_properties::VOLUME)
                                    : mass / bulk.dry_density;
        const double rg = modal ? props(c, k, m, modal_properties::DRY_RADIUS)
                                : bulk.lognormal_dg / 2;
        const double sigma =
            modal ? config->modes[m].geometric_std_dev : bulk.lognormal_sigma;
        const double ln = std::log(sigma);
        set(Q::Mass, mass);
        set(Q::Number, n);
        set(Q::DryVolume, volume);
        set(Q::DryNumberMedianRadius, rg,
            modal && !empty ? A::Derived : A::Configured);
        set(Q::DryScavengingDiameter,
            modal ? 2 * rg : bulk.dry_particle_diameter,
            modal && !empty ? A::Derived : A::Configured);
        set(Q::GeometricStdDev, sigma, A::Configured);
        const long double area = 4 * std::numbers::pi_v<long double> * n * rg *
                                 rg * std::exp(2 * ln * ln);
        require(area <= std::numeric_limits<double>::max(),
                "dry surface area overflow");
        set(Q::DrySurfaceArea, double(area));
        if (modal) {
          const auto &mode = config->modes[m];
          if (mode.dry_geometric_radius_min)
            set(Q::LowerDryRadius, *mode.dry_geometric_radius_min,
                A::Configured);
          if (mode.dry_geometric_radius_max)
            set(Q::UpperDryRadius, *mode.dry_geometric_radius_max,
                A::Configured);
          set(Q::EffectiveDryRadius, rg * std::exp(2.5 * ln * ln));
          if (!empty) {
            set(Q::Kappa, props(c, k, m, modal_properties::KAPPA));
            set(Q::DryDensity, mass / volume);
            set(Q::WetDensity, props(c, k, m, modal_properties::WET_DENSITY));
            const auto wet = props(c, k, m, modal_properties::WET_RADIUS);
            set(Q::WetNumberMedianRadius, wet);
            set(Q::WetScavengingDiameter, 2 * wet);
          }
        } else {
          set(Q::Kappa, bulk.hygroscopicity, A::Configured);
          set(Q::DryDensity, bulk.dry_density, A::Configured);
          const double rh = std::min(env.relative_humidity(c, k), .99),
                       growth_volume = bulk.hygroscopicity * rh / (1 - rh);
          const double growth = std::cbrt(1 + growth_volume);
          set(Q::WetNumberMedianRadius, rg * growth);
          set(Q::WetScavengingDiameter,
              empty ? bulk.dry_particle_diameter * growth : diam(c, k, m, 1));
          set(Q::WetDensity,
              (bulk.dry_density + 1000 * growth_volume) / (1 + growth_volume));
          if (reference_codes(m, 0))
            set(Q::LowerDryRadius, reference(m, 0), A::Configured);
          if (reference_codes(m, 1))
            set(Q::UpperDryRadius, reference(m, 1), A::Configured);
          set(Q::EffectiveDryRadius, reference(m, 2),
              static_cast<A>(reference_codes(m, 2)));
        }
      }
  }
  std::copy(scratch.begin(), scratch.end(), values.data_handle());
  std::copy(codes.begin(), codes.end(), availability.data_handle());
  std::copy(states.begin(), states.end(), status.data_handle());
}
} // namespace exaero
