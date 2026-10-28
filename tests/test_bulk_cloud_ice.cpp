#ifdef NDEBUG
#undef NDEBUG
#endif
#include <array>
#include <cassert>
#include <cmath>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceCoupling_C.h>
#include <stdexcept>

using namespace exaero;
namespace {
struct Package final : IAerosolPackage, IBulkRemovalMetadata,
                       IBulkCloudIceMetadata {
  bool supply_rule = true;
  bool hydrophilic = true;
  void initialize(const std::string &) override {}
  void executeMicrophysics(const EnvironmentalStateView &, View3D<double> &,
                           double) override {}
  void computeDerivedDiagnostics(const EnvironmentalStateView &,
                                 const View3D<const double> &,
                                 View3D<double> &) override {}
  void computeOptics(const EnvironmentalStateView &,
                     const View3D<const double> &, const View1D<const double> &,
                     View4D<double> &) override {}
  void computeCCN(const EnvironmentalStateView &, const View3D<const double> &,
                  const View1D<const double> &, View4D<double> &) override {}
  void computeEmissions(const EnvironmentalStateView &,
                        const EmissionsInputView &, View3D<double> &) override {}
  int getSpeciesIndex(const std::string &) const override { return 0; }
  std::string getSpeciesName(int) const override { return "DU"; }
  std::vector<StateField> stateFields() const override {
    return {{"DU", StateQuantity::Mass, -1, false}};
  }
  BulkSpeciesSize bulkSpeciesSize(std::size_t) const override {
    return {8e-7, 1000, 4e-7, 1.5, 1e-7, 9e-7, 5e-7,
            hydrophilic ? .5 : 0};
  }
  std::vector<IceCategory> bulkIceCategories(std::size_t) const override {
    return {IceCategory::Dust};
  }
  std::vector<IcePopulationRule> bulkIceRules(std::size_t) const override {
    if (!supply_rule)
      return {};
    return {{IceCategory::Dust,
             {IceSizeRuleKind::FixedFraction, .25},
             std::nullopt, true}};
  }
};
bool near(double a, double b) {
  return std::abs(a - b) <= 1e-12 * std::max({1e-30, std::abs(a), std::abs(b)});
}
} // namespace

int main() {
  Package p;
  std::array<double, 1> temperature{280}, pressure{90000}, density{1},
      humidity{.5}, thickness{100}, ambient_mass{1e-12}, cloud_mass{5e-13};
  EnvironmentalStateView env{View2D<const double>(temperature.data(), 1, 1),
                             View2D<const double>(pressure.data(), 1, 1),
                             View2D<const double>(density.data(), 1, 1),
                             View2D<const double>(humidity.data(), 1, 1),
                             View2D<const double>(thickness.data(), 1, 1)};
  View3D<const double> ambient(ambient_mass.data(), 1, 1, 1);
  View3D<const double> cloud(cloud_mass.data(), 1, 1, 1);
  std::array<double, cloud_ice_properties::COUNT> a{}, c{};
  std::array<int, 1> as{}, cs{};
  View4D<double> ambient_out(a.data(), 1, 1, 1,
                             cloud_ice_properties::COUNT);
  View4D<double> cloud_out(c.data(), 1, 1, 1, cloud_ice_properties::COUNT);
  View3D<int> ambient_status(as.data(), 1, 1, 1);
  View3D<int> cloud_status(cs.data(), 1, 1, 1);
  bulkCloudIce(p, env, ambient, &cloud, IceCategory::Dust, false, ambient_out,
               cloud_out, ambient_status, cloud_status);
  const double l = std::log(1.5);
  const double particle_mass =
      std::acos(-1.) / 6 * 1000 * std::pow(4e-7, 3) * std::exp(4.5 * l * l);
  assert(near(a[cloud_ice_properties::NUMBER], 1e-12 / particle_mass));
  assert(near(c[cloud_ice_properties::NUMBER], 5e-13 / particle_mass));
  assert(near(a[cloud_ice_properties::DRY_NUMBER_MEDIAN_RADIUS], 2e-7));
  assert(near(a[cloud_ice_properties::DRY_SCAVENGING_RADIUS], 4e-7));
  assert(near(a[cloud_ice_properties::WET_SCAVENGING_RADIUS],
              4e-7 * std::cbrt(1.5)));
  assert(a[cloud_ice_properties::TYPE_WEIGHT] == 1);
  assert(a[cloud_ice_properties::SIZE_WEIGHT] == .25);
  assert(near(a[cloud_ice_properties::ELIGIBLE_NUMBER],
              .25e-12 / particle_mass));
  assert(a[cloud_ice_properties::HETEROGENEOUS_FREEZING] == 1);
  assert(as[0] == static_cast<int>(CloudIceStatus::Configured));
  assert(cs[0] == static_cast<int>(CloudIceStatus::Configured));

  bulkCloudIce(p, env, ambient, nullptr, IceCategory::Sulfate, false,
               ambient_out, cloud_out, ambient_status, cloud_status);
  assert(a[cloud_ice_properties::TYPE_WEIGHT] == 0);
  assert(a[cloud_ice_properties::ELIGIBLE_NUMBER] == 0);
  assert(as[0] == static_cast<int>(CloudIceStatus::NotConfigured));
  assert(cs[0] == static_cast<int>(CloudIceStatus::Absent));
  for (double v : c)
    assert(v == 0);

  p.supply_rule = false;
  bulkCloudIce(p, env, ambient, nullptr, IceCategory::Dust, false,
               ambient_out, cloud_out, ambient_status, cloud_status);
  assert(a[cloud_ice_properties::TYPE_WEIGHT] == 1);
  assert(a[cloud_ice_properties::SIZE_WEIGHT] == 0);
  assert(as[0] == static_cast<int>(CloudIceStatus::NotConfigured));
  p.supply_rule = true;

  std::array<double, 2> supersaturation{1e-10, 1e-8};
  std::array<double, 4> fractions{};
  View1D<const double> ss(supersaturation.data(), 2);
  View4D<double> activation(fractions.data(), 1, 1, 1, 4);
  bulkActivationFractions(p, env, ambient, ss, activation);
  assert((fractions == std::array<double, 4>{0, 0, 1, 1}));
  p.hydrophilic = false;
  bulkActivationFractions(p, env, ambient, ss, activation);
  assert((fractions == std::array<double, 4>{0, 0, 0, 0}));
  p.hydrophilic = true;

  std::array<char, 256> error{};
  int error_flag = -1;
  exaero_bulk_cloud_ice(&p, 1, 1, 1, temperature.data(), humidity.data(),
                        ambient_mass.data(), cloud_mass.data(),
                        EXAERO_ICE_DUST, 0, a.data(), c.data(), as.data(),
                        cs.data(), error.data(), &error_flag);
  assert(error_flag == 0 && error[0] == 0);
  assert(near(a[cloud_ice_properties::ELIGIBLE_NUMBER],
              .25e-12 / particle_mass));
  const auto prior = a;
  exaero_bulk_cloud_ice(&p, 1, 1, 1, temperature.data(), humidity.data(),
                        ambient_mass.data(), cloud_mass.data(), -1, 0,
                        a.data(), c.data(), as.data(), cs.data(),
                        error.data(), &error_flag);
  assert(error_flag == 1 && a == prior);
  fractions.fill(-9);
  exaero_bulk_activation_fractions(
      &p, 1, 1, 1, 2, temperature.data(), humidity.data(),
      ambient_mass.data(), supersaturation.data(), fractions.data(),
      error.data(), &error_flag);
  assert(error_flag == 0);
  assert((fractions == std::array<double, 4>{0, 0, 1, 1}));
  exaero_bulk_activation_fractions(
      &p, 1, 1, 1, 0, temperature.data(), humidity.data(),
      ambient_mass.data(), supersaturation.data(), fractions.data(),
      error.data(), &error_flag);
  assert(error_flag == 1);
  assert((fractions == std::array<double, 4>{0, 0, 1, 1}));

  const auto old = a;
  ambient_mass[0] = -1;
  bool rejected = false;
  try {
    bulkCloudIce(p, env, ambient, nullptr, IceCategory::Dust, false,
                 ambient_out, cloud_out, ambient_status, cloud_status);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected && a == old);
}
