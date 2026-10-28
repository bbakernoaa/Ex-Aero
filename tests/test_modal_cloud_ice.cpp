#ifdef NDEBUG
#undef NDEBUG
#endif
#include <array>
#include <cassert>
#include <cmath>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/CloudIceCoupling_C.h>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemovalMath.hpp>
#include <stdexcept>

using namespace exaero;
namespace {
struct Package final : IAerosolPackage, IModalCoupling {
  ModalConfig config;
  Package() {
    config.species.resize(2);
    config.species[0].name = "DU";
    config.species[0].density = 1000;
    config.species[0].ice_categories = {IceCategory::Dust};
    config.species[1].name = "SU";
    config.species[1].density = 2000;
    config.species[1].ice_categories = {IceCategory::Sulfate};
    config.modes.resize(1);
    config.modes[0].name = "mixed";
    config.modes[0].geometric_radius = 1e-7;
    config.modes[0].geometric_std_dev = 1.5;
    config.modes[0].dry_geometric_radius_min = 5e-8;
    config.modes[0].dry_geometric_radius_max = 2e-7;
    config.modes[0].ice_rules = {{IceCategory::Dust,
                                 {IceSizeRuleKind::FixedFraction, .4},
                                 IceSizeRule{IceSizeRuleKind::All, 0}, true}};
  }
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
  int getSpeciesIndex(const std::string &) const override { return -1; }
  std::string getSpeciesName(int) const override { return {}; }
  AerosolRepresentation representation() const override {
    return AerosolRepresentation::Modal;
  }
  std::vector<StateField> stateFields() const override {
    return {{"DU_mixed_i", StateQuantity::Mass, 0, false},
            {"SU_mixed_i", StateQuantity::Mass, 0, false},
            {"num_mixed_i", StateQuantity::Number, 0, false}};
  }
  void initializeModal(const ModalConfig &v) override { config = v; }
  const ModalConfig &modalConfig() const override { return config; }
  std::vector<ModalOpticalTableInfo> opticalTables() const override {
    return {};
  }
  void modalProperties(const EnvironmentalStateView &,
                       const View3D<const double> &state, View4D<double> &out,
                       View4D<double> *) const override {
    const double volume = state(0, 0, 0) / 1000 + state(0, 0, 1) / 2000;
    const double number = state(0, 0, 2);
    const double radius =
        number > 0 ? modalRadius(volume, number, 1.5) : 1e-7;
    out(0, 0, 0, modal_properties::DRY_RADIUS) = radius;
    out(0, 0, 0, modal_properties::WET_RADIUS) =
        number > 0 ? radius * 1.2 : 0;
    out(0, 0, 0, modal_properties::KAPPA) = .5;
  }
  void modalUptake(const EnvironmentalStateView &, const View3D<const double> &,
                   std::span<const ModalVapor>, View4D<double> &,
                   View4D<double> *) const override {}
  ModalNucleationTerms nucleationTerms(double) const override { return {}; }
  double emittedNumber(std::size_t, std::span<const double>) const override {
    return 0;
  }
  void setModalOpticsProvider(
      std::shared_ptr<const IModalOpticsProvider>) override {}
};
bool near(double a, double b) {
  return std::abs(a - b) <= 1e-12 * std::max({1e-30, std::abs(a), std::abs(b)});
}
} // namespace

int main() {
  Package p;
  std::array<double, 1> temperature{280}, pressure{90000}, density{1},
      humidity{.5}, thickness{100};
  EnvironmentalStateView env{View2D<const double>(temperature.data(), 1, 1),
                             View2D<const double>(pressure.data(), 1, 1),
                             View2D<const double>(density.data(), 1, 1),
                             View2D<const double>(humidity.data(), 1, 1),
                             View2D<const double>(thickness.data(), 1, 1)};
  std::array<double, 3> ambient_data{1e-12, 1e-12, 1e7};
  std::array<double, 3> cloud_data{2e-12, 0, 5e6};
  View3D<const double> ambient(ambient_data.data(), 1, 1, 3);
  View3D<const double> cloud(cloud_data.data(), 1, 1, 3);
  std::array<double, cloud_ice_properties::COUNT> a{}, c{};
  std::array<int, 1> as{}, cs{};
  View4D<double> ambient_out(a.data(), 1, 1, 1,
                             cloud_ice_properties::COUNT);
  View4D<double> cloud_out(c.data(), 1, 1, 1, cloud_ice_properties::COUNT);
  View3D<int> ambient_status(as.data(), 1, 1, 1);
  View3D<int> cloud_status(cs.data(), 1, 1, 1);
  modalCloudIce(p, env, ambient, &cloud, IceCategory::Dust, false, ambient_out,
                cloud_out, ambient_status, cloud_status);
  const double radius = modalRadius(1.5e-15, 1e7, 1.5);
  assert(near(a[cloud_ice_properties::MASS], 2e-12));
  assert(near(a[cloud_ice_properties::NUMBER], 1e7));
  assert(near(a[cloud_ice_properties::DRY_NUMBER_MEDIAN_RADIUS], radius));
  assert(near(a[cloud_ice_properties::WET_SCAVENGING_RADIUS], radius * 1.2));
  assert(near(a[cloud_ice_properties::TYPE_WEIGHT], .5));
  assert(near(a[cloud_ice_properties::SIZE_WEIGHT], .4));
  assert(near(a[cloud_ice_properties::ELIGIBLE_NUMBER], 2e6));
  assert(near(c[cloud_ice_properties::TYPE_WEIGHT], 1));
  assert(near(c[cloud_ice_properties::ELIGIBLE_NUMBER], 2e6));
  assert(as[0] == static_cast<int>(CloudIceStatus::Configured));
  assert(cs[0] == static_cast<int>(CloudIceStatus::Configured));

  modalCloudIce(p, env, ambient, &cloud, IceCategory::Dust, true, ambient_out,
                cloud_out, ambient_status, cloud_status);
  assert(near(a[cloud_ice_properties::SIZE_WEIGHT], 1));
  assert(near(a[cloud_ice_properties::ELIGIBLE_NUMBER], 5e6));
  modalCloudIce(p, env, ambient, nullptr, IceCategory::Dust, false,
                ambient_out, cloud_out, ambient_status, cloud_status);
  assert(cs[0] == static_cast<int>(CloudIceStatus::Absent));
  for (double v : c)
    assert(v == 0);

  const double kelvin = 4 * .076 * .01801528 / (8.314462618 * 280 * 1000);
  const double critical =
      std::sqrt(4 * std::pow(kelvin, 3) / (27 * .5 * std::pow(2 * radius, 3)));
  std::array<double, 2> supersaturation{critical, 4 * critical};
  std::array<double, 4> fractions{};
  View1D<const double> ss(supersaturation.data(), 2);
  View4D<double> activation(fractions.data(), 1, 1, 1, 4);
  modalActivationFractions(p, env, ambient, ss, activation);
  assert(near(fractions[0], .5));
  assert(fractions[1] > .5 && fractions[1] < 1);
  assert(fractions[2] > fractions[0] && fractions[3] > fractions[1]);

  std::array<char, 256> error{};
  int error_flag = -1;
  exaero_modal_cloud_ice(&p, 1, 1, 3, temperature.data(), humidity.data(),
                         ambient_data.data(), cloud_data.data(),
                         EXAERO_ICE_DUST, 0, a.data(), c.data(), as.data(),
                         cs.data(), error.data(), &error_flag);
  assert(error_flag == 0 && error[0] == 0);
  assert(near(a[cloud_ice_properties::ELIGIBLE_NUMBER], 2e6));
  fractions.fill(-9);
  exaero_modal_activation_fractions(
      &p, 1, 1, 3, 2, temperature.data(), humidity.data(),
      ambient_data.data(), supersaturation.data(), fractions.data(),
      error.data(), &error_flag);
  assert(error_flag == 0 && near(fractions[0], .5));
  ambient_data[2] = 0;
  bool invalid_pair = false;
  try {
    modalActivationFractions(p, env, ambient, ss, activation);
  } catch (const std::invalid_argument &) {
    invalid_pair = true;
  }
  assert(invalid_pair);
  ambient_data[2] = 1e7;

  const auto old = a;
  ambient_data[2] = -1;
  bool rejected = false;
  try {
    modalCloudIce(p, env, ambient, nullptr, IceCategory::Dust, false,
                  ambient_out, cloud_out, ambient_status, cloud_status);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected && a == old);
}
