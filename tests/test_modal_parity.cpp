/// @file test_modal_parity.cpp
/// @brief Source-derived CAM-SIMA modal formulas and host-coupling invariants.
#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <cmath>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemoval.hpp>
#include <exaero/PackageFactory.hpp>
#include <numeric>

namespace {
using namespace exaero;
constexpr double parity_pi = 3.14159265358979323846;
/// @brief Minimal physical mixture, with no chemistry or implicit transfer.
ModalConfig parityConfig(std::size_t nm = 1) {
  ModalConfig cfg;
  cfg.coagulation = false;
  cfg.species = {{"A", 1000, .1, .5, {}}, {"B", 2000, .2, .5, {}}};
  for (std::size_t m = 0; m < nm; ++m) {
    ModalModeConfig mode;
    mode.name = "mode" + std::to_string(m);
    mode.geometric_radius = 1e-7 * (m + 1);
    mode.geometric_std_dev = 1.5;
    mode.dry_geometric_radius_min = 5e-8;
    mode.dry_geometric_radius_max = 2e-7;
    cfg.modes.push_back(mode);
  }
  return cfg;
}
/// @brief Single grid point with a configurable count of mixed modes.
struct ParityState {
  ModalConfig cfg;
  std::unique_ptr<IAerosolPackage> package = createPackage("modal");
  std::vector<double> values;
  std::array<double, 1> t{280}, p{90000}, rho{1.1}, rh{.5}, dz{10};
  explicit ParityState(ModalConfig c = parityConfig())
      : cfg(std::move(c)), values(cfg.modes.size() * 3) {
    initializeModal(*package, cfg);
  }
  EnvironmentalStateView env() const {
    return {View2D<const double>(t.data(), 1, 1),
            View2D<const double>(p.data(), 1, 1),
            View2D<const double>(rho.data(), 1, 1),
            View2D<const double>(rh.data(), 1, 1),
            View2D<const double>(dz.data(), 1, 1)};
  }
  View3D<const double> input() const {
    return View3D<const double>(values.data(), 1, 1, values.size());
  }
  std::vector<double> properties() {
    std::vector<double> result(cfg.modes.size() * modal_properties::COUNT);
    View4D<double> out(result.data(), 1, 1, cfg.modes.size(),
                       modal_properties::COUNT);
    modalCoupling(*package).modalProperties(env(), input(), out);
    return result;
  }
  void populate(double radius = 1e-7, double number = 1e8) {
    const double l = std::log(1.5);
    const double volume = number * 4 * parity_pi / 3 * std::pow(radius, 3) *
                          std::exp(4.5 * l * l);
    values = {volume * .5 * 1000, volume * .5 * 2000, number};
  }
};
} // namespace

TEST_CASE("CAM-SIMA source-derived number limits preserve component mass",
          "[parity][bounds]") {
  ParityState f;
  for (double radius : {2.5e-8, 5e-8, 1e-7, 2e-7, 4e-7}) {
    CAPTURE(radius);
    f.populate(radius);
    const auto original = f.values;
    std::array<double, 3> report{};
    double delta = -99;
    int status = -99;
    View4D<double> rv(report.data(), 1, 1, 1, 3);
    View3D<int> sv(&status, 1, 1, 1);
    inspectModalBounds(*f.package, f.input(), rv, sv);
    REQUIRE(f.values == original);
    const auto expected_status = radius < 5e-8   ? ModalBoundStatus::BelowRadius
                                 : radius > 2e-7 ? ModalBoundStatus::AboveRadius
                                                 : ModalBoundStatus::Within;
    REQUIRE(status == static_cast<int>(expected_status));
    const double expected =
        1e8 * std::pow(radius / std::clamp(radius, 5e-8, 2e-7), 3);
    View3D<double> state(f.values.data(), 1, 1, 3), dv(&delta, 1, 1, 1);
    enforceModalBounds(*f.package, state, dv, sv);
    REQUIRE(f.values[0] == original[0]);
    REQUIRE(f.values[1] == original[1]);
    REQUIRE(f.values[2] == Approx(expected).epsilon(1e-12));
    REQUIRE(delta == Approx(expected - original[2]).margin(1e-4));
    enforceModalBounds(*f.package, state, dv, sv);
    REQUIRE(delta == 0);
  }
  f.values.assign(3, 0);
  double delta = -1;
  int status = -1;
  View3D<double> state(f.values.data(), 1, 1, 3), dv(&delta, 1, 1, 1);
  View3D<int> sv(&status, 1, 1, 1);
  enforceModalBounds(*f.package, state, dv, sv);
  REQUIRE(status == static_cast<int>(ModalBoundStatus::Empty));
  REQUIRE(delta == 0);
  REQUIRE(f.values == std::vector<double>(3, 0));
}

TEST_CASE("Current modal scavenging sizes follow kappa growth and RH cap",
          "[parity][sizes]") {
  ParityState f;
  f.populate();
  const auto original = f.values;
  for (double humidity : {0., .5, .98, 1.}) {
    f.rh[0] = humidity;
    std::array<double, 2> diameters{};
    int populated = -1;
    View4D<double> out(diameters.data(), 1, 1, 1, 2);
    View3D<int> pop(&populated, 1, 1, 1);
    modalScavengingSizes(*f.package, f.env(), f.input(), out, pop);
    const double capped = std::min(humidity, .98);
    REQUIRE(diameters[0] == Approx(2e-7).epsilon(1e-12));
    REQUIRE(diameters[1] ==
            Approx(2e-7 * std::cbrt(1 + .5 * capped / (1 - capped)))
                .epsilon(1e-12));
    REQUIRE(populated == 1);
    REQUIRE(f.values == original);
  }
}

TEST_CASE(
    "Composition changes the modal properties supplied to host scavenging",
    "[parity][composition]") {
  auto cfg = parityConfig();
  cfg.species[0].hygroscopicity = .6;
  cfg.species[1].hygroscopicity = 0;
  ParityState f(cfg);
  f.populate();
  const double volume = f.values[0] / 1000 + f.values[1] / 2000;
  double last_wet = 0;
  for (double fraction : {.75, .25}) {
    f.values[0] = volume * fraction * 1000;
    f.values[1] = volume * (1 - fraction) * 2000;
    const auto original = f.values, props = f.properties();
    REQUIRE(props[modal_properties::KAPPA] ==
            Approx(.6 * fraction).epsilon(1e-12));
    REQUIRE(props[modal_properties::DRY_RADIUS] == Approx(1e-7).epsilon(1e-12));
    REQUIRE(props[modal_properties::WET_RADIUS] ==
            Approx(1e-7 * std::cbrt(1 + .6 * fraction)).epsilon(1e-12));
    if (last_wet > 0)
      REQUIRE(props[modal_properties::WET_RADIUS] < last_wet);
    last_wet = props[modal_properties::WET_RADIUS];
    REQUIRE(f.values == original);
  }
}

TEST_CASE("CAM-SIMA activation reference has distinct number and mass tails",
          "[parity][activation]") {
  ParityState f;
  f.populate();
  const auto original = f.values;
  const auto &constants = f.cfg.constants;
  const double kelvin = 4 * constants.surface_tension *
                        constants.water_molecular_weight /
                        (constants.gas_constant * 280 * 1000);
  const double critical =
      std::sqrt(4 * std::pow(kelvin, 3) / (27 * .5 * std::pow(2e-7, 3)));
  const std::array<double, 3> ss{critical / 2, critical,
                                 std::min(2 * critical, 1.)};
  std::array<double, 6> fractions{};
  View4D<double> out(fractions.data(), 1, 1, 1, 6);
  modalActivationFractions(*f.package, f.env(), f.input(),
                           View1D<const double>(ss.data(), 3), out);
  REQUIRE(fractions[2] == Approx(.5).epsilon(1e-12));
  REQUIRE(fractions[3] ==
          Approx(.5 * (1 + std::erf(1.5 * std::sqrt(2.) * std::log(1.5))))
              .epsilon(1e-12));
  REQUIRE(fractions[0] < fractions[2]);
  REQUIRE(fractions[2] < fractions[4]);
  REQUIRE(fractions[1] < fractions[3]);
  REQUIRE(fractions[3] < fractions[5]);
  REQUIRE(f.values == original);
}

TEST_CASE("Explicit resuspension routes reproduce largest-mode budgets and "
          "read sources once",
          "[parity][resuspension]") {
  for (bool chained : {false, true}) {
    auto cfg = parityConfig(3);
    cfg.resuspension_routes = {{"A", "mode0", chained ? "mode1" : "mode2"},
                               {"A", "mode1", "mode2"}};
    ParityState f(cfg);
    const std::array<double, 9> rates{10, 4, 0, 3, 0, 0, 0, 0, 0};
    std::array<double, 9> routed{};
    View3D<double> out(routed.data(), 1, 1, 9);
    redistributeModalResuspension(
        *f.package, View3D<const double>(rates.data(), 1, 1, 9), out);
    REQUIRE(routed[0] == 0);
    REQUIRE(routed[1] == 4);
    REQUIRE(routed[3] == (chained ? 10 : 0));
    REQUIRE(routed[6] == (chained ? 3 : 13));
    REQUIRE(std::accumulate(routed.begin(), routed.end(), 0.) == 17);
    REQUIRE(routed[2] == 0);
    REQUIRE(routed[5] == 0);
    REQUIRE(routed[8] == 0);
  }
}

TEST_CASE("CAM-SIMA single-material deposited mass CDF depends on number and "
          "closes budget",
          "[parity][rebin]") {
  ParityState f;
  f.populate();
  f.values[0] += f.values[1] / 2;
  f.values[1] = 0;
  const auto original = f.values;
  const std::array<double, 4> edges{1e-9, 1e-7, 2e-7, 4e-7};
  std::array<double, 3> bins{};
  std::array<double, 2> tails{};
  View3D<double> out(bins.data(), 1, 1, 3), tail(tails.data(), 1, 1, 2);
  const auto rebin = [&] {
    rebinModalDeposition(*f.package, f.input(), 0,
                         View1D<const double>(edges.data(), edges.size()), out,
                         tail);
    REQUIRE(std::accumulate(bins.begin(), bins.end(), tails[0] + tails[1]) ==
            Approx(f.values[0]).epsilon(1e-12));
  };
  rebin();
  const auto original_bins = bins;
  const double ls = std::log(1.5);
  const double vmd =
      std::cbrt(6 * f.values[0] / (parity_pi * 1000 * f.values[2])) *
      std::exp(1.5 * ls * ls);
  for (std::size_t j = 0; j < bins.size(); ++j) {
    const double fraction =
        .5 * (std::erf(std::log(edges[j + 1] / vmd) / (std::sqrt(2.) * ls)) -
              std::erf(std::log(edges[j] / vmd) / (std::sqrt(2.) * ls)));
    REQUIRE(bins[j] == Approx(f.values[0] * fraction).epsilon(1e-12));
  }
  for (double &x : f.values)
    x *= 2;
  rebin();
  for (std::size_t j = 0; j < bins.size(); ++j)
    REQUIRE(bins[j] == Approx(2 * original_bins[j]).epsilon(1e-12));
  f.values = original;
  f.values[2] *= 8;
  rebin();
  REQUIRE(bins[0] > original_bins[0]);
}
