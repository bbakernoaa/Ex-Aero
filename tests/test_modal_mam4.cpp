/// @file test_modal_mam4.cpp
/// @brief Generic presets, configurable native water, and MAM4xx parity.
/// @author EX-aero contributors
/// @details History: 2026-10. Host layout_left fixtures; no file/device I/O.
#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <cmath>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalCoupling_C.h>
#include <exaero/ModalPresets.hpp>
#include <exaero/PackageFactory.hpp>
#include <numeric>

namespace {
using namespace exaero;
/// @brief Arbitrary-mode/material host fixture; Kokkos lifetime is in main TU.
struct MamFixture {
  static constexpr int nc = 2, nl = 3, np = nc * nl;
  ModalConfig config;
  std::unique_ptr<IAerosolPackage> package = createPackage("modal");
  std::array<double, np> t{}, p{}, rho{}, rh{}, dz{};
  std::vector<double> state;
  explicit MamFixture(ModalConfig c)
      : config(std::move(c)),
        state(np * config.modes.size() * (config.species.size() + 1)) {
    initializeModal(*package, config);
    for (int point = 0; point < np; ++point) {
      t[point] = 275 + point;
      p[point] = 90000 - 1000 * point;
      rho[point] = 1 + .1 * point;
      rh[point] = .5;
      dz[point] = 100 + 10 * point;
      for (std::size_t m = 0; m < config.modes.size(); ++m) {
        std::vector<double> masses(config.species.size());
        for (std::size_t s = 0; s < masses.size(); ++s) {
          masses[s] = 1e-10 * (point + 1) * (m + 1) * (s + 1);
          state[point + np * (m * (masses.size() + 1) + s)] = masses[s];
        }
        state[point + np * (m * (masses.size() + 1) + masses.size())] =
            modalCoupling(*package).emittedNumber(m, masses);
      }
    }
  }
  /// @brief Environmental host views, layout_left.
  EnvironmentalStateView env() const {
    return {View2D<const double>(t.data(), nc, nl),
            View2D<const double>(p.data(), nc, nl),
            View2D<const double>(rho.data(), nc, nl),
            View2D<const double>(rh.data(), nc, nl),
            View2D<const double>(dz.data(), nc, nl)};
  }
  /// @brief Immutable generic state view.
  View3D<const double> input() const {
    return View3D<const double>(state.data(), nc, nl, state.size() / np);
  }
  /// @brief Standard diagnostics through either package.
  std::vector<double> diagnostics(IAerosolPackage &pkg,
                                  const std::vector<double> &data) const {
    std::vector<double> out(np * diagnostic_indices::NUM_DIAGNOSTICS);
    View3D<double> output(out.data(), nc, nl,
                          diagnostic_indices::NUM_DIAGNOSTICS);
    pkg.computeDerivedDiagnostics(
        env(), View3D<const double>(data.data(), nc, nl, data.size() / np),
        output);
    return out;
  }
  /// @brief Advance generic state atomically.
  void advance(double dt) {
    View3D<double> view(state.data(), nc, nl, state.size() / np);
    package->executeMicrophysics(env(), view, dt);
  }
};
} // namespace

TEST_CASE("MAM4 preset is editable and available without the native backend") {
  auto config = mam4ModalConfig();
  REQUIRE(config.modes.size() == 4);
  REQUIRE(config.species.size() == 7);
  REQUIRE(config.species[1].density == 1770);
  REQUIRE(config.species[5].hygroscopicity == .14);
  auto yaml = createPackage("modal");
  yaml->initialize("backend: modal\npreset: mam4");
  REQUIRE(yaml->stateFields().size() == 32);
  REQUIRE(modalCoupling(*yaml).modalConfig().coagulation_routes.size() == 6);
  yaml->initialize("preset: mam4\nmode_overrides: {accumulation: "
                   "{geometric_std_dev: 1.7, geometric_radius: 6e-8}}");
  REQUIRE(modalCoupling(*yaml).modalConfig().modes[0].geometric_std_dev == 1.7);
  REQUIRE(modalCoupling(*yaml).modalConfig().modes[0].geometric_radius == 6e-8);
  MamFixture f(config);
  const auto before = f.state;
  f.advance(30);
  double before_mass = 0, after_mass = 0;
  for (std::size_t i = 0; i < before.size(); ++i)
    if ((i / MamFixture::np) % 8 != 7) {
      before_mass += before[i];
      after_mass += f.state[i];
    }
  REQUIRE(after_mass == Approx(before_mass).epsilon(1e-12));
  REQUIRE(f.state != before);
  const double wavelength = 550e-9;
  std::vector<double> optical(MamFixture::np * optical_indices::NUM_OPTICS);
  View4D<double> output(optical.data(), MamFixture::nc, MamFixture::nl, 1,
                        optical_indices::NUM_OPTICS);
  f.package->computeOptics(f.env(), f.input(),
                           View1D<const double>(&wavelength, 1), output);
  REQUIRE(optical[0] > 0);

  // Same factory remains arbitrary-mode: add a fifth mode and disable coag.
  config.coagulation = false;
  config.coagulation_routes.clear();
  config.modes.push_back({"extra", 4e-8, 1.4, {}});
  config.modes[0].geometric_std_dev = 1.7;
  MamFixture extended(config);
  REQUIRE(extended.package->stateFields().size() == 40);
  REQUIRE(extended.package->stateFields()[7].name == "num_accumulation_i");
  const auto diag = extended.diagnostics(*extended.package, extended.state);
  REQUIRE(diag[diagnostic_indices::MASS_CONCENTRATION] > 0);
  std::array<char, 256> error{};
  int status = -1;
  auto *handle = exaero_create_package("modal", error.data(), &status);
  REQUIRE(handle);
  exaero_initialize_package(handle,
                            "preset: mam4\nprocesses: {coagulation: none, "
                            "water_uptake: kappa_equilibrium}",
                            error.data(), &status);
  const bool initialized = status == 0;
  const int modes = exaero_get_modal_mode_count(handle);
  exaero_free_package(handle);
  REQUIRE(initialized);
  REQUIRE(modes == 4);
}

TEST_CASE(
    "process configuration rejects unsupported selections transactionally") {
  auto package = createPackage("modal");
  package->initialize("preset: mam4");
  for (const auto *yaml :
       {"preset: invalid", "preset: mam4\nconstants: {gravity: -1}",
        "preset: mam4\nmode_overrides: {missing: {geometric_std_dev: 1.7}}",
        "preset: mam4\nmode_overrides: {aitken: {typo: 1.7}}",
        "preset: mam4\nmode_overrides: {aitken: {geometric_std_dev: 1}}",
        "preset: mam4\nprocesses: {aging: mam4xx}",
        "preset: mam4\nprocesses: {coagulation: invalid}",
        "preset: mam4\nprocesses: {water_uptake: invalid}",
        "preset: mam4\ncoagulation: false\nprocesses: {coagulation: none}"}) {
    REQUIRE_THROWS(package->initialize(yaml));
    REQUIRE(package->stateFields().size() == 32);
  }
  auto bad = mam4ModalConfig();
  bad.modes[0].crystallization_rh = .9;
  bad.modes[0].deliquescence_rh = .8;
  REQUIRE_THROWS(initializeModal(*package, bad));
  if (!packageAvailable("mam4xx")) {
    REQUIRE_THROWS_WITH(package->initialize("preset: mam4xx_reference"),
                        Catch::Matchers::Contains("EXAERO_ENABLE_MAM4XX"));
    REQUIRE(package->stateFields().size() == 32);
    return;
  }
  package->initialize("preset: mam4xx_reference");
  bad = mam4xxReferenceModalConfig();
  bad.modes[0].geometric_std_dev = 1.7;
  REQUIRE_THROWS(initializeModal(*package, bad));
  bad = mam4xxReferenceModalConfig();
  bad.species[0].density *= 2;
  REQUIRE_THROWS(initializeModal(*package, bad));
  bad = mam4xxReferenceModalConfig();
  bad.water_uptake = ModalWaterUptakeModel::KappaEquilibrium;
  REQUIRE_THROWS(initializeModal(*package, bad));
  bad = mam4xxReferenceModalConfig();
  bad.coagulation_routes.push_back({"aitken", "coarse", "coarse"});
  REQUIRE_THROWS(initializeModal(*package, bad));
}

TEST_CASE(
    "generic MAM4 reference matches native water diagnostics and coagulation") {
  if (!packageAvailable("mam4xx"))
    return;
  MamFixture f(mam4xxReferenceModalConfig());
  auto native = createPackage("mam4xx");
  native->initialize("backend: mam4xx");
  std::vector<double> state(f.state.size() * 2, 0);
  std::copy(f.state.begin(), f.state.end(), state.begin());
  for (double rh : {0., .2, .35, .5, .8, .95, .98}) {
    f.rh.fill(rh);
    const auto actual = f.diagnostics(*f.package, f.state);
    const auto expected = f.diagnostics(*native, state);
    for (std::size_t i = 0; i < actual.size(); ++i) {
      INFO(rh);
      INFO(i);
      REQUIRE(actual[i] == Approx(expected[i]).epsilon(2e-7).margin(1e-24));
    }
  }
  const std::array<double, 3> ss{.001, .005, .01};
  std::vector<double> actual_ccn(MamFixture::np * ss.size()),
      expected_ccn(actual_ccn.size());
  View4D<double> ac(actual_ccn.data(), MamFixture::nc, MamFixture::nl,
                    ss.size(), 1);
  View4D<double> ec(expected_ccn.data(), MamFixture::nc, MamFixture::nl,
                    ss.size(), 1);
  const View1D<const double> axis(ss.data(), ss.size());
  f.package->computeCCN(f.env(), f.input(), axis, ac);
  native->computeCCN(
      f.env(),
      View3D<const double>(state.data(), MamFixture::nc, MamFixture::nl, 64),
      axis, ec);
  for (std::size_t i = 0; i < actual_ccn.size(); ++i)
    REQUIRE(actual_ccn[i] == Approx(expected_ccn[i]).epsilon(1e-10));
  f.rh.fill(.6);
  View3D<double> view(state.data(), MamFixture::nc, MamFixture::nl, 64);
  native->executeMicrophysics(f.env(), view, 60);
  f.advance(60);
  for (std::size_t i = 0; i < f.state.size(); ++i)
    REQUIRE(f.state[i] == Approx(state[i]).epsilon(1e-13));
  const auto before = f.state;
  f.advance(0);
  REQUIRE(f.state == before);
  REQUIRE_THROWS(f.advance(-1));
  REQUIRE(f.state == before);
}

TEST_CASE("native water supports arbitrary modes, materials and dry uptake "
          "derivatives") {
  if (!packageAvailable("mam4xx"))
    return;
  auto c = mam4ModalConfig();
  c.coagulation = false;
  c.coagulation_routes.clear();
  c.water_uptake = ModalWaterUptakeModel::Mam4xx;
  c.species.resize(2);
  for (auto &s : c.species)
    s.optics.reset();
  c.modes.resize(1);
  c.modes[0].crystallization_rh = .1;
  c.modes[0].deliquescence_rh = .4;
  MamFixture f(c);
  auto &coupling = modalCoupling(*f.package);
  std::vector<double> properties(MamFixture::np * modal_properties::COUNT, -1);
  View4D<double> out(properties.data(), MamFixture::nc, MamFixture::nl, 1,
                     modal_properties::COUNT);
  f.rh.fill(.05);
  coupling.modalProperties(f.env(), f.input(), out);
  REQUIRE(properties[modal_properties::WATER * MamFixture::np] == 0);
  f.rh.fill(.5);
  coupling.modalProperties(f.env(), f.input(), out);
  REQUIRE(properties[modal_properties::WATER * MamFixture::np] > 0);
  std::vector<double> jac(properties.size() * 3, -1);
  View4D<double> dj(jac.data(), MamFixture::nc, MamFixture::nl, 1,
                    modal_properties::COUNT * 3);
  const auto saved = properties;
  REQUIRE_THROWS(coupling.modalProperties(f.env(), f.input(), out, &dj));
  REQUIRE(properties == saved);
  REQUIRE(
      std::all_of(jac.begin(), jac.end(), [](double v) { return v == -1; }));
  const std::array<ModalVapor, 1> vapor{{{1e-5, .098079, 1}}};
  std::array<double, MamFixture::np> uptake{};
  std::array<double, MamFixture::np * 3> derivative{};
  View4D<double> up(uptake.data(), MamFixture::nc, MamFixture::nl, 1, 1);
  View4D<double> du(derivative.data(), MamFixture::nc, MamFixture::nl, 1, 3);
  coupling.modalUptake(f.env(), f.input(), vapor, up, &du);
  REQUIRE(uptake[0] > 0);
  REQUIRE(std::isfinite(derivative[0]));
  const double original = f.state[0], step = original * 1e-5;
  f.state[0] = original + step;
  coupling.modalUptake(f.env(), f.input(), vapor, up);
  const double plus = uptake[0];
  f.state[0] = original - step;
  coupling.modalUptake(f.env(), f.input(), vapor, up);
  const double minus = uptake[0];
  f.state[0] = original;
  REQUIRE(derivative[0] == Approx((plus - minus) / (2 * step)).epsilon(1e-6));
  for (int m = 1; m < 5; ++m) {
    auto mode = c.modes[0];
    mode.name = "extra" + std::to_string(m);
    mode.geometric_std_dev = 1.2 + .1 * m;
    c.modes.push_back(mode);
  }
  MamFixture extended(c);
  const auto diagnostics =
      extended.diagnostics(*extended.package, extended.state);
  REQUIRE(extended.package->stateFields().size() == 15);
  REQUIRE(
      diagnostics[diagnostic_indices::AEROSOL_LIQUID_WATER * MamFixture::np] >
      0);
}
