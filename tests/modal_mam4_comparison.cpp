/// @file modal_mam4_comparison.cpp
/// @brief Deterministic dry box comparisons of generic closures and MAM4xx.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: executable [--check]. Input: built-in
/// cases. Output: CSV on stdout, redirected by the report driver under DATA.
/// Host-only validation driver; public scientific views use layout_left.
#include "modal_mam4_reference.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <exaero/ModalChemistry.hpp>
#include <exaero/ModalPresets.hpp>
#include <exaero/PackageFactory.hpp>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
using namespace exaero;
constexpr int nm = 4, ns = 7, nf = 32;
constexpr double avogadro = 6.02214076e23, acid_mw = .098079,
                 sulfate_mw = .09606, seed_molecules = 100;
using Fields = std::array<double, nf>;
using Modes = std::array<double, nm>;
/// @brief Fail validation without relying on assertions disabled in Release.
void require(bool ok, std::string_view message) {
  if (!ok)
    throw std::runtime_error(std::string(message));
}
/// @brief Named process switches for controlled and scientific comparisons.
struct Experiment {
  const char *name;
  bool coag = false, condensation = false, nucleation = false, aging = false;
  bool native_coag_control = false;
  double exponent = 1, prefactor_cm = 1e-6;
};
const std::array<Experiment, 9> experiments{{
    {"coagulation_control", true, false, false, false, true},
    {"coagulation", true},
    {"condensation", false, true},
    {"nucleation_linear", false, false, true},
    {"nucleation_quadratic", false, false, true, false, false, 2, 1e-12},
    {"nucleation_example", false, false, true, false, false, 1.5, 1e-8},
    {"aging", false, false, false, true},
    {"combined", true, true, true, true, false, 1.5, 1e-8},
    {"combined_wang_linear", true, true, true, true},
}};
/// @brief Uniform dry environment; three cases vary transport and acid supply.
struct Case {
  double temperature, pressure, acid_cm, size_scale;
};
const std::array<Case, 3> cases{
    {{250, 60000, 1e5, .7}, {280, 90000, 1e6, 1}, {300, 101325, 1e7, 1.5}}};
/// @brief Host box with equal initial number, mass and material basis in both
/// paths. No cloudborne state, optical calculations, or wet growth is enabled.
struct Box {
  ModalConfig config;
  std::unique_ptr<IAerosolPackage> package = createPackage("modal");
  std::unique_ptr<IAerosolPackage> nucleator;
  Fields state{};
  double temperature, pressure, density, rh = 0, dz = 100, gas;
  bool native;
  /// @brief Construct a four-mode box; primary carbon starts with POM/BC only.
  Box(const Case &c, const Experiment &e, bool reference)
      : config((reference || e.native_coag_control)
                   ? mam4xxReferenceModalConfig()
                   : mam4ModalConfig()),
        temperature(c.temperature), pressure(c.pressure),
        density(pressure / (287.05 * temperature)),
        gas(c.acid_cm * 1e6 * sulfate_mw / avogadro), native(reference) {
    config.coagulation = e.coag;
    config.transfer_owner = ModalTransferOwner::Host;
    config.aging_routes = {
        {"primary_carbon", "accumulation", {"SO4"}, 8 * 4.76e-10}};
    initializeModal(*package, config);
    if (e.nucleation && !reference) {
      // EX-aero's empirical helper uses actual pure-H2SO4 seed mass. Keep that
      // basis explicit, then map its acid loss to the public SO4 equivalents.
      auto nucleation_config = config;
      nucleation_config.coagulation = false;
      nucleation_config.species[1].molecular_weight = acid_mw;
      const double width = std::log(config.modes[1].geometric_std_dev);
      const double volume =
          seed_molecules * acid_mw / avogadro / config.species[1].density;
      const double radius = std::cbrt(
          volume / (4 * std::numbers::pi / 3 * std::exp(4.5 * width * width)));
      std::vector<double> fractions(ns, 0);
      fractions[1] = 1;
      nucleation_config.nucleation = ModalNucleationConfig{
          "aitken",   e.prefactor_cm * std::pow(1e6, 1 - e.exponent),
          e.exponent, radius,
          fractions,  acid_mw};
      nucleator = createPackage("modal");
      initializeModal(*nucleator, nucleation_config);
    }
    constexpr Modes numbers{2e8, 8e8, 1e5, 1e8};
    for (int m = 0; m < nm; ++m) {
      const double width = std::log(config.modes[m].geometric_std_dev);
      const double radius = config.modes[m].geometric_radius * c.size_scale;
      const double volume = numbers[m] * 4 * std::numbers::pi / 3 *
                            std::pow(radius, 3) * std::exp(4.5 * width * width);
      if (m == 3) {
        state[m * 8 + 2] = .8 * volume * config.species[2].density;
        state[m * 8 + 3] = .2 * volume * config.species[3].density;
      } else {
        state[m * 8 + 1] = volume * config.species[1].density;
      }
      state[m * 8 + 7] = numbers[m];
    }
  }
  /// @brief One-cell layout_left environment for public EX-aero calls.
  EnvironmentalStateView env() const {
    return {View2D<const double>(&temperature, 1, 1),
            View2D<const double>(&pressure, 1, 1),
            View2D<const double>(&density, 1, 1),
            View2D<const double>(&rh, 1, 1), View2D<const double>(&dz, 1, 1)};
  }
  /// @brief Immutable one-cell aerosol state view.
  View3D<const double> input() const {
    return View3D<const double>(state.data(), 1, 1, nf);
  }
  /// @brief Modal dry properties from the same public diagnostic interface.
  std::array<double, nm * modal_properties::COUNT> properties() const {
    std::array<double, nm * modal_properties::COUNT> result{};
    View4D<double> out(result.data(), 1, 1, nm, modal_properties::COUNT);
    modalCoupling(*package).modalProperties(env(), input(), out);
    return result;
  }
  /// @brief Per-species material budget including gas on the SO4-equivalent
  /// basis. Species mass/number conversion is deliberately outside chemistry.
  std::array<double, ns> budget() const {
    std::array<double, ns> result{};
    result[1] = gas;
    for (int m = 0; m < nm; ++m)
      for (int s = 0; s < ns; ++s)
        result[s] += state[m * 8 + s];
    return result;
  }
};
/// @brief Dry uptake: native two-point beta-adjusted quadrature versus the
/// generic EX-aero transition-regime integral. Use identical diffusivity/alpha.
Modes uptake(Box &b) {
  const double diffusivity =
      modal_comparison::diffusivity(b.temperature, b.pressure);
  Modes result{};
  if (b.native) {
    const auto data = b.properties();
    const View4D<const double> p(data.data(), 1, 1, nm,
                                 modal_properties::COUNT);
    Modes radius{}, width{};
    for (int m = 0; m < nm; ++m) {
      radius[m] = p(0, 0, m, modal_properties::DRY_RADIUS);
      width[m] = b.config.modes[m].geometric_std_dev;
    }
    result =
        modal_comparison::uptake(b.temperature, diffusivity, radius, width);
    for (int m = 0; m < nm; ++m)
      result[m] *= b.state[m * 8 + 7];
  } else {
    const std::array<ModalVapor, 1> vapor{{{diffusivity, acid_mw, .65}}};
    View4D<double> out(result.data(), 1, 1, nm, 1);
    modalCoupling(*b.package).modalUptake(b.env(), b.input(), vapor, out);
  }
  return result;
}
/// @brief Cluster formation J from native Wang first/second-order option.
/// Survival/growth to Aitken is excluded; both paths use a controlled 100-unit
/// seed so rate and material bookkeeping can be compared independently.
double nativeJ(double acid_cm, int choice) {
  return modal_comparison::nucleation(acid_cm, choice);
}
/// @brief Actual configurable EX-aero nucleation helper with seed-yield check.
double genericJ(Box &b) {
  const auto terms =
      modalCoupling(*b.nucleator).nucleationTerms(b.gas * acid_mw / sulfate_mw);
  require(std::abs(terms.precursor_loss -
                   terms.number_rate * seed_molecules * acid_mw / avogadro) <=
              1e-12 * std::max(terms.precursor_loss, 1e-30),
          "generic nucleation seed material conversion failed");
  return terms.number_rate;
}
/// @brief Map public dry radii to the isolated native finite-step aging kernel.
void nativeAging(Box &b, const Fields &condensation, const Fields &coag_gain,
                 double dt) {
  const auto data = b.properties();
  const View4D<const double> p(data.data(), 1, 1, nm, modal_properties::COUNT);
  Modes radius{};
  for (int m = 0; m < nm; ++m)
    radius[m] = p(0, 0, m, modal_properties::DRY_RADIUS);
  modal_comparison::aging(b.temperature, b.pressure, b.state, radius,
                          condensation, coag_gain, dt);
}
/// @brief Generic shell loss frequency from the public chemistry interface.
double agingFrequency(Box &b, const Fields &condensation,
                      const Fields &coag_gain, double dt) {
  Fields rates{}, coag_rate{};
  for (int i = 0; i < nf; ++i)
    coag_rate[i] = coag_gain[i] / dt;
  View3D<double> out(rates.data(), 1, 1, nf);
  modalAgingTerms(*b.package, b.env(), b.input(),
                  View3D<const double>(condensation.data(), 1, 1, nf),
                  View3D<const double>(coag_rate.data(), 1, 1, nf), out);
  const double n = b.state[31];
  return n > 0 ? -rates[31] / n : 0;
}
/// @brief Evaluate the frozen stage RHS; optional MKPP build executes generated
/// YAML kernels and checks every physical row against the independent host RHS.
/// Returns mass/number rates in the common physical state and gas rate.
std::array<double, nf + 1> stageRhs(const Box &b, const Modes &k, double j,
                                    double age) {
  std::array<double, nf + 1> expected{};
  for (int m = 0; m < nm; ++m) {
    expected[m * 8 + 1] += k[m] * b.gas;
    expected[nf] -= k[m] * b.gas;
  }
  const double nucleated = j * seed_molecules * sulfate_mw / avogadro;
  expected[8 + 1] += nucleated;
  expected[8 + 7] += j;
  expected[nf] -= nucleated;
  for (int s = 0; s < 8; ++s) {
    const double transferred = age * b.state[24 + s];
    expected[24 + s] -= transferred;
    expected[s] += transferred;
  }
#ifdef EXAERO_COMPARISON_MKPP
  std::array<double, 39> state{};
  state[0] = b.gas * avogadro / sulfate_mw / 1e6;
  for (int m = 0; m < nm; ++m)
    for (int s = 0; s < 8; ++s) {
      const double factor =
          s == 7 ? 1 : avogadro / b.config.species[s].molecular_weight;
      state[1 + m * 8 + s] = b.state[m * 8 + s] * factor / 1e6;
    }
  for (int m = 0; m < nm; ++m)
    state[33 + m] = k[m];
  state[37] = j / 1e6;
  state[38] = age;
  const auto rates = modal_comparison::generatedRhs(state, b.temperature);
  std::array<double, nf + 1> result{};
  result[nf] = rates[0] * 1e6 * sulfate_mw / avogadro;
  for (int m = 0; m < nm; ++m)
    for (int s = 0; s < 8; ++s) {
      const double factor =
          s == 7 ? 1 : b.config.species[s].molecular_weight / avogadro;
      result[m * 8 + s] = rates[1 + m * 8 + s] * factor * 1e6;
    }
  for (int i = 0; i < nf + 1; ++i)
    require(std::abs(result[i] - expected[i]) <=
                1e-11 * std::max(std::abs(expected[i]), 1e-30),
            "generated MKPP RHS differs from independent host rate");
  for (int i = 33; i < 39; ++i)
    require(rates[i] == 0, "generated MKPP modified a fixed coefficient");
  return result;
#else
  return expected;
#endif
}
/// @brief One first-order split step. Identical Euler chemistry integration
/// isolates closure differences. Coagulation uses each package's native finite
/// update; both cases pass its accepted gross gains to their aging closure.
void advance(Box &b, const Experiment &e, double dt) {
  Fields gains{}, condensation{};
  if (e.coag) {
    if (b.native) {
      const auto properties = b.properties();
      const View4D<const double> p(properties.data(), 1, 1, nm,
                                   modal_properties::COUNT);
      Modes radius{}, density{};
      for (int m = 0; m < nm; ++m) {
        radius[m] = p(0, 0, m, modal_properties::DRY_RADIUS);
        density[m] = p(0, 0, m, modal_properties::WET_DENSITY);
      }
      gains = modal_comparison::coagulation(b.temperature, b.pressure, b.state,
                                            radius, density, dt);
    } else {
      View3D<double> state(b.state.data(), 1, 1, nf),
          gain(gains.data(), 1, 1, nf);
      advanceModalWithCoating(*b.package, b.env(), state, dt, gain);
    }
  }
  const auto k = e.condensation ? uptake(b) : Modes{};
  for (int m = 0; m < nm; ++m)
    condensation[m * 8 + 1] = k[m] * b.gas;
  // Aging-only case uses a prescribed, budgeted sulfate coating source.
  if (e.aging && !e.condensation)
    condensation[25] = 1e-15;
  const double acid_cm = b.gas * avogadro / sulfate_mw / 1e6;
  const double j = !e.nucleation ? 0
                   : b.native    ? nativeJ(acid_cm, e.exponent == 2 ? 12 : 11)
                                 : genericJ(b);
  // Condensation/nucleation first, then coating aging. Both paths use exactly
  // the same ordering and explicit chemistry substep.
  const auto chemical = stageRhs(b, k, j, 0);
  require(b.gas + dt * chemical[nf] >= 0,
          "gas exhausted: reduce comparison dt");
  for (int i = 0; i < nf; ++i)
    b.state[i] += dt * chemical[i];
  b.gas += dt * chemical[nf];
  if (e.aging && !e.condensation)
    b.state[25] += dt * condensation[25];
  if (e.aging) {
    if (b.native) {
      nativeAging(b, condensation, gains, dt);
    } else {
      const double frequency = agingFrequency(b, condensation, gains, dt);
      require(dt * frequency <= 1,
              "aging exhausted source: reduce comparison dt");
      const auto aging = stageRhs(b, Modes{}, 0, frequency);
      for (int i = 0; i < nf; ++i)
        b.state[i] += dt * aging[i];
    }
  }
  require(std::isfinite(b.gas) && b.gas >= 0, "invalid gas state");
  for (double value : b.state)
    require(std::isfinite(value) && value >= 0, "invalid aerosol state");
}
/// @brief Emit paired values and units; no rounded pass/fail science
/// thresholds.
void row(const Experiment &e, int c, double time, std::string_view field,
         std::string_view unit, double generic, double native) {
  std::cout << e.name << ',' << c << ',' << time << ',' << field << ',' << unit
            << ',' << generic << ',' << native << ','
#ifdef EXAERO_COMPARISON_MKPP
            << "mkpp_generated"
#else
            << "host_reference"
#endif
            << '\n';
}
/// @brief Export all transported fields plus gas and modal dry radii.
void snapshot(const Experiment &e, int c, double time, const Box &g,
              const Box &n) {
  const auto fields = g.package->stateFields();
  for (int i = 0; i < nf; ++i)
    row(e, c, time, fields[i].name, i % 8 == 7 ? "particles/m3" : "kg/m3",
        g.state[i], n.state[i]);
  row(e, c, time, "H2SO4_SO4_equivalent", "kg/m3", g.gas, n.gas);
  const auto gp = g.properties(), np = n.properties();
  const View4D<const double> gv(gp.data(), 1, 1, nm, modal_properties::COUNT),
      nv(np.data(), 1, 1, nm, modal_properties::COUNT);
  for (int m = 0; m < nm; ++m)
    row(e, c, time, "dry_radius_" + g.config.modes[m].name, "m",
        gv(0, 0, m, modal_properties::DRY_RADIUS),
        nv(0, 0, m, modal_properties::DRY_RADIUS));
}
/// @brief Validate every material's budget; prescribed aging coating is the
/// only external source. Also enforce matched-law parity in controlled cases.
void validate(const Box &b, const std::array<double, ns> &initial,
              const Experiment &e, double time) {
  const auto final = b.budget();
  for (int s = 0; s < ns; ++s) {
    const double expected =
        initial[s] + (e.aging && !e.condensation && s == 1 ? time * 1e-15 : 0);
    require(std::abs(final[s] - expected) <= 1e-10 * std::max(expected, 1e-25),
            "material budget failed");
  }
}
} // namespace

/// @brief Kokkos lifetime and bounded regression/report driver. CPU host only.
int runComparison(int argc, char **argv) {
  try {
    const bool check = argc == 2 && std::string_view(argv[1]) == "--check";
    const bool half = argc == 2 && std::string_view(argv[1]) == "--half-step";
    require(argc == 1 || check || half,
            "usage: exaero_modal_mam4_comparison [--check|--half-step]");
    const double dt = half ? .5 : 1., end = check ? 10 : 600;
    std::cout
        << std::setprecision(17)
        << "experiment,case,time_s,field,unit,generic,mam4xx,chemistry_rhs\n";
    for (const auto &e : experiments)
      for (int c = 0; c < static_cast<int>(cases.size()); ++c) {
        Box generic(cases[c], e, false), native(cases[c], e, true);
        require(generic.state == native.state && generic.gas == native.gas,
                "comparison initial states differ");
        const auto initial = generic.budget();
        snapshot(e, c, 0, generic, native);
        for (int step = 1; step <= static_cast<int>(end / dt); ++step) {
          advance(generic, e, dt);
          advance(native, e, dt);
          const double time = step * dt;
          validate(generic, initial, e, time);
          validate(native, initial, e, time);
          if (e.native_coag_control || e.exponent == 1 || e.exponent == 2) {
            if (e.native_coag_control || (e.nucleation && !e.condensation)) {
              for (int i = 0; i < nf; ++i)
                require(std::abs(generic.state[i] - native.state[i]) <=
                            1e-10 * std::max(native.state[i], 1e-25),
                        "controlled matching-law comparison failed");
            }
          }
          if (time == end ||
              (!check && (time == 1 || time == 60 || time == 300)))
            snapshot(e, c, time, generic, native);
        }
      }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "FATAL ERROR: " << e.what() << '\n';
    return 1;
  }
}
