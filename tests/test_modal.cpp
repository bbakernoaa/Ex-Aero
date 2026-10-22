/// @file test_modal.cpp
/// @brief Modal conservation, analytic references, optics, and MKPP contracts.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: CTest; input: analytic fixtures; no output
/// files.
#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <cmath>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalCoupling_C.h>
#include <exaero/ModalMechanismMap.hpp>
#include <exaero/PackageFactory.hpp>
#include <fstream>
#include <limits>
#include <numeric>

namespace {
using namespace exaero;
constexpr double pi = 3.14159265358979323846;
/// @brief Deterministic mode configuration with explicit collision routes.
ModalConfig config(int nm = 1) {
  ModalConfig c;
  c.species = {{"acid", 1800, .098079, .6, {}}, {"organic", 1500, .15, .1, {}}};
  for (int m = 0; m < nm; ++m)
    c.modes.push_back({"mode" + std::to_string(m),
                       3e-8 * std::pow(3., m),
                       1.5,
                       {},
                       ModalOpticalTreatment::ComponentTables});
  for (int i = 0; i < nm; ++i)
    for (int j = i + 1; j < nm; ++j)
      c.coagulation_routes.push_back(
          {c.modes[i].name, c.modes[j].name, c.modes[j].name});
  return c;
}
/// @brief Multi-column host fixture with independent values at every level.
struct Fixture {
  static constexpr int nc = 2, nl = 3, np = nc * nl;
  ModalConfig cfg;
  std::unique_ptr<IAerosolPackage> package = createPackage("modal");
  std::array<double, np> t{}, p{}, rho{}, rh{}, dz{};
  std::vector<double> data;
  int nf, nm;
  explicit Fixture(ModalConfig c = config())
      : cfg(std::move(c)),
        data(np * cfg.modes.size() * (cfg.species.size() + 1)),
        nf(cfg.modes.size() * (cfg.species.size() + 1)), nm(cfg.modes.size()) {
    initializeModal(*package, cfg);
    auto s = state();
    for (int k = 0; k < nl; ++k)
      for (int cell = 0; cell < nc; ++cell) {
        const int point = cell + nc * k;
        t[point] = 280 + point;
        p[point] = 90000;
        rho[point] = 1.1;
        rh[point] = .2;
        dz[point] = 100 + point * 10;
        for (int m = 0; m < nm; ++m) {
          const std::array<double, 2> masses{7e-10 * (point + 1) * (m + 1),
                                             3e-10 * (point + 1) * (m + 1)};
          s(cell, k, m * 3) = masses[0];
          s(cell, k, m * 3 + 1) = masses[1];
          s(cell, k, m * 3 + 2) = coupling().emittedNumber(m, masses);
        }
      }
  }
  IModalCoupling &coupling() { return modalCoupling(*package); }
  EnvironmentalStateView env() const {
    return {View2D<const double>(t.data(), nc, nl),
            View2D<const double>(p.data(), nc, nl),
            View2D<const double>(rho.data(), nc, nl),
            View2D<const double>(rh.data(), nc, nl),
            View2D<const double>(dz.data(), nc, nl)};
  }
  View3D<double> state() { return View3D<double>(data.data(), nc, nl, nf); }
  View3D<const double> input() const {
    return View3D<const double>(data.data(), nc, nl, nf);
  }
  std::vector<double> properties() {
    std::vector<double> out(np * nm * modal_properties::COUNT);
    View4D<double> view(out.data(), nc, nl, nm, modal_properties::COUNT);
    coupling().modalProperties(env(), input(), view);
    return out;
  }
  void advance(double dt) {
    auto s = state();
    package->executeMicrophysics(env(), s, dt);
  }
  std::vector<double> optics(double wavelength = 5e-7) {
    std::vector<double> out(np * optical_indices::NUM_OPTICS);
    View4D<double> v(out.data(), nc, nl, 1, optical_indices::NUM_OPTICS);
    package->computeOptics(env(), input(), View1D<const double>(&wavelength, 1),
                           v);
    return out;
  }
};
/// @brief Test table semantics are explicitly declared, never inferred.
ModalTableBinding
binding(ModalTableKind kind = ModalTableKind::LognormalAverage) {
  ModalTableBinding b;
  b.kind = kind;
  b.geometric_std_dev = 1.5;
  b.curve.source_label = "MODAL_TEST";
  b.data_file = std::string(EXAERO_TEST_DATA_DIR) + "/mie_modal_reference.txt";
  b.wavelength_nodes = {5e-7, 1e-6};
  return b;
}
/// @brief Independent scalar Fuchs flux-matching kernel reference.
double referenceKernel(double r1, double r2, double density, double t,
                       double pressure) {
  const double kbt = 1.380649e-23 * t,
               mu = 1.716e-5 * std::pow(t / 273.15, 1.5) * 384.15 / (t + 111);
  std::array<double, 2> d{}, v{}, g{};
  const std::array<double, 2> r{r1, r2};
  for (int i = 0; i < 2; ++i) {
    const double kn = 6.65e-8 * (t / 293.15) * (101325 / pressure) / r[i];
    d[i] = kbt / (6 * pi * mu * r[i]) *
           (1 + kn * (1.257 + .4 * std::exp(-1.1 / kn)));
    v[i] =
        std::sqrt(8 * kbt / (pi * (4 * pi / 3 * std::pow(r[i], 3) * density)));
    const double path = 8 * d[i] / (pi * v[i]);
    g[i] = (std::pow(2 * r[i] + path, 3) -
            std::pow(4 * r[i] * r[i] + path * path, 1.5)) /
               (6 * r[i] * path) -
           2 * r[i];
  }
  const double rs = r1 + r2, ds = d[0] + d[1];
  return 4 * pi * rs * ds /
         (rs / (rs + std::hypot(g[0], g[1])) +
          4 * ds / (rs * std::hypot(v[0], v[1])));
}
} // namespace

TEST_CASE("generic mode count, metadata and geometry are configuration data") {
  for (int n : {1, 3, 4, 7}) {
    Fixture f(config(n));
    REQUIRE(f.package->representation() == AerosolRepresentation::Modal);
    REQUIRE(f.package->stateFields().size() == static_cast<std::size_t>(n * 3));
    const auto values = f.properties();
    View4D<const double> v(values.data(), f.nc, f.nl, n,
                           modal_properties::COUNT);
    for (int m = 0; m < n; ++m) {
      REQUIRE(f.package->getSpeciesIndex("num_mode" + std::to_string(m) +
                                         "_i") == 3 * m + 2);
      REQUIRE(v(1, 2, m, modal_properties::DRY_RADIUS) ==
              Approx(f.cfg.modes[m].geometric_radius));
      REQUIRE(v(1, 2, m, modal_properties::WET_RADIUS) >
              v(1, 2, m, modal_properties::DRY_RADIUS));
    }
  }
}

TEST_CASE(
    "YAML and structured modal initialization agree and fail transactionally") {
  Fixture f;
  auto p = createPackage("modal");
  p->initialize(R"(backend: modal
species:
  - {name: acid, density: 1800, molecular_weight: 0.098079, hygroscopicity: 0.6}
  - {name: organic, density: 1500, molecular_weight: 0.15, hygroscopicity: 0.1}
modes:
  - {name: mode0, geometric_radius: 3e-8, geometric_std_dev: 1.5}
)");
  REQUIRE(p->getSpeciesName(1) == f.package->getSpeciesName(1));
  const auto fields = p->stateFields();
  REQUIRE_THROWS(p->initialize("backend: modal\nunknown: true"));
  REQUIRE(p->stateFields().size() == fields.size());
  auto c = config(3);
  c.coagulation_routes.clear();
  REQUIRE_THROWS(initializeModal(*p, c));
  c = config(3);
  c.transfer_routes = {{"mode0", "mode1", 1e-7}, {"mode1", "mode0", 1e-7}};
  REQUIRE_THROWS(initializeModal(*p, c));
  c = config();
  c.modes[0].geometric_std_dev = 1;
  REQUIRE_THROWS(initializeModal(*p, c));
  c = config();
  c.species.push_back(c.species[0]);
  REQUIRE_THROWS(initializeModal(*p, c));
  c = config(2);
  c.modes[1].species = {"organic"};
  REQUIRE_THROWS(initializeModal(*p, c));
}

TEST_CASE("modal properties and uptake provide analytic derivatives") {
  Fixture f;
  std::vector<double> props(f.np * 9), jac(f.np * 9 * f.nf);
  View4D<double> pv(props.data(), f.nc, f.nl, 1, 9),
      jv(jac.data(), f.nc, f.nl, 1, 9 * f.nf);
  f.coupling().modalProperties(f.env(), f.input(), pv, &jv);
  const std::array<ModalVapor, 1> gases{{{1e-5, .098079, .7}}};
  std::vector<double> uptake(f.np), ujac(f.np * f.nf);
  View4D<double> uv(uptake.data(), f.nc, f.nl, 1, 1),
      uvj(ujac.data(), f.nc, f.nl, 1, f.nf);
  f.coupling().modalUptake(f.env(), f.input(), gases, uv, &uvj);
  REQUIRE(uv(0, 0, 0, 0) > 0);
  for (int field = 0; field < f.nf; ++field) {
    auto state = f.state();
    const double original = state(0, 0, field), step = original * 1e-5;
    state(0, 0, field) = original + step;
    const auto plus = f.properties();
    std::vector<double> up(f.np);
    View4D<double> upv(up.data(), f.nc, f.nl, 1, 1);
    f.coupling().modalUptake(f.env(), f.input(), gases, upv);
    state(0, 0, field) = original - step;
    const auto minus = f.properties();
    std::vector<double> um(f.np);
    View4D<double> umv(um.data(), f.nc, f.nl, 1, 1);
    f.coupling().modalUptake(f.env(), f.input(), gases, umv);
    state(0, 0, field) = original;
    for (int q = 0; q < 9; ++q) {
      const double fd = (plus[f.np * q] - minus[f.np * q]) / (2 * step);
      REQUIRE(jv(0, 0, 0, q * f.nf + field) ==
              Approx(fd).epsilon(2e-5).margin(1e-14));
    }
    const double fd = (up[0] - um[0]) / (2 * step);
    REQUIRE(uvj(0, 0, 0, field) == Approx(fd).epsilon(2e-5));
  }
}

TEST_CASE(
    "single-mode coagulation matches independent distribution integration") {
  auto cfg = config();
  cfg.species[0].density = cfg.species[1].density = 1800;
  Fixture f(cfg);
  f.rh.fill(0);
  auto state = f.state();
  const double n = state(0, 0, 2), rg = cfg.modes[0].geometric_radius;
  double beta = 0;
  constexpr int bins = 180;
  const double dx = 16. / bins;
  for (int i = 0; i < bins; ++i)
    for (int j = 0; j < bins; ++j) {
      const double x = -8 + (i + .5) * dx, y = -8 + (j + .5) * dx;
      beta += dx * dx / (2 * pi) * std::exp(-.5 * (x * x + y * y)) *
              referenceKernel(rg * std::pow(1.5, x), rg * std::pow(1.5, y),
                              1800, 280, 90000);
    }
  f.advance(.1);
  const double numerical = (n - state(0, 0, 2)) / .1;
  REQUIRE(numerical == Approx(.5 * beta * n * n).epsilon(.003));
}

TEST_CASE(
    "coagulation conserves every material and populates an empty destination") {
  auto cfg = config(3);
  for (auto &r : cfg.coagulation_routes)
    r.destination = "mode2";
  Fixture f(cfg);
  auto s = f.state();
  for (int k = 0; k < f.nl; ++k)
    for (int c = 0; c < f.nc; ++c)
      for (int q = 0; q < 3; ++q)
        s(c, k, 6 + q) = 0;
  const auto before = f.data;
  f.advance(60);
  View3D<const double> old(before.data(), f.nc, f.nl, f.nf);
  for (int k = 0; k < f.nl; ++k)
    for (int c = 0; c < f.nc; ++c) {
      for (int species = 0; species < 2; ++species) {
        double a = 0, b = 0;
        for (int m = 0; m < 3; ++m) {
          a += old(c, k, 3 * m + species);
          b += s(c, k, 3 * m + species);
        }
        REQUIRE(b == Approx(a).epsilon(1e-12));
      }
      REQUIRE(s(c, k, 8) > 0);
      REQUIRE(s(c, k, 2) + s(c, k, 5) + s(c, k, 8) <
              old(c, k, 2) + old(c, k, 5));
    }
  for (double a : f.data)
    REQUIRE(a >= 0);
}

TEST_CASE("coagulation converges with timestep refinement") {
  Fixture coarse, medium, fine;
  coarse.advance(400);
  for (int i = 0; i < 8; ++i)
    medium.advance(50);
  for (int i = 0; i < 80; ++i)
    fine.advance(5);
  const double truth = fine.state()(0, 0, 2);
  REQUIRE(std::abs(medium.state()(0, 0, 2) - truth) <
          std::abs(coarse.state()(0, 0, 2) - truth));
}

TEST_CASE("growth transfer uses separate number and volume tails and preserves "
          "totals") {
  auto cfg = config(3);
  cfg.coagulation = false;
  cfg.transfer_routes = {
      {"mode1", "mode2", 1e-7},
      {"mode0", "mode1", 2e-8}}; // deliberately reverse order
  Fixture f(cfg);
  const auto before = f.data;
  f.advance(1);
  const auto s = f.state();
  View3D<const double> old(before.data(), f.nc, f.nl, f.nf);
  const double ln = std::log(1.5), x = std::log(2e-8 / 3e-8);
  const double fn = .5 * std::erfc(x / (std::sqrt(2.) * ln)),
               fm = .5 * std::erfc((x - 3 * ln * ln) / (std::sqrt(2.) * ln));
  REQUIRE(s(0, 0, 2) == Approx(old(0, 0, 2) * (1 - fn)));
  REQUIRE(s(0, 0, 0) == Approx(old(0, 0, 0) * (1 - fm)));
  for (int q = 0; q < 3; ++q)
    REQUIRE(s(0, 0, q) + s(0, 0, 3 + q) + s(0, 0, 6 + q) ==
            Approx(old(0, 0, q) + old(0, 0, 3 + q) + old(0, 0, 6 + q)));
}

TEST_CASE(
    "invalid inputs leave state and output unchanged; zero step is exact") {
  Fixture f;
  const auto before = f.data;
  f.advance(0);
  REQUIRE(f.data == before);
  f.rh[0] = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS(f.advance(1));
  REQUIRE(f.data == before);
  f.rh[0] = .2;
  f.state()(0, 0, 2) = 0;
  const auto bad = f.data;
  REQUIRE_THROWS(f.advance(1));
  REQUIRE(f.data == bad);
  f.data = before;
  f.state()(0, 0, 0) = -1;
  REQUIRE_THROWS(f.properties());
  f.data.assign(f.data.size(), 0);
  const auto empty = f.properties();
  REQUIRE(empty[0] == Approx(f.cfg.modes[0].geometric_radius));
}

TEST_CASE("diagnostics, CCN, and emissions respect modal concentration units") {
  Fixture f;
  std::vector<double> diag(f.np * diagnostic_indices::NUM_DIAGNOSTICS);
  View3D<double> dv(diag.data(), f.nc, f.nl,
                    diagnostic_indices::NUM_DIAGNOSTICS);
  f.package->computeDerivedDiagnostics(f.env(), f.input(), dv);
  REQUIRE(dv(0, 0, diagnostic_indices::MASS_CONCENTRATION) == Approx(1e-9));
  REQUIRE(dv(0, 0, diagnostic_indices::AEROSOL_LIQUID_WATER) > 0);
  double column = 0;
  for (int k = 0; k < f.nl; ++k)
    column += dv(0, k, 0) * f.dz[f.nc * k];
  REQUIRE(dv(0, 0, diagnostic_indices::COLUMN_MASS) == Approx(column));
  const std::array<double, 4> ss{0, .001, .01, 1};
  std::vector<double> ccn(f.np * ss.size());
  View4D<double> cv(ccn.data(), f.nc, f.nl, ss.size(), 1);
  f.package->computeCCN(f.env(), f.input(),
                        View1D<const double>(ss.data(), ss.size()), cv);
  REQUIRE(cv(0, 0, 0, 0) == 0);
  for (int b = 1; b < 4; ++b) {
    REQUIRE(cv(0, 0, b, 0) >= cv(0, 0, b - 1, 0));
    REQUIRE(cv(0, 0, b, 0) <= f.state()(0, 0, 2));
  }
  std::vector<double> rates(f.data.size());
  View3D<double> rv(rates.data(), f.nc, f.nl, f.nf);
  f.package->computeEmissions(f.env(), {f.input(), FluxType::AREA_FLUX}, rv);
  REQUIRE(rv(0, 0, 2) == Approx(f.state()(0, 0, 2) / f.dz[0]));
  f.dz[0] = 0;
  const auto saved = rates;
  REQUIRE_THROWS(f.package->computeEmissions(
      f.env(), {f.input(), FluxType::AREA_FLUX}, rv));
  REQUIRE(rates == saved);
}

TEST_CASE("oversized emission extents fail before touching host memory") {
  Fixture f;
  const auto before = f.data;
  const auto huge = static_cast<std::size_t>(std::numeric_limits<int>::max());
  View3D<const double> flux(f.data.data(), huge, 2, f.nf);
  View3D<double> output(f.data.data(), huge, 2, f.nf);
  REQUIRE_THROWS(f.package->computeEmissions(
      f.env(), {flux, FluxType::MASS_CONCENTRATION_RATE}, output));
  REQUIRE(f.data == before);
}

TEST_CASE("MKPP mappings distinguish material molecules from particle number") {
  auto cfg = config();
  cfg.nucleation =
      ModalNucleationConfig{"mode0", 1e-12, 1, 1e-9, {1, 0}, .098079};
  Fixture f(cfg);
  const std::array<std::string, 4> mechanism{"H2SO4_g", "AERO_ORG", "N_MODE",
                                             "AERO_ACID"};
  const std::array<ModalMechanismBinding, 3> bindings{
      {{"acid_mode0_i", "AERO_ACID"},
       {"organic_mode0_i", "AERO_ORG"},
       {"num_mode0_i", "N_MODE"}}};
  const auto map = modalMechanismMap(*f.package, mechanism, bindings);
  REQUIRE(map[0].mechanism_index == 3);
  REQUIRE(map[0].to_mechanism == Approx(6.02214076e17 / .098079));
  REQUIRE(map[2].to_mechanism == 1e-6);
  REQUIRE(map[2].quantity == StateQuantity::Number);
  auto bad = bindings;
  bad[2].mechanism_variable = "AERO_ACID";
  REQUIRE_THROWS(modalMechanismMap(*f.package, mechanism, bad));
  const double gas = 1e-12;
  const auto terms = f.coupling().nucleationTerms(gas);
  REQUIRE(terms.number_rate == Approx(1e-12 * gas / .098079 * 6.02214076e23));
  REQUIRE(terms.precursor_loss == Approx(terms.mass_rates[0]));
  REQUIRE(terms.mass_rates[1] == 0);
  // Mock chemistry owner applies matching gas loss/material gain once.
  const double aerosol = f.state()(0, 0, 0), number = f.state()(0, 0, 2),
               dt = 1;
  f.state()(0, 0, 0) += terms.mass_rates[0] * dt;
  f.state()(0, 0, 2) += terms.number_rate * dt;
  REQUIRE(gas - terms.precursor_loss * dt + f.state()(0, 0, 0) ==
          Approx(gas + aerosol));
  REQUIRE(f.state()(0, 0, 2) > number);
  // Condensation exchanges mass, retaining independent number.
  const auto n = f.state()(0, 0, 2);
  const auto radius = f.properties()[0];
  f.state()(0, 0, 0) += 1e-10;
  REQUIRE(f.state()(0, 0, 2) == n);
  REQUIRE(f.properties()[0] > radius);
}

TEST_CASE(
    "component optics reuse source coefficients and scattering-weight g") {
  auto cfg = config();
  cfg.species[0].optics = binding();
  cfg.species[1].optics = binding();
  auto &b = cfg.species[1].optics->curve;
  b.radius_nodes = {1e-10, 1e-4};
  b.overrides = {{AttributeCategory::SpectralOptical,
                  spectral_optical_indices::MASS_EXTINCTION,
                  {500, 500}},
                 {AttributeCategory::SpectralOptical,
                  spectral_optical_indices::MASS_SCATTERING,
                  {100, 100}},
                 {AttributeCategory::SpectralOptical,
                  spectral_optical_indices::ASYMMETRY_FACTOR,
                  {.8, .8}}};
  Fixture f(cfg);
  f.rh.fill(0);
  REQUIRE(f.package->capabilities() & Optics);
  const auto out = f.optics();
  View4D<const double> ov(out.data(), f.nc, f.nl, 1,
                          optical_indices::NUM_OPTICS);
  const double ext = 7e-10 * 1000 + 3e-10 * 500,
               sca = 7e-10 * 500 + 3e-10 * 100;
  REQUIRE(ov(0, 0, 0, optical_indices::EXTINCTION_COEFF) == Approx(ext));
  REQUIRE(ov(0, 0, 0, optical_indices::ASYMMETRY_FACTOR) ==
          Approx((7e-10 * 500 * .3 + 3e-10 * 100 * .8) / sca));
  f.rh.fill(.98);
  const auto wet = f.optics();
  REQUIRE(wet[0] == Approx(7e-10 * 2000 + 3e-10 * 500));
  // A second package initialization cannot alter this package's table snapshot.
  Fixture unrelated;
  REQUIRE(f.optics()[0] == Approx(wet[0]));
  const auto before = wet;
  REQUIRE_THROWS(f.optics(2e-6));
}

TEST_CASE("particle-resolved optical integration matches lognormal analytic "
          "moments") {
  auto cfg = config();
  for (auto &s : cfg.species) {
    s.optics = binding(ModalTableKind::ParticleResolved);
    s.optics->curve.radius_nodes = {1e-10, 1e-4};
    s.optics->curve.overrides = {{AttributeCategory::SpectralOptical,
                                  spectral_optical_indices::MASS_EXTINCTION,
                                  {1, 1e6}},
                                 {AttributeCategory::SpectralOptical,
                                  spectral_optical_indices::MASS_SCATTERING,
                                  {.5, 5e5}}};
  }
  Fixture f(cfg);
  f.rh.fill(0);
  const auto out = f.optics();
  const double ln = std::log(1.5), expected = 1e-9 * 1e10 *
                                              cfg.modes[0].geometric_radius *
                                              std::exp(3.5 * ln * ln);
  REQUIRE(out[0] == Approx(expected).epsilon(2e-4));
}

TEST_CASE("unknown optical semantics and mismatched distribution widths fail") {
  auto cfg = config();
  cfg.species[0].optics = binding(ModalTableKind::Unknown);
  auto package = createPackage("modal");
  REQUIRE_THROWS(initializeModal(*package, cfg));
  cfg.species[0].optics = binding();
  cfg.species[0].optics->geometric_std_dev = 1.8;
  REQUIRE_THROWS(initializeModal(*package, cfg));
  cfg.species[0].optics = binding();
  cfg.species[0].optics->wavelength_nodes.clear();
  REQUIRE_THROWS(initializeModal(*package, cfg));
  Fixture no_tables;
  REQUIRE(!(no_tables.package->capabilities() & Optics));
  REQUIRE_THROWS(no_tables.optics());
}

TEST_CASE(
    "baked GEOSmie curves carry known semantics and support evolving spreads") {
  auto cfg = config();
  cfg.species[0].density = 1700;
  cfg.species[0].hygroscopicity = .507;
  cfg.species[1].density = 1800;
  cfg.species[1].hygroscopicity = .14;
  cfg.species[0].optics = ModalTableBinding{};
  cfg.species[0].optics->curve.source_label = "MODAL_SU";
  cfg.species[1].optics = ModalTableBinding{};
  cfg.species[1].optics->curve.source_label = "MODAL_OC";
  cfg.modes[0].geometric_std_dev = 1.6;
  Fixture f(cfg);
  REQUIRE(f.package->capabilities() & Optics);
  const auto before = f.optics(5.5e-7);
  REQUIRE(before[0] > 0);
  REQUIRE(before[f.np] <= before[0]);
  for (double &x : f.data)
    x *= 2;
  const auto doubled = f.optics(5.5e-7);
  REQUIRE(doubled[0] == Approx(2 * before[0]).epsilon(1e-10));
  REQUIRE(doubled[3 * f.np] == Approx(before[3 * f.np]));
  auto bad = cfg;
  bad.species[0].density = 999;
  REQUIRE_THROWS(initializeModal(*f.package, bad));
  bad = cfg;
  bad.modes[0].geometric_radius = 1e-4;
  Fixture outside(bad);
  REQUIRE_THROWS(outside.optics(5.5e-7));
}

TEST_CASE("documented examples compute optics with dust and ammonium nitrate") {
  for (const char *file : {"modal_three_modes.yaml", "modal_four_modes.yaml"}) {
    auto package = createPackage("modal");
    std::ifstream stream(std::string(EXAERO_TEST_SOURCE_DIR) + "/examples/" +
                         file);
    REQUIRE(stream.good());
    package->initialize(
        std::string(std::istreambuf_iterator<char>(stream), {}));
    auto &coupling = modalCoupling(*package);
    const auto &cfg = coupling.modalConfig();
    REQUIRE(cfg.species.size() == 6);
    REQUIRE(package->capabilities() & Optics);
    const auto sources = coupling.opticalTables();
    REQUIRE(std::any_of(sources.begin(), sources.end(), [](const auto &s) {
      return s.source_label == "MODAL_DU";
    }));
    REQUIRE(std::any_of(sources.begin(), sources.end(), [](const auto &s) {
      return s.source_label == "MODAL_NI";
    }));
    for (const auto &s : sources)
      REQUIRE(s.source_version.find("c535388") != std::string::npos);
    const auto ns = cfg.species.size(), nf = package->stateFields().size();
    std::vector<double> state(nf), mass(ns, 1e-10);
    for (std::size_t m = 0; m < cfg.modes.size(); ++m) {
      for (std::size_t s = 0; s < ns; ++s)
        state[m * (ns + 1) + s] = mass[s];
      state[m * (ns + 1) + ns] = coupling.emittedNumber(m, mass);
    }
    const double t = 280, p = 90000, density = 1.1, rh = .8, dz = 100,
                 wavelength = 5.5e-7;
    EnvironmentalStateView env{
        View2D<const double>(&t, 1, 1), View2D<const double>(&p, 1, 1),
        View2D<const double>(&density, 1, 1), View2D<const double>(&rh, 1, 1),
        View2D<const double>(&dz, 1, 1)};
    std::array<double, optical_indices::NUM_OPTICS> optics{};
    View4D<double> out(optics.data(), 1, 1, 1, optics.size());
    package->computeOptics(env, View3D<const double>(state.data(), 1, 1, nf),
                           View1D<const double>(&wavelength, 1), out);
    REQUIRE(out(0, 0, 0, 0) > 0);
    REQUIRE(out(0, 0, 0, 1) <= out(0, 0, 0, 0));
  }
}

namespace {
/// @brief Supplied mixed-particle provider, with explicit domain validation.
class Provider final : public IModalOpticsProvider {
public:
  ModalOpticalResult evaluate(const ModalConfig &, std::size_t, double, double,
                              double wavelength, std::span<const double> masses,
                              std::span<const double>) const override {
    if (wavelength > 1e-6)
      throw std::invalid_argument("FATAL ERROR: provider wavelength domain");
    const double e = 100 * std::accumulate(masses.begin(), masses.end(), 0.);
    return {e, e * .5, .6, e * .01, e * .02, e * .8, e * .4, e * .9, e * .45};
  }
};
} // namespace
TEST_CASE("internally mixed optics require a supplied compatible provider") {
  auto cfg = config();
  cfg.modes[0].optical_treatment = ModalOpticalTreatment::InternalMixture;
  Fixture f(cfg);
  REQUIRE_THROWS(f.optics());
  f.coupling().setModalOpticsProvider(std::make_shared<Provider>());
  REQUIRE(f.package->capabilities() & Optics);
  REQUIRE(f.optics()[0] == Approx(1e-7));
  REQUIRE_THROWS(f.optics(2e-6));
  f.coupling().setModalOpticsProvider(nullptr);
  REQUIRE(!(f.package->capabilities() & Optics));
}

TEST_CASE("C coupling ABI checks dimensions, units and failure atomicity") {
  Fixture f;
  auto handle = static_cast<exaero_package_t>(f.package.get());
  char error[256]{};
  int status = -1;
  REQUIRE(exaero_get_modal_mode_count(handle) == 1);
  REQUIRE(exaero_get_modal_mode_count(nullptr) == -1);
  std::vector<double> out(f.np * 9, -7), jac(f.np * 9 * f.nf);
  exaero_modal_properties(handle, f.nc, f.nl, f.nf, 1, f.t.data(), f.rh.data(),
                          f.data.data(), out.data(), jac.data(), error,
                          &status);
  REQUIRE(status == 0);
  REQUIRE(out[0] == Approx(3e-8));
  const auto saved = out;
  exaero_modal_properties(handle, f.nc, f.nl, f.nf, 2, f.t.data(), f.rh.data(),
                          f.data.data(), out.data(), nullptr, error, &status);
  REQUIRE(status == 1);
  REQUIRE(out == saved);
  exaero_modal_properties(handle, f.nc, f.nl, f.nf, 1, nullptr, f.rh.data(),
                          f.data.data(), out.data(), nullptr, error, &status);
  REQUIRE(status == 1);
  REQUIRE(out == saved);
  double number = -1;
  const std::array<double, 2> mass{7e-10, 3e-10};
  exaero_modal_emitted_number(handle, 0, 2, mass.data(), &number, error,
                              &status);
  REQUIRE(status == 0);
  REQUIRE(number == Approx(f.state()(0, 0, 2)));
}
