/// @file test_modal_chemistry.cpp
/// @brief Independent budgets, Jacobian finite differences and coupling
/// lifecycle.
/// @author EX-aero contributors
/// @details History: 2026-10. Multi-column host views; Kokkos lifetime in main
/// TU.
#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <cmath>
#include <exaero/ModalChemistry.hpp>
#include <exaero/ModalChemistry_C.h>
#include <exaero/ModalMechanismMap.hpp>
#include <exaero/ModalPresets.hpp>
#include <exaero/PackageFactory.hpp>
#include <limits>
#include <numeric>

namespace {
using namespace exaero;
/// @brief Two-material two-mode configuration with explicit coating-shell
/// closure.
ModalConfig chemistryConfig() {
  ModalConfig c;
  c.species = {{"acid", 1700, .098079, .5, {}}, {"organic", 1000, .15, .1, {}}};
  c.modes = {{"fresh", 3e-8, 1.6, {}}, {"aged", 1e-7, 1.8, {}}};
  c.coagulation_routes = {{"fresh", "aged", "aged"}};
  c.aging_routes = {{"fresh", "aged", {"acid"}, 1e-9}};
  c.transfer_owner = ModalTransferOwner::Host;
  return c;
}
/// @brief Nonuniform host fixture, flat vectors interpreted exclusively through
/// mdspan.
struct ChemistryFixture {
  static constexpr int nc = 2, nl = 3, np = nc * nl;
  ModalConfig config;
  int ns, nm, nf;
  std::unique_ptr<IAerosolPackage> package = createPackage("modal");
  std::array<double, np> t{}, p{}, rho{}, rh{}, dz{};
  std::vector<double> data;
  explicit ChemistryFixture(ModalConfig c = chemistryConfig())
      : config(std::move(c)), ns(config.species.size()),
        nm(config.modes.size()), nf(nm * (ns + 1)), data(np * nf) {
    initializeModal(*package, config);
    auto state = view();
    for (int k = 0; k < nl; ++k)
      for (int cell = 0; cell < nc; ++cell) {
        const int point = cell + nc * k;
        t[point] = 280 + point;
        p[point] = 90000 - 1000 * point;
        rho[point] = 1 + .1 * point;
        rh[point] = .5;
        dz[point] = 100 + 10 * point;
        for (int m = 0; m < nm; ++m) {
          std::vector<double> mass(ns);
          for (int s = 0; s < ns; ++s)
            state(cell, k, m * (ns + 1) + s) = mass[s] =
                2e-10 * (point + 1) * (m + 1) * (s + 1);
          state(cell, k, m * (ns + 1) + ns) =
              modalCoupling(*package).emittedNumber(m, mass);
        }
      }
  }
  /// @brief Host environmental views.
  EnvironmentalStateView env() const {
    return {View2D<const double>(t.data(), nc, nl),
            View2D<const double>(p.data(), nc, nl),
            View2D<const double>(rho.data(), nc, nl),
            View2D<const double>(rh.data(), nc, nl),
            View2D<const double>(dz.data(), nc, nl)};
  }
  /// @brief Mutable stage state.
  View3D<double> view() { return View3D<double>(data.data(), nc, nl, nf); }
  /// @brief Immutable stage state.
  View3D<const double> input() const {
    return View3D<const double>(data.data(), nc, nl, nf);
  }
};
/// @brief Coupled exchange/aging RHS with explicit linear test thermochemistry.
struct ChemicalRHS {
  std::vector<double> aerosol, gas, jac, aging, aging_jac;
};
/// @brief Evaluate exchange plus aging, including complete
/// thermochemical/budget derivatives.
ChemicalRHS evaluate(ChemistryFixture &f, const std::vector<double> &gas,
                     bool jacobian = true) {
  constexpr int nv = 2;
  const int nq = f.nf + nv;
  ChemicalRHS out{std::vector<double>(f.np * f.nf),
                  std::vector<double>(f.np * nv),
                  std::vector<double>(jacobian ? f.np * nq * nq : 0),
                  std::vector<double>(f.np * f.nf),
                  std::vector<double>(jacobian ? f.np * f.nf * nq : 0)};
  const std::array<ModalCondensable, 2> bindings{
      {{"acid", {1e-5, .098079, 1}}, {"organic", {8e-6, .15, .7}}}};
  std::vector<double> equilibrium(f.np * f.nm * nv), de(f.np * f.nm * nv * nq),
      coag(f.np * f.nf), dg(f.np * f.nf * nq);
  View4D<double> ce(equilibrium.data(), f.nc, f.nl, f.nm, nv),
      eqjac(de.data(), f.nc, f.nl, f.nm * nv, nq),
      gj(dg.data(), f.nc, f.nl, f.nf, nq);
  View3D<double> budget(coag.data(), f.nc, f.nl, f.nf);
  const View3D<const double> g(gas.data(), f.nc, f.nl, nv);
  for (int k = 0; k < f.nl; ++k)
    for (int c = 0; c < f.nc; ++c) {
      for (int m = 0; m < f.nm; ++m)
        for (int v = 0; v < nv; ++v) {
          const int field = m * (f.ns + 1) + v;
          ce(c, k, m, v) = (m + 1) * 1e-12 + .001 * f.input()(c, k, field) +
                           .1 * g(c, k, 1 - v);
          eqjac(c, k, m * nv + v, field) = .001;
          eqjac(c, k, m * nv + v, f.nf + 1 - v) = .1;
        }
      budget(c, k, 0) = 1e-4 * f.input()(c, k, 0) + .001 * g(c, k, 0);
      gj(c, k, 0, 0) = 1e-4;
      gj(c, k, 0, f.nf) = .001;
    }
  View3D<double> ar(out.aerosol.data(), f.nc, f.nl, f.nf),
      gr(out.gas.data(), f.nc, f.nl, nv),
      age(out.aging.data(), f.nc, f.nl, f.nf);
  View4D<double> j(out.jac.data(), f.nc, f.nl, nq, nq),
      aj(out.aging_jac.data(), f.nc, f.nl, f.nf, nq);
  const View4D<const double> dej(de.data(), f.nc, f.nl, f.nm * nv, nq);
  modalPhaseExchange(
      *f.package, f.env(), f.input(), bindings, g,
      View4D<const double>(equilibrium.data(), f.nc, f.nl, f.nm, nv), ar, gr,
      jacobian ? &j : nullptr, jacobian ? &dej : nullptr);
  const View4D<const double> cj(out.jac.data(), f.nc, f.nl, nq, nq),
      coagj(dg.data(), f.nc, f.nl, f.nf, nq);
  modalAgingTerms(*f.package, f.env(), f.input(),
                  View3D<const double>(out.aerosol.data(), f.nc, f.nl, f.nf),
                  View3D<const double>(coag.data(), f.nc, f.nl, f.nf), age,
                  jacobian ? &aj : nullptr, jacobian ? &cj : nullptr,
                  jacobian ? &coagj : nullptr);
  return out;
}
} // namespace

TEST_CASE("phase exchange and coating aging conserve every material and "
          "particle number") {
  ChemistryFixture f;
  std::vector<double> gas(f.np * 2, 4e-12);
  const auto before = f.data;
  const auto out = evaluate(f, gas);
  REQUIRE(f.data == before);
  const View3D<const double> ar(out.aerosol.data(), f.nc, f.nl, f.nf),
      gr(out.gas.data(), f.nc, f.nl, 2),
      age(out.aging.data(), f.nc, f.nl, f.nf);
  const int nq = f.nf + 2;
  const View4D<const double> j(out.jac.data(), f.nc, f.nl, nq, nq),
      aj(out.aging_jac.data(), f.nc, f.nl, f.nf, nq);
  for (int k = 0; k < f.nl; ++k)
    for (int c = 0; c < f.nc; ++c) {
      for (int s = 0; s < 2; ++s) {
        REQUIRE(ar(c, k, s) + ar(c, k, 3 + s) + gr(c, k, s) ==
                Approx(0).margin(1e-25));
        REQUIRE(age(c, k, s) + age(c, k, 3 + s) == Approx(0).margin(1e-25));
        for (int col = 0; col < nq; ++col) {
          REQUIRE(j(c, k, s, col) + j(c, k, 3 + s, col) +
                      j(c, k, f.nf + s, col) ==
                  Approx(0).margin(1e-12));
          REQUIRE(aj(c, k, s, col) + aj(c, k, 3 + s, col) ==
                  Approx(0).margin(1e-12));
        }
      }
      REQUIRE(ar(c, k, 2) == 0);
      REQUIRE(ar(c, k, 5) == 0);
      REQUIRE(age(c, k, 2) < 0);
      REQUIRE(age(c, k, 2) + age(c, k, 5) == 0);
      const double n = f.input()(c, k, 2),
                   l = std::log(f.config.modes[0].geometric_std_dev),
                   r = f.config.modes[0].geometric_radius, h = 1e-9;
      const double shell = 4 * std::acos(-1.) * n *
                           (r * r * std::exp(2 * l * l) * h +
                            r * std::exp(.5 * l * l) * h * h + h * h * h / 3);
      const double budget =
          (ar(c, k, 0) + 1e-4 * f.input()(c, k, 0) + .001 * gas[c + f.nc * k]) /
          1700;
      REQUIRE(-age(c, k, 2) / n == Approx(budget / shell).epsilon(1e-12));
    }
}

TEST_CASE("full exchange and aging Jacobians match finite differences "
          "including thermochemistry") {
  ChemistryFixture f;
  std::vector<double> gas(f.np * 2, 4e-12);
  const int nq = f.nf + 2;
  const auto baseline = evaluate(f, gas);
  const View4D<const double> j(baseline.jac.data(), f.nc, f.nl, nq, nq),
      aj(baseline.aging_jac.data(), f.nc, f.nl, f.nf, nq);
  for (int column = 0; column < nq; ++column) {
    double &value =
        column < f.nf ? f.data[f.np * column] : gas[f.np * (column - f.nf)];
    const double saved = value, step = saved * 1e-5;
    value = saved + step;
    const auto plus = evaluate(f, gas, false);
    value = saved - step;
    const auto minus = evaluate(f, gas, false);
    value = saved;
    for (int row = 0; row < nq; ++row) {
      const double p =
          row < f.nf ? plus.aerosol[f.np * row] : plus.gas[f.np * (row - f.nf)];
      const double m = row < f.nf ? minus.aerosol[f.np * row]
                                  : minus.gas[f.np * (row - f.nf)];
      INFO(row);
      INFO(column);
      REQUIRE(j(0, 0, row, column) ==
              Approx((p - m) / (2 * step)).epsilon(2e-6).margin(1e-22));
      if (row < f.nf)
        REQUIRE(aj(0, 0, row, column) ==
                Approx((plus.aging[f.np * row] - minus.aging[f.np * row]) /
                       (2 * step))
                    .epsilon(2e-6)
                    .margin(1e-22));
    }
  }
}

TEST_CASE(
    "evaporation is blocked for absent material and outputs fail atomically") {
  ChemistryFixture f;
  auto state = f.view();
  state(0, 0, 0) = 0;
  std::array<ModalCondensable, 1> binding{{{"acid", {1e-5, .098079, 1}}}};
  std::vector<double> gas(f.np, 0), eq(f.np * f.nm, 1e-12), ar(f.np * f.nf, -1),
      gr(f.np, -1), jac(f.np * (f.nf + 1) * (f.nf + 1), -1);
  View3D<double> a(ar.data(), f.nc, f.nl, f.nf), g(gr.data(), f.nc, f.nl, 1);
  View4D<double> j(jac.data(), f.nc, f.nl, f.nf + 1, f.nf + 1);
  modalPhaseExchange(*f.package, f.env(), f.input(), binding,
                     View3D<const double>(gas.data(), f.nc, f.nl, 1),
                     View4D<const double>(eq.data(), f.nc, f.nl, f.nm, 1), a, g,
                     &j);
  REQUIRE(a(0, 0, 0) == 0);
  for (int column = 0; column <= f.nf; ++column)
    REQUIRE(j(0, 0, 0, column) == 0);
  REQUIRE(a(0, 0, 3) < 0); // Existing aged material can evaporate.
  const auto saved_ar = ar, saved_gr = gr, saved_jac = jac;
  gas[0] = std::numeric_limits<double>::quiet_NaN();
  REQUIRE_THROWS(modalPhaseExchange(
      *f.package, f.env(), f.input(), binding,
      View3D<const double>(gas.data(), f.nc, f.nl, 1),
      View4D<const double>(eq.data(), f.nc, f.nl, f.nm, 1), a, g, &j));
  REQUIRE(ar == saved_ar);
  REQUIRE(gr == saved_gr);
  REQUIRE(jac == saved_jac);
  gas[0] = 0;
  View3D<double> alias(f.data.data(), f.nc, f.nl, f.nf);
  REQUIRE_THROWS(modalPhaseExchange(
      *f.package, f.env(), f.input(), binding,
      View3D<const double>(gas.data(), f.nc, f.nl, 1),
      View4D<const double>(eq.data(), f.nc, f.nl, f.nm, 1), alias, g));
  auto bad = f.config;
  bad.aging_routes[0].destination = "fresh";
  REQUIRE_THROWS(initializeModal(*f.package, bad));
  bad = f.config;
  bad.aging_routes.push_back({"aged", "fresh", {"acid"}, 1e-9});
  REQUIRE_THROWS(initializeModal(*f.package, bad));
  std::vector<double> cond(f.data.size()), coag(f.data.size()),
      aging(f.data.size(), -9);
  coag[f.np * 2] = 1;
  View3D<double> output(aging.data(), f.nc, f.nl, f.nf);
  REQUIRE_THROWS(modalAgingTerms(
      *f.package, f.env(), f.input(),
      View3D<const double>(cond.data(), f.nc, f.nl, f.nf),
      View3D<const double>(coag.data(), f.nc, f.nl, f.nf), output));
  REQUIRE(std::all_of(aging.begin(), aging.end(),
                      [](double x) { return x == -9; }));
}

TEST_CASE("gross coagulation coating gains are integrated once and renaming "
          "has one owner") {
  ChemistryFixture f;
  const auto before = f.data;
  std::vector<double> gain(f.data.size(), -1);
  auto state = f.view();
  View3D<double> out(gain.data(), f.nc, f.nl, f.nf);
  advanceModalWithCoating(*f.package, f.env(), state, 100, out);
  for (int k = 0; k < f.nl; ++k)
    for (int c = 0; c < f.nc; ++c)
      for (int s = 0; s < 2; ++s) {
        REQUIRE(out(c, k, s) == 0);
        REQUIRE(out(c, k, 3 + s) ==
                Approx(before[c + f.nc * k + f.np * s] - state(c, k, s))
                    .epsilon(1e-10)
                    .margin(1e-25));
        REQUIRE(out(c, k, 3 + s) > 0);
      }
  REQUIRE(out(0, 0, 2) == 0);
  REQUIRE(out(0, 0, 5) == 0);
  const auto after = f.data;
  advanceModalWithCoating(*f.package, f.env(), state, 0, out);
  REQUIRE(f.data == after);
  REQUIRE(
      std::all_of(gain.begin(), gain.end(), [](double x) { return x == 0; }));
  auto config = chemistryConfig();
  config.coagulation = false;
  config.transfer_routes = {{"fresh", "aged", 2e-8}};
  ChemistryFixture renamed(config);
  const auto initial = renamed.data;
  auto view = renamed.view();
  renamed.package->executeMicrophysics(renamed.env(), view, 1);
  REQUIRE(renamed.data == initial);
  renameModalModes(*renamed.package, renamed.env(), view);
  REQUIRE(renamed.data != initial);
  for (int point = 0; point < f.np; ++point)
    for (int field = 0; field < 3; ++field)
      REQUIRE(renamed.data[point + f.np * field] +
                  renamed.data[point + f.np * (3 + field)] ==
              Approx(initial[point + f.np * field] +
                     initial[point + f.np * (3 + field)])
                  .epsilon(1e-12));
  config.transfer_owner = ModalTransferOwner::Microphysics;
  ChemistryFixture automatic(config);
  auto av = automatic.view();
  automatic.package->executeMicrophysics(automatic.env(), av, 1);
  REQUIRE(automatic.data == renamed.data);
  REQUIRE_THROWS(renameModalModes(*automatic.package, automatic.env(), av));
  REQUIRE_THROWS(
      advanceModalWithCoating(*automatic.package, automatic.env(), av, 1, out));
}

TEST_CASE(
    "native MAM4 coagulation exports positive primary-carbon coating gains") {
  if (!packageAvailable("mam4xx"))
    return;
  ChemistryFixture f(mam4xxReferenceModalConfig());
  const auto before = f.data;
  auto state = f.view();
  std::vector<double> gain(f.data.size());
  View3D<double> out(gain.data(), f.nc, f.nl, f.nf);
  advanceModalWithCoating(*f.package, f.env(), state, 60, out);
  REQUIRE(out(0, 0, 3 * 8 + 1) >
          0); // sulfate transferred from Aitken into primary carbon.
  for (int point = 0; point < f.np; ++point)
    for (int s = 0; s < 7; ++s) {
      const double accumulation =
          f.data[point + f.np * s] - before[point + f.np * s];
      REQUIRE(gain[point + f.np * s] ==
              Approx(accumulation).epsilon(1e-12).margin(1e-25));
    }
  for (int m = 0; m < 4; ++m)
    REQUIRE(out(0, 0, m * 8 + 7) == 0);
}

TEST_CASE("chemical coupling C API preserves layouts and rejects invalid "
          "dimensions") {
  ChemistryFixture f;
  const auto handle = static_cast<exaero_package_t>(f.package.get());
  std::array<char, 256> error{};
  int status = -1;
  const int s = 0;
  const double diff = 1e-5, mw = .098079, alpha = 1;
  std::vector<double> gas(f.np, 3e-12), eq(f.np * f.nm), ar(f.np * f.nf, -1),
      gr(f.np, -1), jac(f.np * (f.nf + 1) * (f.nf + 1), -1);
  exaero_modal_phase_exchange(
      handle, f.nc, f.nl, f.nf, f.nm, 1, f.nf + 1, &s, &diff, &mw, &alpha,
      f.t.data(), f.rh.data(), f.data.data(), gas.data(), eq.data(), ar.data(),
      gr.data(), jac.data(), nullptr, error.data(), &status);
  INFO(error.data());
  REQUIRE(status == 0);
  REQUIRE(ar[0] > 0);
  REQUIRE(gr[0] < 0);
  REQUIRE(ar[0] + ar[f.np * 3] + gr[0] == Approx(0).margin(1e-25));
  const auto saved_ar = ar, saved_gr = gr, saved_jac = jac;
  const int bad = -1;
  exaero_modal_phase_exchange(
      handle, f.nc, f.nl, f.nf, f.nm, 1, f.nf + 1, &bad, &diff, &mw, &alpha,
      f.t.data(), f.rh.data(), f.data.data(), gas.data(), eq.data(), ar.data(),
      gr.data(), jac.data(), nullptr, error.data(), &status);
  REQUIRE(status == 1);
  REQUIRE(ar == saved_ar);
  REQUIRE(gr == saved_gr);
  REQUIRE(jac == saved_jac);
  std::vector<double> coag(f.data.size()), rates(f.data.size(), -1),
      aging_jac(f.np * f.nf * (f.nf + 1), -1);
  exaero_modal_aging_terms(handle, f.nc, f.nl, f.nf, f.nf + 1, f.t.data(),
                           f.rh.data(), f.data.data(), ar.data(), coag.data(),
                           rates.data(), aging_jac.data(), jac.data(), f.nf + 1,
                           nullptr, 0, error.data(), &status);
  INFO(error.data());
  REQUIRE(status == 0);
  REQUIRE(rates[0] < 0);
  REQUIRE(rates[0] + rates[f.np * 3] == 0);
  const auto saved = rates;
  exaero_modal_aging_terms(
      handle, f.nc, f.nl, f.nf, std::numeric_limits<int>::max(), f.t.data(),
      f.rh.data(), f.data.data(), ar.data(), coag.data(), rates.data(),
      aging_jac.data(), nullptr, 0, nullptr, 0, error.data(), &status);
  REQUIRE(status == 1);
  REQUIRE(rates == saved);
  std::vector<double> gain(f.data.size(), -1);
  exaero_modal_microphysics_with_coating(handle, f.nc, f.nl, f.nf, f.t.data(),
                                         f.p.data(), f.rho.data(), f.rh.data(),
                                         f.dz.data(), f.data.data(), 0,
                                         gain.data(), error.data(), &status);
  REQUIRE(status == 0);
  REQUIRE(
      std::all_of(gain.begin(), gain.end(), [](double v) { return v == 0; }));
  const auto before = f.data;
  exaero_modal_rename_modes(handle, f.nc, f.nl, f.nf, f.t.data(), f.p.data(),
                            f.rho.data(), f.rh.data(), f.dz.data(),
                            f.data.data(), error.data(), &status);
  REQUIRE(status == 0);
  REQUIRE(f.data == before);
}

TEST_CASE("combined MKPP mapping preserves equivalent-material gas units and "
          "Jacobian scaling") {
  auto config = chemistryConfig();
  config.species[0].molecular_weight = .09606;
  ChemistryFixture f(config);
  std::vector<std::string> names;
  std::vector<ModalMechanismBinding> bindings;
  const auto fields = f.package->stateFields();
  for (std::size_t i = 0; i < fields.size(); ++i) {
    names.push_back("CHEM_" + std::to_string(i));
    bindings.push_back({fields[i].name, names.back()});
  }
  names.push_back("H2SO4_g");
  std::array<ModalVaporMechanismBinding, 1> gas{{{"acid", "H2SO4_g", 1}}};
  const auto map = modalChemicalMechanismMap(*f.package, names, bindings, gas);
  REQUIRE(map.size() == static_cast<std::size_t>(f.nf + 1));
  REQUIRE(map.back().local_index == f.nf);
  REQUIRE(map.back().to_mechanism == Approx(6.02214076e17 / .09606));
  REQUIRE(map[0].to_mechanism == map.back().to_mechanism);
  REQUIRE(map[2].to_mechanism == 1e-6);
  REQUIRE(modalChemicalJacobianScale(map[0], map[2]) ==
          Approx(map[0].to_mechanism / 1e-6));
  gas[0].mechanism_variable = names[0];
  REQUIRE_THROWS(modalChemicalMechanismMap(*f.package, names, bindings, gas));
}

TEST_CASE("host chemical coordinates enter rectangular exchange and aging "
          "Jacobians") {
  ChemistryFixture f;
  const int nq = f.nf + 1, ni = nq + 1;
  std::array<ModalCondensable, 1> binding{{{"acid", {1e-5, .098079, 1}}}};
  std::vector<double> gas(f.np, 4e-12), eq(f.np * f.nm, 1e-12),
      de(f.np * f.nm * ni), ar(f.np * f.nf), gr(f.np), jac(f.np * nq * ni),
      coag(f.np * f.nf), age(f.np * f.nf), aj(f.np * f.nf * ni);
  View4D<double> ej(de.data(), f.nc, f.nl, f.nm, ni);
  for (int m = 0; m < f.nm; ++m)
    for (int k = 0; k < f.nl; ++k)
      for (int c = 0; c < f.nc; ++c)
        ej(c, k, m, ni - 1) = .2;
  View3D<double> a(ar.data(), f.nc, f.nl, f.nf), g(gr.data(), f.nc, f.nl, 1),
      aging(age.data(), f.nc, f.nl, f.nf);
  View4D<double> j(jac.data(), f.nc, f.nl, nq, ni),
      agingj(aj.data(), f.nc, f.nl, f.nf, ni);
  const View4D<const double> eqj(de.data(), f.nc, f.nl, f.nm, ni);
  modalPhaseExchange(*f.package, f.env(), f.input(), binding,
                     View3D<const double>(gas.data(), f.nc, f.nl, 1),
                     View4D<const double>(eq.data(), f.nc, f.nl, f.nm, 1), a, g,
                     &j, &eqj);
  const View4D<const double> cj(jac.data(), f.nc, f.nl, nq, ni);
  modalAgingTerms(*f.package, f.env(), f.input(),
                  View3D<const double>(ar.data(), f.nc, f.nl, f.nf),
                  View3D<const double>(coag.data(), f.nc, f.nl, f.nf), aging,
                  &agingj, &cj);
  const double expected = j(0, 0, 0, ni - 1),
               expected_age = agingj(0, 0, 0, ni - 1), step = 1e-16;
  const auto evaluate_extra = [&](double delta) {
    eq[0] = 1e-12 + .2 * delta;
    modalPhaseExchange(*f.package, f.env(), f.input(), binding,
                       View3D<const double>(gas.data(), f.nc, f.nl, 1),
                       View4D<const double>(eq.data(), f.nc, f.nl, f.nm, 1), a,
                       g);
    modalAgingTerms(*f.package, f.env(), f.input(),
                    View3D<const double>(ar.data(), f.nc, f.nl, f.nf),
                    View3D<const double>(coag.data(), f.nc, f.nl, f.nf), aging);
    return std::array<double, 2>{ar[0], age[0]};
  };
  const auto plus = evaluate_extra(step), minus = evaluate_extra(-step);
  REQUIRE(expected == Approx((plus[0] - minus[0]) / (2 * step)).epsilon(1e-6));
  REQUIRE(expected_age ==
          Approx((plus[1] - minus[1]) / (2 * step)).epsilon(1e-6));
}

TEST_CASE("an implicit host step uses the coupled Jacobian without losing "
          "material mass") {
  auto config = chemistryConfig();
  config.species.resize(1);
  config.modes.resize(1);
  config.coagulation = false;
  config.coagulation_routes.clear();
  config.aging_routes.clear();
  auto package = createPackage("modal");
  initializeModal(*package, config);
  const double t = 280, p = 90000, rho = 1.1, rh = .5, dz = 100, dt = 1e4;
  const EnvironmentalStateView env{
      View2D<const double>(&t, 1, 1), View2D<const double>(&p, 1, 1),
      View2D<const double>(&rho, 1, 1), View2D<const double>(&rh, 1, 1),
      View2D<const double>(&dz, 1, 1)};
  const double initial_aerosol = 1e-9, initial_gas = 1e-9,
               total = initial_aerosol + initial_gas;
  const std::array<double, 1> mass{{initial_aerosol}};
  const double number = modalCoupling(*package).emittedNumber(0, mass);
  std::array<double, 2> state{{initial_aerosol, number}}, ar{};
  std::array<double, 9> jac{};
  double gas = initial_gas, gr = 0, eq = 0;
  const std::array<ModalCondensable, 1> binding{{{"acid", {1e-5, .098079, 1}}}};
  View3D<double> aout(ar.data(), 1, 1, 2), gout(&gr, 1, 1, 1);
  View4D<double> jout(jac.data(), 1, 1, 3, 3);
  const auto evaluate_step = [&] {
    state[0] = total - gas;
    modalPhaseExchange(
        *package, env, View3D<const double>(state.data(), 1, 1, 2), binding,
        View3D<const double>(&gas, 1, 1, 1),
        View4D<const double>(&eq, 1, 1, 1, 1), aout, gout, &jout);
  };
  for (int iteration = 0; iteration < 30; ++iteration) {
    evaluate_step();
    const double residual = gas - initial_gas - dt * gr;
    if (std::abs(residual) < 1e-23)
      break;
    const double derivative = 1 - dt * (jout(0, 0, 2, 2) - jout(0, 0, 2, 0));
    gas = std::clamp(gas - residual / derivative, 0., initial_gas);
  }
  evaluate_step();
  REQUIRE(gas > 0);
  REQUIRE(gas < initial_gas);
  REQUIRE(state[0] > initial_aerosol);
  REQUIRE(state[1] == number);
  REQUIRE(state[0] + gas == Approx(total).epsilon(1e-14));
  REQUIRE(gas - initial_gas == Approx(dt * gr).epsilon(1e-12));
  REQUIRE(state[0] - initial_aerosol == Approx(dt * ar[0]).epsilon(1e-12));
}
