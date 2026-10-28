/// @file test_mam4xx.cpp
/// @brief Native-kernel parity, conservation, units, layout and C ABI tests.
#define CATCH_CONFIG_RUNNER
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <exaero/AerosolIndices.hpp>
#include <exaero/IAerosolPackage_C.h>
#include <exaero/Mam4xxConfig.hpp>
#include <exaero/OpticsBandCoupling.hpp>
#include <limits>
#include <mam4xx/coagulation.hpp>
#include <mam4xx/conversions.hpp>
#include <memory>
#include <numeric>
#include <vector>

namespace {
namespace layout = exaero::mam4xx_layout;
constexpr int nc = 2, nl = 3, points = nc * nl;
/// @brief Host-only column-major test fixture with distinct cells and levels.
struct Fixture {
  std::unique_ptr<void, decltype(&exaero_free_package)> pkg{
      nullptr, exaero_free_package};
  std::array<double, points> t{}, p{}, rho{}, rh{}, dz{};
  std::vector<double> state = std::vector<double>(points * layout::fields);
  std::array<char, 256> error{};
  int status = -1;
  /// @brief Initialize a modal package and two aerosol phases at every grid
  /// point.
  Fixture() {
    pkg.reset(exaero_create_package("mam4xx", error.data(), &status));
    REQUIRE(pkg);
    REQUIRE(status == 0);
    exaero_initialize_package(pkg.get(), "backend: mam4xx", error.data(),
                              &status);
    REQUIRE(status == 0);
    for (int point = 0; point < points; ++point) {
      t[point] = 280 + point;
      p[point] = 90000 - 1000 * point;
      rho[point] = 1 + 0.1 * point;
      rh[point] = 0;
      dz[point] = 100 + 10 * point;
      for (int phase = 0; phase < 2; ++phase)
        for (int m = 0; m < 4; ++m) {
          const double number = 1e8 * (m + 1) * (point + 1) * (phase ? 0.1 : 1);
          const auto mode = mam4::modes(m);
          state[point + points * layout::number(m, phase)] = number;
          state[point +
                points * layout::mass(m, static_cast<int>(mam4::AeroId::SO4),
                                      phase)] =
              number * mam4::defaults::mam4_density_so4 *
              mam4::conversions::mean_particle_volume_from_diameter(
                  mode.nom_diameter, mode.mean_std_dev);
        }
    }
  }
  /// @brief Advance through the public C ABI.
  void advance(double dt) {
    exaero_execute_microphysics(pkg.get(), nc, nl, layout::fields, t.data(),
                                p.data(), rho.data(), rh.data(), dz.data(),
                                state.data(), dt, error.data(), &status);
  }
  /// @brief Diagnostics through the same generic C wrapper used by GOCART.
  std::vector<double> diagnostics() {
    std::vector<double> out(
        points * exaero::diagnostic_indices::NUM_DIAGNOSTICS, -1);
    exaero_compute_diagnostics(pkg.get(), nc, nl, layout::fields, t.data(),
                               p.data(), rho.data(), rh.data(), dz.data(),
                               state.data(), out.data(), error.data(), &status);
    REQUIRE(status == 0);
    return out;
  }
};
} // namespace

TEST_CASE("modal metadata preserves independent mass, number and phase") {
  Fixture f;
  REQUIRE(exaero_get_representation(f.pkg.get()) == 1);
  REQUIRE((exaero_get_capabilities(f.pkg.get()) & 4) == 0);
  for (int phase = 0; phase < 2; ++phase)
    for (int m = 0; m < 4; ++m) {
      int quantity = -1, mode = -1, cloud = -1;
      REQUIRE(exaero_get_state_field_info(f.pkg.get(), layout::number(m, phase),
                                          &quantity, &mode, &cloud) == 0);
      REQUIRE(quantity == 1);
      REQUIRE(mode == m);
      REQUIRE(cloud == phase);
      REQUIRE(exaero_get_state_field_info(f.pkg.get(),
                                          layout::mass(m, 1, phase), &quantity,
                                          &mode, &cloud) == 0);
      REQUIRE(quantity == 0);
    }
  REQUIRE(exaero_get_species_index(f.pkg.get(), "num_aitken_i") ==
          layout::number(1));
  REQUIRE(exaero_get_species_index(f.pkg.get(), "SO4_coarse_c") ==
          layout::mass(2, 1, true));
  REQUIRE(exaero_get_species_index(f.pkg.get(), "invalid") == -1);
}

TEST_CASE("MAM4xx delegates native band metadata and preserves unsupported output") {
  Fixture f;
  exaero_initialize_package(f.pkg.get(),
                            "optics: {treatment: component_tables}",
                            f.error.data(), &f.status);
  REQUIRE(f.status==0);
  auto *base=static_cast<exaero::IAerosolPackage *>(f.pkg.get());
  auto *bands=dynamic_cast<exaero::IOpticsBandCoupling *>(base);
  REQUIRE(bands);
  REQUIRE((base->capabilities() & exaero::OpticsBands)!=0);
  REQUIRE(bands->radiationBandInfo(exaero::RadiationBandSet::RRTMG_30).size()==30);
  const auto availability=bands->bandAvailability(exaero::RadiationBandSet::RRTMG_30);
  REQUIRE(availability.size()==static_cast<std::size_t>(layout::modes*layout::species));
  for (const auto &component:availability)
    for (const auto status:component.status)
      REQUIRE(status==exaero::RadiationBandStatus::NotInSource);
  exaero::EnvironmentalStateView env{
    exaero::View2D<const double>(f.t.data(),nc,nl),
    exaero::View2D<const double>(f.p.data(),nc,nl),
    exaero::View2D<const double>(f.rho.data(),nc,nl),
    exaero::View2D<const double>(f.rh.data(),nc,nl),
    exaero::View2D<const double>(f.dz.data(),nc,nl)};
  exaero::View3D<const double> state(f.state.data(),nc,nl,layout::fields);
  std::vector<double> output(points*exaero::optical_indices::NUM_OPTICS,91.);
  exaero::View4D<double> out(output.data(),nc,nl,1,
                             exaero::optical_indices::NUM_OPTICS);
  const int selected[]={2};
  REQUIRE_THROWS(bands->computeBandOptics(env,state,
      exaero::RadiationBandSet::RRTMG_30,selected,out));
  for (double value:output) REQUIRE(value==91.);
}

TEST_CASE(
    "coagulation matches the native MAM4xx kernel on every cell and level") {
  Fixture f;
  const auto before = f.state;
  constexpr double dt = 60;
  f.advance(dt);
  INFO(f.error.data());
  REQUIRE(f.status == 0);
  for (int point = 0; point < points; ++point) {
    const double air = f.p[point] / (mam4::Constants::r_gas * f.t[point]);
    double number[4], mass[7][4]{}, diameter[4], density[4], transfer[7][1]{};
    double total_before = 0, total_after = 0, n_before = 0, n_after = 0;
    for (int m = 0; m < 4; ++m) {
      number[m] = before[point + points * layout::number(m)] / air;
      diameter[m] = mam4::modes(m).nom_diameter;
      density[m] = mam4::defaults::mam4_density_so4;
      for (int s = 0; s < 7; ++s)
        mass[s][m] = before[point + points * layout::mass(m, s)];
    }
    mam4::coagulation::mam_coag_1subarea(dt, f.t[point], f.p[point], air,
                                         diameter, density, number, mass,
                                         transfer);
    for (int m = 0; m < 4; ++m) {
      REQUIRE(f.state[point + points * layout::number(m)] ==
              Approx(number[m] * air).epsilon(1e-12));
      n_before += before[point + points * layout::number(m)];
      n_after += f.state[point + points * layout::number(m)];
      for (int s = 0; s < 7; ++s) {
        REQUIRE(f.state[point + points * layout::mass(m, s)] ==
                Approx(mass[s][m]).epsilon(1e-12));
        total_before += before[point + points * layout::mass(m, s)];
        total_after += f.state[point + points * layout::mass(m, s)];
      }
    }
    REQUIRE(total_after == Approx(total_before).epsilon(1e-13));
    REQUIRE(n_after < n_before);
  }
  REQUIRE(std::equal(before.begin() + points * layout::fields_per_phase,
                     before.end(),
                     f.state.begin() + points * layout::fields_per_phase));
}

TEST_CASE("modal diagnostics integrate both phases with layout_left indexing") {
  Fixture f;
  const auto d = f.diagnostics();
  using namespace exaero::diagnostic_indices;
  for (int point = 0; point < points; ++point) {
    double mass = 0, number = 0;
    for (int phase = 0; phase < 2; ++phase)
      for (int m = 0; m < 4; ++m) {
        mass += f.state[point + points * layout::mass(m, 1, phase)];
        number += f.state[point + points * layout::number(m, phase)];
      }
    REQUIRE(d[point + points * MASS_CONCENTRATION] ==
            Approx(mass).epsilon(1e-13));
    REQUIRE(d[point + points * NUMBER_CONCENTRATION] ==
            Approx(number).epsilon(1e-13));
    REQUIRE(d[point + points * AEROSOL_LIQUID_WATER] == 0);
    REQUIRE(d[point + points * PM10_CONCENTRATION] >=
            d[point + points * PM2_5_CONCENTRATION]);
  }
  for (int c = 0; c < nc; ++c) {
    double column = 0;
    for (int k = 0; k < nl; ++k)
      column += d[c + nc * k + points * MASS_CONCENTRATION] * f.dz[c + nc * k];
    REQUIRE(d[c + points * COLUMN_MASS] == Approx(column).epsilon(1e-13));
  }
  f.rh.fill(0.9);
  const auto wet = f.diagnostics();
  REQUIRE(wet[points * AEROSOL_LIQUID_WATER] > 0);
  REQUIRE(wet[points * SURFACE_AREA_DENSITY] >
          d[points * SURFACE_AREA_DENSITY]);
}

TEST_CASE("CCN spectra are bounded by prognostic number and monotonic") {
  Fixture f;
  const double axis[] = {0, 0.0002, 0.001, 0.01};
  std::vector<double> out(points * 4, -1);
  exaero_compute_ccn(f.pkg.get(), nc, nl, 4, layout::fields, axis, f.t.data(),
                     f.rh.data(), f.state.data(), out.data(), f.error.data(),
                     &f.status);
  INFO(f.error.data());
  REQUIRE(f.status == 0);
  const auto d = f.diagnostics();
  for (int point = 0; point < points; ++point) {
    REQUIRE(out[point] == 0);
    for (int s = 1; s < 4; ++s) {
      REQUIRE(out[point + points * s] >= out[point + points * (s - 1)]);
      REQUIRE(
          out[point + points * s] <=
          d[point + points * exaero::diagnostic_indices::NUMBER_CONCENTRATION]);
    }
  }
}

TEST_CASE("invalid modal inputs fail before mutating state and unsupported "
          "optics errors") {
  Fixture f;
  const auto before = f.state;
  f.advance(0);
  REQUIRE(f.status == 0);
  REQUIRE(f.state == before);
  f.advance(-1);
  REQUIRE(f.status == 1);
  REQUIRE(f.state == before);
  f.t[3] = std::numeric_limits<double>::quiet_NaN();
  f.advance(10);
  REQUIRE(f.status == 1);
  REQUIRE(f.state == before);
  f.t[3] = 283;
  exaero_initialize_package(f.pkg.get(), "species: []", f.error.data(),
                            &f.status);
  REQUIRE(f.status == 1);
  std::vector<double> optics(points * exaero::optical_indices::NUM_OPTICS, -1);
  const double wavelength = 550e-9;
  exaero_compute_optics(f.pkg.get(), nc, nl, 1, layout::fields, &wavelength,
                        f.t.data(), f.p.data(), f.rho.data(), f.rh.data(),
                        f.dz.data(), f.state.data(), optics.data(),
                        f.error.data(), &f.status);
  REQUIRE(f.status == 1);
  REQUIRE(std::all_of(optics.begin(), optics.end(),
                      [](double x) { return x == -1; }));
}

TEST_CASE("modal emissions preserve mass and number units and divide area flux "
          "by layer depth") {
  Fixture f;
  std::vector<double> flux(points * layout::fields), out(flux.size(), -1);
  for (int field = 0; field < layout::fields; ++field)
    for (int point = 0; point < points; ++point)
      flux[point + points * field] = (field + 1) * (point + 1) * 1e-12;
  exaero_compute_emissions(f.pkg.get(), nc, nl, layout::fields, layout::fields,
                           0, nullptr, flux.data(), out.data(), f.error.data(),
                           &f.status);
  REQUIRE(f.status == 0);
  REQUIRE(out == flux);
  exaero_compute_emissions(f.pkg.get(), nc, nl, layout::fields, layout::fields,
                           1, f.dz.data(), flux.data(), out.data(),
                           f.error.data(), &f.status);
  REQUIRE(f.status == 0);
  for (int field = 0; field < layout::fields; ++field)
    for (int point = 0; point < points; ++point)
      REQUIRE(
          out[point + points * field] ==
          Approx(flux[point + points * field] / f.dz[point]).epsilon(1e-13));
  const auto before = out;
  f.dz[2] = 0;
  exaero_compute_emissions(f.pkg.get(), nc, nl, layout::fields, layout::fields,
                           1, f.dz.data(), flux.data(), out.data(),
                           f.error.data(), &f.status);
  REQUIRE(f.status == 1);
  REQUIRE(out == before);
}

TEST_CASE("empty and inconsistent modes have explicit behavior") {
  Fixture f;
  std::fill(f.state.begin(), f.state.end(), 0);
  f.advance(60);
  REQUIRE(f.status == 0);
  const auto d = f.diagnostics();
  REQUIRE(std::all_of(d.begin(), d.end(), [](double x) { return x == 0; }));
  f.state[layout::number(0) * points] = 1;
  const auto before = f.state;
  f.advance(60);
  REQUIRE(f.status == 1);
  REQUIRE(f.state == before);
  exaero_execute_microphysics(f.pkg.get(), nc, nl, layout::fields - 1,
                              f.t.data(), f.p.data(), f.rho.data(), f.rh.data(),
                              f.dz.data(), f.state.data(), 60, f.error.data(),
                              &f.status);
  REQUIRE(f.status == 1);
  REQUIRE(f.state == before);
}

/// @brief One Kokkos lifetime for the full test session.
int main(int argc, char **argv) {
  const Kokkos::ScopeGuard guard(argc, argv);
  return Catch::Session().run(argc, argv);
}

TEST_CASE("MAM4xx component optics use native mass, number and phase layout") {
  Fixture f;
  const auto configure = [&](const char *yaml) {
    exaero_initialize_package(f.pkg.get(), yaml, f.error.data(), &f.status);
    INFO(f.error.data());
    REQUIRE(f.status == 0);
  };
  const auto compute = [&](double wavelength = 550e-9, bool success = true) {
    std::vector<double> out(points * exaero::optical_indices::NUM_OPTICS, -1);
    exaero_compute_optics(f.pkg.get(), nc, nl, 1, layout::fields, &wavelength,
                          f.t.data(), f.p.data(), f.rho.data(), f.rh.data(),
                          f.dz.data(), f.state.data(), out.data(),
                          f.error.data(), &f.status);
    INFO(f.error.data());
    REQUIRE(f.status == (success ? 0 : 1));
    if (!success)
      REQUIRE(std::all_of(out.begin(), out.end(),
                          [](double x) { return x == -1; }));
    return out;
  };
  configure("backend: mam4xx\noptics: {treatment: component_tables}");
  REQUIRE((exaero_get_capabilities(f.pkg.get()) & 4) != 0);
  const auto original = f.state;
  const auto baseline = compute();
  // Independent Simpson integration of the committed GEOSmie sulfate curves:
  // 200001 log-radius samples, native widths, geometric radii and fixture N.
  // This checks native species ordering, density, radius/diameter and units.
  REQUIRE(baseline[exaero::optical_indices::EXTINCTION_COEFF * points] ==
          Approx(0.004020626577863798).epsilon(1e-7));
  REQUIRE(baseline[exaero::optical_indices::SCATTERING_COEFF * points] ==
          Approx(0.004020626577845975).epsilon(1e-7));
  REQUIRE(baseline[exaero::optical_indices::LIDAR_BACKSCATTER * points] ==
          Approx(0.00011513250778826159).epsilon(1e-7));
  REQUIRE(baseline[exaero::optical_indices::ASYMMETRY_FACTOR * points] ==
          Approx(0.5430516938975144).epsilon(1e-7));
  for (int c = 0; c < nc; ++c) {
    double column = 0;
    for (int k = 0; k < nl; ++k)
      column += baseline[c + nc * k] * f.dz[c + nc * k];
    REQUIRE(baseline[c + points * exaero::optical_indices::EXTINCTION_AOT] ==
            Approx(column).epsilon(1e-11));
  }
  REQUIRE(f.state == original);
  for (double value : baseline)
    REQUIRE(std::isfinite(value));
  // Increasing both mass and number preserves size and doubles extinction.
  for (auto &value : f.state)
    value *= 2;
  const auto doubled = compute();
  for (int point = 0; point < points; ++point) {
    REQUIRE(baseline[point] > 0);
    REQUIRE(doubled[point] == Approx(2 * baseline[point]).epsilon(1e-11));
  }
  f.state = original;
  configure("optics: {treatment: component_tables, include_cloudborne: true}");
  const auto both = compute();
  for (int point = 0; point < points; ++point)
    REQUIRE(both[point] == Approx(1.1 * baseline[point]).epsilon(1e-11));
  configure("optics: {treatment: component_tables}");
  // A different cloudborne population is ignored by the default optical phase.
  for (int field = layout::fields_per_phase; field < layout::fields; ++field)
    for (int point = 0; point < points; ++point)
      f.state[point + points * field] *= 10;
  REQUIRE(compute() == baseline);
  compute(100e-9, false);
  // Invalid replacement configurations preserve the previously usable optics.
  for (const auto *yaml :
       {"optics: {treatment: internal_mixture}",
        "optics: {treatment: component_tables, unknown: true}",
        "optics: {treatment: component_tables, sources: {NO3: MODAL_NI}}",
        "optics: {treatment: component_tables, sources: {SO4: MODAL_SU}}",
        "optics: {treatment: component_tables, sources: {SOA: missing}}"}) {
    exaero_initialize_package(f.pkg.get(), yaml, f.error.data(), &f.status);
    REQUIRE(f.status == 1);
    REQUIRE(compute() == baseline);
  }
  configure("backend: mam4xx");
  REQUIRE((exaero_get_capabilities(f.pkg.get()) & 4) == 0);
  compute(550e-9, false);
}

TEST_CASE(
    "GEOSmie MAM4xx optics support all seven materials at dry and humid RH") {
  Fixture f;
  exaero_initialize_package(f.pkg.get(),
                            "optics: {treatment: component_tables}",
                            f.error.data(), &f.status);
  INFO(f.error.data());
  REQUIRE(f.status == 0);
  const auto species = mam4::default_aero_species();
  const double wavelength = 550e-9;
  for (int s = 0; s < layout::species; ++s) {
    std::fill(f.state.begin(), f.state.end(), 0);
    for (int m = 0; m < layout::modes; ++m)
      for (int point = 0; point < points; ++point) {
        const auto mode = mam4::modes(m);
        const double number = 1e8;
        f.state[point + points * layout::number(m)] = number;
        f.state[point + points * layout::mass(m, s)] =
            number * species(s).density *
            mam4::conversions::mean_particle_volume_from_diameter(
                mode.nom_diameter, mode.mean_std_dev);
      }
    for (const double rh : {0.0, 0.8, 0.98}) {
      f.rh.fill(rh);
      std::vector<double> out(points * exaero::optical_indices::NUM_OPTICS, -1);
      exaero_compute_optics(f.pkg.get(), nc, nl, 1, layout::fields, &wavelength,
                            f.t.data(), f.p.data(), f.rho.data(), f.rh.data(),
                            f.dz.data(), f.state.data(), out.data(),
                            f.error.data(), &f.status);
      INFO(s);
      INFO(rh);
      INFO(f.error.data());
      REQUIRE(f.status == 0);
      for (int point = 0; point < points; ++point)
        REQUIRE(out[point] > 0);
      for (double value : out)
        REQUIRE(std::isfinite(value));
    }
  }
}
