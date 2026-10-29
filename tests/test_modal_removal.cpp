#include <algorithm>
#include <array>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include "removal_tendency_cases.hpp"
#include <cassert>
#include <cmath>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemoval.hpp>
#include <exaero/ModalRemovalMath.hpp>
#include <exaero/ModalRemoval_C.h>
#include <stdexcept>

namespace {
using namespace exaero;

struct Package final : IAerosolPackage, IModalCoupling {
  ModalConfig config;
  Package() {
    config.species.resize(2);
    config.species[0].name = "A";
    config.species[0].density = 1000;
    config.species[1].name = "B";
    config.species[1].density = 2000;
    config.modes.resize(2);
    config.modes[0].name = "small";
    config.modes[1].name = "large";
    config.modes[0].geometric_radius = 1e-7;
    config.modes[1].geometric_radius = 2e-7;
    config.modes[0].geometric_std_dev = config.modes[1].geometric_std_dev = 1.5;
    config.modes[0].dry_geometric_radius_min = 5e-8;
    config.modes[0].dry_geometric_radius_max = 2e-7;
    config.resuspension_routes = {{"A", "small", "large"}};
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
                        const EmissionsInputView &, View3D<double> &) override {
  }
  int getSpeciesIndex(const std::string &) const override { return -1; }
  std::string getSpeciesName(int) const override { return {}; }
  AerosolRepresentation representation() const override {
    return AerosolRepresentation::Modal;
  }
  void initializeModal(const ModalConfig &v) override { config = v; }
  const ModalConfig &modalConfig() const override { return config; }
  std::vector<ModalOpticalTableInfo> opticalTables() const override {
    return {};
  }
  void modalProperties(const EnvironmentalStateView &,
                       const View3D<const double> &, View4D<double> &out,
                       View4D<double> *) const override {
    for (std::size_t m = 0; m < config.modes.size(); ++m) {
      out(0, 0, m, modal_properties::DRY_RADIUS) =
          config.modes[m].geometric_radius;
      out(0, 0, m, modal_properties::WET_RADIUS) =
          1.5 * config.modes[m].geometric_radius;
    }
  }
  void modalUptake(const EnvironmentalStateView &, const View3D<const double> &,
                   std::span<const ModalVapor>, View4D<double> &,
                   View4D<double> *) const override {}
  ModalNucleationTerms nucleationTerms(double) const override { return {}; }
  double emittedNumber(std::size_t, std::span<const double>) const override {
    return 0;
  }
  void
  setModalOpticsProvider(std::shared_ptr<const IModalOpticsProvider>) override {
  }
  void setModalBandOpticsProvider(
      std::shared_ptr<const IModalBandOpticsProvider>) override {}
};

bool near(double x, double y, double tolerance = 1e-12) {
  return std::abs(x - y) <=
         tolerance * std::max({1e-30, std::abs(x), std::abs(y)});
}
} // namespace

int main() {
  Package p;
  // Distinct species losses in a mode, with independently supplied particle
  // loss.
  {
    std::array<double, 6> x{2e-9, 3e-9, 5e7, 0, 0, 0},
        rates{5e-4, 2e-4, 0, 0, 0, 0}, sink{};
    std::array<double, 2> nr{1e-4, 0}, nsink{};
    View3D<const double> xv(x.data(), 1, 1, 6), rv(rates.data(), 1, 1, 6),
        nv(nr.data(), 1, 1, 2);
    View3D<double> mv(sink.data(), 1, 1, 6), no(nsink.data(), 1, 1, 2);
    modalSpeciesRemovalRateTendencies(p, xv, rv, nv, mv, no);
    assert(near(sink[0], -1e-12) && near(sink[1], -6e-13) &&
           near(nsink[0], -5e3));
    assert(sink[2] == 0 && sink[5] == 0 && x[0] == 2e-9);
    for (double invalid : {1., -1., std::numeric_limits<double>::quiet_NaN()}) {
      rates[5] = invalid;
      sink.fill(77);
      nsink.fill(88);
      bool failed = false;
      try {
        modalSpeciesRemovalRateTendencies(p, xv, rv, nv, mv, no);
      } catch (const std::exception &) {
        failed = true;
      }
      assert(failed &&
             std::all_of(sink.begin(), sink.end(),
                         [](double v) { return v == 77; }) &&
             nsink[0] == 88);
    }
    rates[5] = 0;
    bool failed = false;
    View3D<const double> wrong(rates.data(), 1, 1, 2);
    try {
      modalSpeciesRemovalRateTendencies(p, xv, wrong, nv, mv, no);
    } catch (const std::exception &) {
      failed = true;
    }
    assert(failed && sink[0] == 77 && nsink[0] == 88);
    View3D<double> alias(x.data(), 1, 1, 6);
    failed = false;
    try {
      modalSpeciesRemovalRateTendencies(p, xv, rv, nv, alias, no);
    } catch (const std::exception &) {
      failed = true;
    }
    assert(failed && x[0] == 2e-9);
  }
  removal_test::tendencies(p, true);
  const double n = 1e8;
  const double mass = n * modalNumberVolume(4e-7, 1.5) * 1000;
  std::array<double, 6> raw{mass, 0, n, 0, 0, 0};
  View3D<double> state(raw.data(), 1, 1, 6);
  View3D<const double> input(raw.data(), 1, 1, 6);
  std::array<double, 6> report_data{};
  std::array<int, 2> status_data{};
  View4D<double> report(report_data.data(), 1, 1, 2, 3);
  View3D<int> status(status_data.data(), 1, 1, 2);
  inspectModalBounds(p, input, report, status);
  assert(status(0, 0, 0) == static_cast<int>(ModalBoundStatus::AboveRadius));
  assert(status(0, 0, 1) == static_cast<int>(ModalBoundStatus::Empty));
  assert(near(report(0, 0, 0, 0), 4e-7));
  assert(near(report(0, 0, 0, 1), 8 * n));
  assert(raw[2] == n);

  std::array<double, 2> delta_data{};
  View3D<double> delta(delta_data.data(), 1, 1, 2);
  enforceModalBounds(p, state, delta, status);
  assert(raw[0] == mass && raw[1] == 0 && raw[2] == 8 * n);
  assert(near(delta(0, 0, 0), 7 * n));
  assert(delta(0, 0, 1) == 0);
  enforceModalBounds(p, state, delta, status);
  assert(near(delta(0, 0, 0), 0));

  const auto before = raw;
  raw[2] = 0;
  bool failed = false;
  try {
    enforceModalBounds(p, state, delta, status);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && raw[0] == before[0] && raw[2] == 0);
  raw = before;

  std::array<double, 2> sizes_data{77, 77};
  std::array<double, 2> wet_data{77, 77};
  std::array<double, 4> diameters_data{};
  std::array<int, 2> populated_data{};
  View4D<double> diameters(diameters_data.data(), 1, 1, 2, 2);
  View3D<int> populated(populated_data.data(), 1, 1, 2);
  EnvironmentalStateView env{View2D<const double>(sizes_data.data(), 1, 1),
                             View2D<const double>(nullptr, 1, 1),
                             View2D<const double>(nullptr, 1, 1),
                             View2D<const double>(wet_data.data(), 1, 1),
                             View2D<const double>(nullptr, 1, 1)};
  modalScavengingSizes(p, env, input, diameters, populated);
  assert(populated(0, 0, 0) == 1 && populated(0, 0, 1) == 0);
  assert(near(diameters(0, 0, 0, 0), 2e-7));
  assert(near(diameters(0, 0, 0, 1), 3e-7));
  assert(diameters(0, 0, 1, 0) == 0);
  const auto original_diameters = diameters_data;
  const auto original_populated = populated_data;
  auto aliased_env = env;
  aliased_env.temperature = View2D<const double>(diameters_data.data(), 1, 1);
  failed = false;
  try {
    modalScavengingSizes(p, aliased_env, input, diameters, populated);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && diameters_data == original_diameters &&
         populated_data == original_populated);

  std::array<double, 6> rates_data{2, 3, 0, 5, 7, 0};
  std::array<double, 6> routed_data{};
  View3D<const double> rates(rates_data.data(), 1, 1, 6);
  View3D<double> routed(routed_data.data(), 1, 1, 6);
  redistributeModalResuspension(p, rates, routed);
  assert(routed_data == (std::array<double, 6>{0, 3, 0, 7, 7, 0}));
  rates_data[2] = 1;
  failed = false;
  try {
    redistributeModalResuspension(p, rates, routed);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && routed_data == (std::array<double, 6>{0, 3, 0, 7, 7, 0}));

  std::array<double, 6> budget_data{mass, 0, n, 0, 0, 0};
  std::array<double, 3> edges_data{0, 8e-7, INFINITY};
  std::array<double, 2> bin_data{}, tail_data{};
  View3D<const double> budget(budget_data.data(), 1, 1, 6);
  View1D<const double> edges(edges_data.data(), 3);
  View3D<double> bins(bin_data.data(), 1, 1, 2),
      tails(tail_data.data(), 1, 1, 2);
  rebinModalDeposition(p, budget, 0, edges, bins, tails);
  assert(near(bin_data[0] + bin_data[1] + tail_data[0] + tail_data[1], mass));
  assert(near(bin_data[0] / mass, modalMassFractionBelow(8e-7, 4e-7, 1.5)));
  assert(tail_data[0] == 0 && tail_data[1] == 0);
  edges_data = {4e-7, 8e-7, 16e-7};
  rebinModalDeposition(p, budget, 0, edges, bins, tails);
  assert(tail_data[0] > 0 && tail_data[1] > 0);
  assert(near(bin_data[0] + bin_data[1] + tail_data[0] + tail_data[1], mass));
  const auto valid_bins = bin_data;
  const auto valid_tails = tail_data;
  edges_data[1] = edges_data[0];
  failed = false;
  try {
    rebinModalDeposition(p, budget, 0, edges, bins, tails);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && bin_data == valid_bins && tail_data == valid_tails);
  edges_data = {8e-6, 1e-5, 1.2e-5};
  rebinModalDeposition(p, budget, 0, edges, bins, tails);
  assert(bin_data[0] > 0 && bin_data[1] > 0);
  assert(near(bin_data[0] / mass,
              modalMassFractionAbove(edges_data[0], 4e-7, 1.5) -
                  modalMassFractionAbove(edges_data[1], 4e-7, 1.5)));
  edges_data = {0, 8e-7, INFINITY};

  char error[256]{};
  int flag = -1;
  exaero_modal_inspect_bounds(static_cast<exaero_package_t>(&p), 1, 1, 6, 2,
                              raw.data(), report_data.data(),
                              status_data.data(), error, &flag);
  assert(flag == 0 && error[0] == '\0');
  assert(status_data[0] == EXAERO_MODAL_BOUNDS_WITHIN);
  exaero_modal_enforce_bounds(static_cast<exaero_package_t>(&p), 1, 1, 6, 2,
                              raw.data(), delta_data.data(), status_data.data(),
                              error, &flag);
  assert(flag == 0 && delta_data[0] == 0);
  exaero_modal_scavenging_sizes(static_cast<exaero_package_t>(&p), 1, 1, 6, 2,
                                sizes_data.data(), wet_data.data(), raw.data(),
                                diameters_data.data(), populated_data.data(),
                                error, &flag);
  assert(flag == 0 && populated_data[0] == 1);
  rates_data[2] = 0;
  exaero_modal_resuspension(static_cast<exaero_package_t>(&p), 1, 1, 6,
                            rates_data.data(), routed_data.data(), error,
                            &flag);
  assert(flag == 0 && routed_data[3] == 7);
  exaero_modal_rebin_deposition(
      static_cast<exaero_package_t>(&p), 1, 1, 6, 0, 2, budget_data.data(),
      edges_data.data(), bin_data.data(), tail_data.data(), error, &flag);
  assert(flag == 0 && near(bin_data[0] + bin_data[1], mass));
  exaero_modal_inspect_bounds(static_cast<exaero_package_t>(&p), 1, 1, 5, 2,
                              raw.data(), report_data.data(),
                              status_data.data(), error, &flag);
  assert(flag == 1 && error[0] != '\0');
}
