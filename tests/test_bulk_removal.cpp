#include <exaero/BulkRemoval.hpp>
#include <exaero/BulkRemoval_C.h>
#include <exaero/ModalRemovalMath.hpp>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <array>
#include <cassert>
#include <cmath>

namespace {
using namespace exaero;
struct Package final : IAerosolPackage, IBulkRemovalMetadata {
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
  std::vector<StateField> stateFields() const override {
    return {{"A", StateQuantity::Mass, -1, false},
            {"B", StateQuantity::Mass, -1, false}};
  }
  BulkSpeciesSize bulkSpeciesSize(std::size_t s) const override {
    return s == 0 ? BulkSpeciesSize{8e-7, 1000, 8e-7, 1.5, 2e-7, 8e-7, 4e-7, .5}
                  : BulkSpeciesSize{2e-6, 2000, 2e-6, 1, 0, 0, 0, 0};
  }
  void computeWetSizing(const View1D<const double> &rh, int s, double radius,
                        double density, View1D<double> &wet,
                        View1D<double> &wet_density,
                        View1D<int> *status) override {
    (void)rh;
    (void)s;
    (void)radius;
    (void)density;
    (void)wet;
    (void)wet_density;
    (void)status;
    throw std::logic_error(
        "bulk utility must not query radiation-backed sizing");
  }
};
bool near(double x, double y) {
  return std::abs(x - y) <= 1e-12 * std::max({1e-30, std::abs(x), std::abs(y)});
}
} // namespace

int main() {
  Package p;
  std::array<double, 4> state_data{2, 0, 4, 5};
  std::array<double, 2> rh_data{.5, .2};
  std::array<double, 8> diameter_data{};
  std::array<int, 4> size_status_data{};
  View3D<const double> state(state_data.data(), 2, 1, 2);
  View2D<const double> rh(rh_data.data(), 2, 1);
  View4D<double> diameters(diameter_data.data(), 2, 1, 2, 2);
  View3D<int> status(size_status_data.data(), 2, 1, 2);
  std::array<double, 4> number_data{};
  View3D<double> number(number_data.data(), 2, 1, 2);
  bulkFixedNumber(p, state, number);
  const double particle_mass = 1000 * modalNumberVolume(4e-7, 1.5);
  assert(near(number(0, 0, 0), 2 / particle_mass));
  assert(number(1, 0, 0) == 0);
  const auto original_number = number_data;
  state_data[0] = -1;
  bool failed = false;
  try {
    bulkFixedNumber(p, state, number);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && number_data == original_number);
  state_data[0] = 2;
  std::array<double, 6> radii_data{};
  std::array<int, 6> radii_status_data{};
  View2D<double> radii(radii_data.data(), 2, 3);
  View2D<int> radii_status(radii_status_data.data(), 2, 3);
  bulkReferenceRadii(p, radii, radii_status);
  assert(near(radii(0, 0), 2e-7) && near(radii(0, 1), 8e-7));
  assert(near(radii(0, 2), 4e-7));
  assert(radii_status(1, 0) == static_cast<int>(BulkRadiusStatus::Unavailable));
  assert(radii_status(1, 2) == static_cast<int>(BulkRadiusStatus::Derived));
  assert(near(radii(1, 2), 1e-6));
  bulkScavengingSizes(p, rh, state, diameters, status);
  assert(near(diameters(0, 0, 0, 0), 8e-7));
  assert(near(diameters(0, 0, 0, 1), 8e-7 * std::cbrt(1.5)));
  assert(status(1, 0, 0) == BULK_SIZE_EMPTY && diameters(1, 0, 0, 1) == 0);
  assert(status(0, 0, 1) == BULK_SIZE_CONFIGURED);
  assert(near(diameters(0, 0, 1, 1), 2e-6));
  state_data[1] = 2;
  rh_data[1] = 1;
  bulkScavengingSizes(p, rh, state, diameters, status);
  assert(near(diameters(1, 0, 0, 1), 8e-7 * std::cbrt(1 + .5 * .99 / .01)));
  state_data[1] = 0;
  rh_data[1] = .2;

  const std::array<BulkResuspensionRoute, 1> routes{{{0, 1}}};
  std::array<double, 4> rates_data{2, 3, 4, 5}, routed_data{};
  View3D<const double> rates(rates_data.data(), 2, 1, 2);
  View3D<double> routed(routed_data.data(), 2, 1, 2);
  redistributeBulkResuspension(p, rates, routes, routed);
  assert(routed_data == (std::array<double, 4>{0, 0, 6, 8}));
  rates_data[0] = -1;
  failed = false;
  try {
    redistributeBulkResuspension(p, rates, routes, routed);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && routed_data == (std::array<double, 4>{0, 0, 6, 8}));

  std::array<double, 4> budget_data{2, 3, 4, 5};
  std::array<double, 3> edges_data{0, 8e-7, INFINITY};
  std::array<double, 4> bins_data{}, tails_data{};
  View3D<const double> budget(budget_data.data(), 2, 1, 2);
  View1D<const double> edges(edges_data.data(), 3);
  View3D<double> bins(bins_data.data(), 2, 1, 2),
      tails(tails_data.data(), 2, 1, 2);
  rebinBulkDeposition(p, budget, 0, edges, bins, tails);
  const double half = modalMassFractionBelow(8e-7, 4e-7, 1.5);
  assert(near(bins(0, 0, 0), 2 * half));
  assert(near(bins(1, 0, 0), 3 * half));
  assert(near(bins(0, 0, 0) + bins(0, 0, 1), 2));
  assert(tails(0, 0, 0) == 0 && tails(0, 0, 1) == 0);
  edges_data = {1e-6, 2e-6, 3e-6};
  rebinBulkDeposition(p, budget, 1, edges, bins, tails);
  assert(bins(0, 0, 0) == 0 && bins(0, 0, 1) == 4);
  assert(tails(0, 0, 0) == 0 && tails(0, 0, 1) == 0);
  const auto original_bins = bins_data;
  edges_data[1] = edges_data[0];
  failed = false;
  try {
    rebinBulkDeposition(p, budget, 1, edges, bins, tails);
  } catch (const std::invalid_argument &) {
    failed = true;
  }
  assert(failed && bins_data == original_bins);
  edges_data[1] = 2e-6;
  char error[256]{};
  int flag = -1;
  exaero_bulk_fixed_number(static_cast<exaero_package_t>(&p), 2, 1, 2,
                           state_data.data(), number_data.data(), error, &flag);
  assert(flag == 0 && near(number_data[0], 2 / particle_mass));
  exaero_bulk_reference_radii(static_cast<exaero_package_t>(&p), 2,
                              radii_data.data(), radii_status_data.data(),
                              error, &flag);
  assert(flag == 0 && near(radii_data[0], 2e-7));
  exaero_bulk_scavenging_sizes(static_cast<exaero_package_t>(&p), 2, 1, 2,
                               rh_data.data(), state_data.data(),
                               diameter_data.data(), size_status_data.data(),
                               error, &flag);
  assert(flag == 0 && error[0] == '\0');
  const int sources[1]{0}, destinations[1]{1};
  rates_data[0] = 2;
  exaero_bulk_resuspension(static_cast<exaero_package_t>(&p), 2, 1, 2, 1,
                           sources, destinations, rates_data.data(),
                           routed_data.data(), error, &flag);
  assert(flag == 0 && routed_data[2] == 6);
  exaero_bulk_rebin_deposition(
      static_cast<exaero_package_t>(&p), 2, 1, 2, 1, 2, budget_data.data(),
      edges_data.data(), bins_data.data(), tails_data.data(), error, &flag);
  assert(flag == 0 && bins_data[2] == 4);
  exaero_bulk_scavenging_sizes(static_cast<exaero_package_t>(&p), 2, 1, 3,
                               rh_data.data(), state_data.data(),
                               diameter_data.data(), size_status_data.data(),
                               error, &flag);
  assert(flag == 1 && error[0] != '\0');
}
