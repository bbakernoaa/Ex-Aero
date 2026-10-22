#pragma once
/// @file Mam4xxPackage.hpp
/// @brief Experimental modal package; this translation-unit boundary excludes
/// Kokkos.
#include <exaero/IAerosolPackage.hpp>
#include <memory>

namespace exaero {
/// @brief MAM4xx clear-air coagulation, water uptake, and modal diagnostics.
/// Host views use layout_left. Component-table optics are opt-in.
class Mam4xxPackage final : public IAerosolPackage {
  bool initialized_ =
      false; ///< Successful configuration required before compute.
  std::unique_ptr<IAerosolPackage> optics_; ///< Optional shared modal optics.
  bool include_cloudborne_ = false; ///< Treat carried particles as aerosols.
  /// @brief Validate host state and environmental shapes, pointers and values.
  void validate(const EnvironmentalStateView &, const View3D<const double> &,
                bool full_environment = true) const;

public:
  /// @brief Modal representation with independent number fields.
  AerosolRepresentation representation() const override {
    return AerosolRepresentation::Modal;
  }
  /// @brief Supported surfaces; optics available after explicit configuration.
  unsigned capabilities() const override {
    return Microphysics | Diagnostics | CCN | Emissions |
           (optics_ ? Optics : 0u);
  }
  /// @brief All 64 state slots in public mam4xx_layout order.
  std::vector<StateField> stateFields() const override;
  /// @brief Configure native physics and optional component-table optics.
  void initialize(const std::string &config_yaml) override;
  /// @brief Advance interstitial coagulation; cloudborne state is conserved
  /// unchanged.
  void executeMicrophysics(const EnvironmentalStateView &, View3D<double> &,
                           double dt) override;
  /// @brief Dry mass, number, wet sizing, lognormal size cuts and column
  /// integrals.
  void computeDerivedDiagnostics(const EnvironmentalStateView &,
                                 const View3D<const double> &,
                                 View3D<double> &) override;
  /// @brief Integrate GEOSmie component tables over prognostic MAM4 modes.
  void computeOptics(const EnvironmentalStateView &,
                     const View3D<const double> &, const View1D<const double> &,
                     View4D<double> &) override;
  /// @brief Evaluate modal kappa-Kohler spectra at arbitrary nonnegative
  /// supersaturations.
  void computeCCN(const EnvironmentalStateView &, const View3D<const double> &,
                  const View1D<const double> &, View4D<double> &) override;
  /// @brief Identity mapping of already speciated modal mass/number fluxes.
  void computeEmissions(const EnvironmentalStateView &,
                        const EmissionsInputView &, View3D<double> &) override;
  /// @brief Look up a packed field name; return -1 if absent.
  int getSpeciesIndex(const std::string &name) const override;
  /// @brief Packed field name; throws on invalid index.
  std::string getSpeciesName(int index) const override;
};
} // namespace exaero
