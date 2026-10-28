#pragma once
/// @file ModalPackage.hpp
/// @brief Validated host adapter for the generic modal kernels.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: createPackage("modal"). No file output.
#include <exaero/ModalCoupling.hpp>
#include <exaero/OpticsBandCoupling.hpp>
#include <modal/ModalSolver.hpp>

namespace exaero {
/// @brief Two-moment modal backend with optional table/provider optics.
class ModalPackage final : public IAerosolPackage, public IModalCoupling,
                           public IOpticsBandCoupling {
  ModalConfig config_;
  modal::Parameters parameters_;
  std::vector<modal::Table> tables_;
  std::vector<StateField> fields_;
  bool initialized_ = false;
  mutable std::shared_ptr<modal::Solver> solver_;
  std::shared_ptr<const IModalOpticsProvider> optics_provider_;
  std::shared_ptr<const IModalBandOpticsProvider> band_optics_provider_;
  /// @brief Validate host arrays and material/number consistency atomically.
  void validate(const EnvironmentalStateView &, const View3D<const double> &,
                bool full_environment = true) const;
  /// @brief Lazily upload immutable parameters after Kokkos initialization.
  const modal::Solver &solver() const;

public:
  /// @brief Modal representation with independently prognosed numbers.
  AerosolRepresentation representation() const override {
    return AerosolRepresentation::Modal;
  }
  /// @brief State-axis metadata in configuration mode/species order.
  std::vector<StateField> stateFields() const override { return fields_; }
  /// @brief Base surfaces plus optics when all modes have their dependencies.
  unsigned capabilities() const override;
  /// @brief Parse YAML and use the same transactional validation as C++ config.
  void initialize(const std::string &) override;
  void initializeModal(const ModalConfig &) override;
  const ModalConfig &modalConfig() const override { return config_; }
  std::vector<ModalOpticalTableInfo> opticalTables() const override;
  void executeMicrophysics(const EnvironmentalStateView &, View3D<double> &,
                           double) override;
  /// @brief Validate a read-only host chemistry stage with temperature/RH.
  void validateChemistry(const EnvironmentalStateView &e,
                         const View3D<const double> &s) const {
    validate(e, s, false);
  }
  /// @brief Atomic microphysics and optional accepted coagulation gains.
  void advanceWithCoating(const EnvironmentalStateView &, View3D<double> &,
                          double, View3D<double> *);
  /// @brief Host-owned post-chemistry growth-tail renaming only.
  void renameModes(const EnvironmentalStateView &, View3D<double> &);
  void computeDerivedDiagnostics(const EnvironmentalStateView &,
                                 const View3D<const double> &,
                                 View3D<double> &) override;
  void computeCCN(const EnvironmentalStateView &, const View3D<const double> &,
                  const View1D<const double> &, View4D<double> &) override;
  void computeOptics(const EnvironmentalStateView &,
                     const View3D<const double> &, const View1D<const double> &,
                     View4D<double> &) override;
  std::vector<RadiationBandInfo>
  radiationBandInfo(RadiationBandSet) const override;
  std::vector<RadiationBandComponentAvailability>
  bandAvailability(RadiationBandSet) const override;
  void computeBandOptics(const EnvironmentalStateView &,
                         const View3D<const double> &, RadiationBandSet,
                         std::span<const int>, View4D<double> &) override;
  void computeEmissions(const EnvironmentalStateView &,
                        const EmissionsInputView &, View3D<double> &) override;
  int getSpeciesIndex(const std::string &) const override;
  std::string getSpeciesName(int i) const override {
    return fields_.at(i).name;
  }
  void modalProperties(const EnvironmentalStateView &,
                       const View3D<const double> &, View4D<double> &,
                       View4D<double> *) const override;
  void modalUptake(const EnvironmentalStateView &, const View3D<const double> &,
                   std::span<const ModalVapor>, View4D<double> &,
                   View4D<double> *) const override;
  ModalNucleationTerms nucleationTerms(double) const override;
  double emittedNumber(std::size_t, std::span<const double>) const override;
  void setModalOpticsProvider(
      std::shared_ptr<const IModalOpticsProvider> provider) override {
    optics_provider_ = std::move(provider);
  }
  void setModalBandOpticsProvider(
      std::shared_ptr<const IModalBandOpticsProvider> provider) override {
    band_optics_provider_ = std::move(provider);
  }
};
} // namespace exaero
