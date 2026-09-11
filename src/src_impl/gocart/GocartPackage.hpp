#pragma once
#include <exaero/IAerosolPackage.hpp>
#include <gocart/GocartSpeciesParams.hpp>
#include <string>
#include <vector>

namespace exaero {

struct GocartSolverState;

class GocartPackage : public IAerosolPackage {
private:
  int num_species_ = 0;
  std::vector<GocartSpeciesParams> h_species_params_;
  std::vector<std::string>
      species_names_; // Stores parsed YAML species names dynamically
  // Host staging for the single flat curve pool uploaded to device once
  // (ADR-003 R10).
  std::vector<double> h_curve_pool_;
  GocartSolverState *solver_state_ = nullptr;

  // Shared orchestration for both entry points (design 2026-09-10, Approach A):
  // species params -> store reset -> activation -> curve configs -> pool wiring
  // -> emissions -> solver-state upload. Reads config ONLY through the reader.
  void initializeImpl(struct PackageConfigReader &reader);

public:
  GocartPackage();
  ~GocartPackage() override;

  void initialize(const std::string &config_yaml) override;

  /// @brief Structured (no-YAML) entry point (design 2026-09-10): applies an
  /// in-memory GocartConfig through the same orchestration as the YAML path.
  void initialize(const GocartConfig &config) override;

  void executeMicrophysics(const EnvironmentalStateView &env,
                           View3D<double> &state,
                           double delta_time_sec) override;

  void computeDerivedDiagnostics(const EnvironmentalStateView &env,
                                 const View3D<const double> &state,
                                 View3D<double> &diagnostics_out) override;

  void computeEmissions(const EnvironmentalStateView &env,
                        const EmissionsInputView &emissions_in,
                        View3D<double> &emissions_out) override;

  void computeOptics(const EnvironmentalStateView &env,
                     const View3D<const double> &state,
                     const View1D<const double> &wavelengths,
                     View4D<double> &optics_out) override;

  void computeCCN(const EnvironmentalStateView &env,
                  const View3D<const double> &state,
                  const View1D<const double> &supersaturations,
                  View4D<double> &ccn_out) override;

  // Dynamic Species-to-Index mapping query APIs
  int getSpeciesIndex(const std::string &name) const override;
  std::string getSpeciesName(int index) const override;

  // Attribute surface (ADR-003): resolve curve mapping + activation through the
  // store.
  void computeAttributes(const EnvironmentalStateView &env,
                         const View3D<const double> &state, int species_index,
                         AttributeCategory category,
                         const View1D<const double> &wavelengths,
                         View3D<double> &attributes_out,
                         View3D<int> *status_out = nullptr) override;
  void
  setSpeciesCurveConfig(const std::vector<SpeciesCurveConfig> &curves) override;
  void
  setAttributeActivation(const std::vector<std::string> &species,
                         int categories_mask,
                         const std::string &runtime_file_path = "") override;
  AttributeStatus
  queryAttribute(int species_index, AttributeCategory category,
                 int attribute_index, double rh, double wavelength_m,
                 double *value_out,
                 ProvenanceInfo *provenance_out = nullptr) const override;
  void momentCounts(int species_index, int *num_pol_out,
                    int *num_moment_out) const override;

  // Accessors for testing
  int get_num_species() const { return num_species_; }
  GocartSpeciesParams get_species_params(int i) const {
    return h_species_params_[i];
  }
};

} // namespace exaero
