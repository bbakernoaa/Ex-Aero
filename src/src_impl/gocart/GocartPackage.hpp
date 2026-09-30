#pragma once
/// @file GocartPackage.hpp
/// @brief GOCART/GEOSmie aerosol package — the primary IAerosolPackage
/// implementation.
///
/// Owns the host-side species parameter set, the single flat device curve
/// pool (all RH/spectral/polarization tables uploaded once at init), and
/// the opaque solver state that the Kokkos kernels read on the hot path.
/// All physics steps (diagnostics, optics, CCN, emissions, attribute
/// queries) are implemented over these members in GocartPackage.cpp and
/// GocartSolver.cpp.
///
/// @section dataflow Curve-pool data flow
/// YAML or GocartConfig -> resolved species params + curve bindings ->
/// one flat @c double pool staged in @c h_curve_pool_ -> single H2D
/// upload -> device views addressed by the offsets/extents recorded in
/// each GocartSpeciesParams. Timestep kernels never touch host memory.
#include <exaero/IAerosolPackage.hpp>
#include <gocart/GocartSpeciesParams.hpp>
#include <string>
#include <vector>

namespace exaero {

struct GocartSolverState;

/// @brief Concrete GOCART aerosol package (GEOSmie MIE-table backed).
///
/// Implements the full IAerosolPackage contract, including the additive
/// GEOSmie attribute surface (computeAttributes / queryAttribute /
/// activation / curve config / moment counts). See IAerosolPackage.hpp
/// for the per-routine parameter and unit contracts overridden here.
class GocartPackage : public IAerosolPackage {
private:
  int num_species_ = 0; ///< Active species count (state axis length).
  std::vector<GocartSpeciesParams> h_species_params_;
  std::vector<std::string>
      species_names_; ///< Parsed config species labels, index-aligned with
                      ///< the state axis.
  /// @brief Host staging for the single flat curve pool uploaded to device
  /// once at init (zero H2D in the timestep loop).
  std::vector<double> h_curve_pool_;

  /// @brief Opaque device-resident solver state (Kokkos views, pool
  /// pointer). Lifetime managed by this package; forward-declared to
  /// keep Kokkos out of this header.
  GocartSolverState *solver_state_ = nullptr;

  /// @brief Shared orchestration for both initialization entry points.
  ///
  /// Order: species params -> store reset -> activation -> curve configs
  /// -> pool wiring -> emissions -> solver-state upload. Reads config
  /// ONLY through the reader so the YAML and structured paths stay
  /// identical (config curves > runtime file > baked-in).
  /// @param reader Config reader abstraction over YAML or GocartConfig.
  /// @throws std::exception with a "FATAL ERROR:" message on any invalid
  /// input; a failed init leaves the package unusable by design.
  void initializeImpl(struct PackageConfigReader &reader);

public:
  /// @brief Construct an empty (uninitialized) package.
  GocartPackage();
  /// @brief Release host staging and the device solver state.
  ~GocartPackage() override;

  /// @brief YAML entry point: parses text, then delegates to
  /// initializeImpl().
  /// @throws std::exception on malformed/missing required config input.
  void initialize(const std::string &config_yaml) override;

  /// @brief Structured (no-YAML) entry point applies an
  /// in-memory GocartConfig through the same orchestration as the YAML path.
  void initialize(const GocartConfig &config) override;

  /// @brief Passive microphysics step (no-op placeholder for GOCART:
  /// transport is owned by the host dynamical core).
  /// @see IAerosolPackage::executeMicrophysics for the contract.
  void executeMicrophysics(const EnvironmentalStateView &env,
                           View3D<double> &state,
                           double delta_time_sec) override;

  /// @brief Device diagnostics kernel launch (see interface contract).
  /// @see IAerosolPackage::computeDerivedDiagnostics.
  void computeDerivedDiagnostics(const EnvironmentalStateView &env,
                                 const View3D<const double> &state,
                                 View3D<double> &diagnostics_out) override;

  /// @brief Map raw CECE emissions onto package species (see contract).
  /// @see IAerosolPackage::computeEmissions.
  void computeEmissions(const EnvironmentalStateView &env,
                        const EmissionsInputView &emissions_in,
                        View3D<double> &emissions_out) override;

  /// @brief Device optics kernel launch (see interface contract).
  /// @see IAerosolPackage::computeOptics.
  void computeOptics(const EnvironmentalStateView &env,
                     const View3D<const double> &state,
                     const View1D<const double> &wavelengths,
                     View4D<double> &optics_out) override;

  /// @brief Device Köhler CCN-activation kernel launch (see contract).
  /// @see IAerosolPackage::computeCCN.
  void computeCCN(const EnvironmentalStateView &env,
                  const View3D<const double> &state,
                  const View1D<const double> &supersaturations,
                  View4D<double> &ccn_out) override;

  /// @brief Dynamic species name -> packed index lookup.
  /// @return Index or -1 when the label is unknown.
  int getSpeciesIndex(const std::string &name) const override;
  /// @brief Dynamic packed index -> species name lookup.
  /// @throws std::out_of_range on an index outside [0, num_species_).
  std::string getSpeciesName(int index) const override;

  // Attribute surface resolve curve mapping + activation through the
  // store.
  /// @brief Device-resident bulk MIE attribute query (see contract).
  /// Curves resolve through the store: config override > runtime file >
  /// baked-in. Slot layout mirrors the C header's band-major contract.
  /// @see IAerosolPackage::computeAttributes.
  void computeAttributes(const EnvironmentalStateView &env,
                         const View3D<const double> &state, int species_index,
                         AttributeCategory category,
                         const View1D<const double> &wavelengths,
                         View3D<double> &attributes_out,
                         View3D<int> *status_out = nullptr) override;
  /// @brief Install config-species curve bindings (resolved at
  /// initialize(); bad config aborts with "FATAL ERROR:").
  void
  setSpeciesCurveConfig(const std::vector<SpeciesCurveConfig> &curves) override;
  /// @brief Set category/species activation and optional runtime MIE file
  /// (a file failing validation aborts initialization, no fallback).
  void
  setAttributeActivation(const std::vector<std::string> &species,
                         int categories_mask,
                         const std::string &runtime_file_path = "") override;
  /// @brief Host-side scalar attribute query with provenance; writes
  /// @p value_out only on an available/interpolated status (never 0).
  AttributeStatus
  queryAttribute(int species_index, AttributeCategory category,
                 int attribute_index, double rh, double wavelength_m,
                 double *value_out,
                 ProvenanceInfo *provenance_out = nullptr) const override;
  /// @brief Report the species' polarized (element, moment) slot counts
  /// from its bound curve (data-sized; zeros when no moments).
  void momentCounts(int species_index, int *num_pol_out,
                    int *num_moment_out) const override;

  /// @brief Active species count (test accessor).
  int get_num_species() const { return num_species_; }
  /// @brief Host species parameter block for @p i (test accessor).
  GocartSpeciesParams get_species_params(int i) const {
    return h_species_params_[i];
  }
};

} // namespace exaero
