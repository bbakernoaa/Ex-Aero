#pragma once
/// @file CloudIceCoupling.hpp
/// @brief Host-owned cloud and ice candidate diagnostics for bulk and modal.
#include <exaero/BulkRemoval.hpp>
#include <exaero/CloudIceConfig.hpp>
#include <exaero/IAerosolPackage.hpp>

namespace exaero {
namespace cloud_ice_properties {
inline constexpr int MASS = 0;
inline constexpr int NUMBER = 1;
inline constexpr int DRY_NUMBER_MEDIAN_RADIUS = 2;
inline constexpr int DRY_SCAVENGING_RADIUS = 3;
inline constexpr int WET_SCAVENGING_RADIUS = 4;
inline constexpr int LOWER_RADIUS = 5;
inline constexpr int UPPER_RADIUS = 6;
inline constexpr int EFFECTIVE_RADIUS = 7;
inline constexpr int TYPE_WEIGHT = 8;
inline constexpr int SIZE_WEIGHT = 9;
inline constexpr int ELIGIBLE_NUMBER = 10;
inline constexpr int HETEROGENEOUS_FREEZING = 11;
inline constexpr int COUNT = 12;
} // namespace cloud_ice_properties

enum class CloudIceStatus : int {
  Absent = -2,
  Empty = -1,
  NotConfigured = 0,
  Configured = 1
};

/// @brief Additive bulk metadata; no IAerosolPackage vtable change.
class IBulkCloudIceMetadata {
public:
  virtual ~IBulkCloudIceMetadata() = default;
  virtual std::vector<IceCategory>
  bulkIceCategories(std::size_t species_index) const = 0;
  virtual std::vector<IcePopulationRule>
  bulkIceRules(std::size_t species_index) const = 0;
};

/// @brief Common population diagnostics; package selects bulk or modal.
/// @details Same host-owned reservoir, shape, status and atomicity contract as
/// bulkCloudIce/modalCloudIce; population count comes from hostCouplingInfo.
void cloudIce(IAerosolPackage &, const EnvironmentalStateView &,
              const View3D<const double> &ambient,
              const View3D<const double> *cloudborne, IceCategory,
              bool preexisting_ice, View4D<double> &ambient_out,
              View4D<double> &cloudborne_out, View3D<int> &ambient_status,
              View3D<int> &cloudborne_status);
/// @brief Common number/mass activation fractions (cell,level,P,2*num_S).
void activationFractions(IAerosolPackage &, const EnvironmentalStateView &,
                         const View3D<const double> &,
                         const View1D<const double> &, View4D<double> &);

/// @brief Diagnose both populations from host-owned packed states.
/// @details Output slabs are (cell,level,species-or-mode,COUNT); status slabs
/// are (cell,level,species-or-mode). Missing cloudborne input gives zero output
/// and Absent status. No state is mutated.
void bulkCloudIce(IAerosolPackage &package, const EnvironmentalStateView &env,
                  const View3D<const double> &ambient,
                  const View3D<const double> *cloudborne, IceCategory category,
                  bool preexisting_ice, View4D<double> &ambient_out,
                  View4D<double> &cloudborne_out, View3D<int> &ambient_status,
                  View3D<int> &cloudborne_status);
void modalCloudIce(IAerosolPackage &package, const EnvironmentalStateView &env,
                   const View3D<const double> &ambient,
                   const View3D<const double> *cloudborne, IceCategory category,
                   bool preexisting_ice, View4D<double> &ambient_out,
                   View4D<double> &cloudborne_out, View3D<int> &ambient_status,
                   View3D<int> &cloudborne_status);

/// @brief Interstitial number and mass activation fractions for each S.
/// @details Output is (cell,level,species-or-mode,2*num_S), with number
/// fraction at 2*S and mass fraction at 2*S+1. Cloudborne aerosol is excluded.
void bulkActivationFractions(IAerosolPackage &package,
                             const EnvironmentalStateView &env,
                             const View3D<const double> &ambient,
                             const View1D<const double> &supersaturation,
                             View4D<double> &fractions);
void modalActivationFractions(IAerosolPackage &package,
                              const EnvironmentalStateView &env,
                              const View3D<const double> &ambient,
                              const View1D<const double> &supersaturation,
                              View4D<double> &fractions);
} // namespace exaero
