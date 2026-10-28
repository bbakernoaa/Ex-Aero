#pragma once
/// @file OpticsBandCoupling.hpp
/// @brief Additive band optics surface that does not change package vtables.
#include <exaero/IAerosolPackage.hpp>
#include <exaero/RadiationBands.hpp>
#include <span>

namespace exaero {

class IOpticsBandCoupling {
public:
  virtual ~IOpticsBandCoupling() = default;
  /// Metadata order follows the selected set. Empty band_indices selects all
  /// bands; otherwise indices are zero-based positions in that set and retain
  /// caller order. Output shape is (cell, level, selected-band, optical-index)
  /// over host-accessible View4D storage. Layer diagnostics occupy every
  /// level; column AOT values are stored at level zero. Implementations
  /// validate the complete request first and leave output unchanged on error.
  virtual std::vector<RadiationBandInfo>
  radiationBandInfo(RadiationBandSet) const = 0;
  virtual std::vector<RadiationBandComponentAvailability>
  bandAvailability(RadiationBandSet) const = 0;
  virtual void computeBandOptics(const EnvironmentalStateView &,
                                 const View3D<const double> &,
                                 RadiationBandSet, std::span<const int>,
                                 View4D<double> &) = 0;
};

/// @brief Checked access to the optional band optics extension interface.
inline IOpticsBandCoupling &opticsBandCoupling(IAerosolPackage &package) {
  auto *coupling = dynamic_cast<IOpticsBandCoupling *>(&package);
  if (!coupling)
    throw std::logic_error("FATAL ERROR: package has no radiation-band optics");
  return *coupling;
}

} // namespace exaero
