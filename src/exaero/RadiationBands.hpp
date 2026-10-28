#pragma once
/// @file RadiationBands.hpp
/// @brief Explicit metadata for supported radiation-band coordinates.
#include <array>
#include <cstddef>
#include <stdexcept>
#include <span>
#include <string>
#include <vector>

namespace exaero {

enum class RadiationBandSet : int { RRTMG_30 = 0 };
enum class RadiationBandRegion : int { Shortwave = 0, Longwave = 1 };
enum class RadiationBandStatus : int { Available = 0, NotInSource = 1 };

/// @brief A single source band's immutable coordinate metadata.
struct RadiationBandInfo {
  int source_id;
  RadiationBandRegion region;
  double lower_wavenumber_cm1;
  double upper_wavenumber_cm1;
};

/// @brief Per-component support for every canonical band.
struct RadiationBandComponentAvailability {
  std::size_t component_index;
  std::string component_name;
  std::vector<RadiationBandStatus> status;
};

/// @brief Return the canonical 14 SW + 16 LW RRTMG source-band metadata.
inline const std::array<RadiationBandInfo, 30> &
radiationBandInfo(RadiationBandSet set) {
  if (set != RadiationBandSet::RRTMG_30)
    throw std::invalid_argument("unknown radiation band set");
  static constexpr std::array<RadiationBandInfo, 30> bands{{
      {1, RadiationBandRegion::Shortwave, 2600., 3250.},
      {2, RadiationBandRegion::Shortwave, 3250., 4000.},
      {3, RadiationBandRegion::Shortwave, 4000., 4650.},
      {4, RadiationBandRegion::Shortwave, 4650., 5150.},
      {5, RadiationBandRegion::Shortwave, 5150., 6150.},
      {6, RadiationBandRegion::Shortwave, 6150., 7700.},
      {7, RadiationBandRegion::Shortwave, 7700., 8050.},
      {8, RadiationBandRegion::Shortwave, 8050., 12850.},
      {9, RadiationBandRegion::Shortwave, 12850., 16000.},
      {10, RadiationBandRegion::Shortwave, 16000., 22650.},
      {11, RadiationBandRegion::Shortwave, 22650., 29000.},
      {12, RadiationBandRegion::Shortwave, 29000., 38000.},
      {13, RadiationBandRegion::Shortwave, 38000., 50000.},
      {14, RadiationBandRegion::Shortwave, 820., 2600.},
      {15, RadiationBandRegion::Longwave, 10., 350.},
      {16, RadiationBandRegion::Longwave, 350., 500.},
      {17, RadiationBandRegion::Longwave, 500., 630.},
      {18, RadiationBandRegion::Longwave, 630., 700.},
      {19, RadiationBandRegion::Longwave, 700., 820.},
      {20, RadiationBandRegion::Longwave, 820., 980.},
      {21, RadiationBandRegion::Longwave, 980., 1080.},
      {22, RadiationBandRegion::Longwave, 1080., 1180.},
      {23, RadiationBandRegion::Longwave, 1180., 1390.},
      {24, RadiationBandRegion::Longwave, 1390., 1480.},
      {25, RadiationBandRegion::Longwave, 1480., 1800.},
      {26, RadiationBandRegion::Longwave, 1800., 2080.},
      {27, RadiationBandRegion::Longwave, 2080., 2250.},
      {28, RadiationBandRegion::Longwave, 2250., 2380.},
      {29, RadiationBandRegion::Longwave, 2380., 2600.},
      {30, RadiationBandRegion::Longwave, 2600., 3250.},
  }};
  return bands;
}

/// @brief Select zero-based canonical bands, preserving caller order.
inline std::vector<RadiationBandInfo>
selectRadiationBands(RadiationBandSet set, std::span<const int> band_indices) {
  const auto &all = radiationBandInfo(set);
  if (band_indices.empty())
    return {all.begin(), all.end()};
  std::vector<RadiationBandInfo> selected;
  selected.reserve(band_indices.size());
  std::array<bool, 30> seen{};
  for (const int index : band_indices) {
    if (index < 0 || index >= static_cast<int>(all.size()))
      throw std::invalid_argument("radiation band index is out of range");
    if (seen[static_cast<std::size_t>(index)])
      throw std::invalid_argument("radiation band indices must be unique");
    seen[static_cast<std::size_t>(index)] = true;
    selected.push_back(all[static_cast<std::size_t>(index)]);
  }
  return selected;
}

} // namespace exaero
