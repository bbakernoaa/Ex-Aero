#pragma once
/// @file BulkRemoval.hpp
/// @brief Explicit host removal utilities for mass-only GOCART bulk aerosol.
#include <exaero/IAerosolPackage.hpp>
#include <span>

namespace exaero {
/// @brief Configured fixed dry distribution for one bulk species.
struct BulkSpeciesSize {
  double dry_particle_diameter = 0; ///< Caller-anchored dry diameter [m].
  double dry_density = 0;           ///< Material density [kg/m3].
  double lognormal_dg = 0;          ///< Number-median dry diameter [m].
  double lognormal_sigma =
      0; ///< Geometric standard deviation; 1 is a point size.
  double dry_radius_lower = 0;     ///< Physical lower radius [m]; 0 if absent.
  double dry_radius_upper = 0;     ///< Physical upper radius [m]; 0 if absent.
  double dry_effective_radius = 0; ///< [m]; 0 derives from lognormal moments.
  double hygroscopicity = 0; ///< Kappa used for table-independent wet sizing.
};

/// @brief Additive metadata surface for packages with mass-only bulk fields.
class IBulkRemovalMetadata {
public:
  virtual ~IBulkRemovalMetadata() = default;
  virtual BulkSpeciesSize bulkSpeciesSize(std::size_t species_index) const = 0;
};

/// @brief Host-owned route between configured bulk species indices.
struct BulkResuspensionRoute {
  std::size_t source = 0, destination = 0;
};

/// @brief Sentinel in size_status for a species with zero ambient mass.
inline constexpr int BULK_SIZE_EMPTY = -1;
inline constexpr int BULK_SIZE_CONFIGURED = 0;
/// @brief Availability of independently configured physical size values.
enum class BulkRadiusStatus : int {
  Unavailable = 0,
  Configured = 1,
  Derived = 2
};

/// @brief Diagnose particles/m3 for each mass-only bulk species.
/// @details Uses the same fixed lognormal mass-to-number conversion as GOCART
/// diagnostics: dry density, geometric mean diameter, and lognormal width.
void bulkFixedNumber(IAerosolPackage &package,
                     const View3D<const double> &state, View3D<double> &number);

/// @brief Query physical dry lower, upper, and effective radii by species.
/// @details radii and status have shape (species,3), in that order. Each
/// value has its own BulkRadiusStatus; unavailable values remain zero. This
/// never queries optical or radiation tables.
void bulkReferenceRadii(IAerosolPackage &package, View2D<double> &radii,
                        View2D<int> &status);

/// @brief Query dry/wet diameter (cell,level,species,2) from current bulk mass.
/// @details size_status is BULK_SIZE_EMPTY or BULK_SIZE_CONFIGURED. Wet growth
/// uses configured kappa and RH, independently of optical or radiation tables.
void bulkScavengingSizes(IAerosolPackage &package,
                         const View2D<const double> &relative_humidity,
                         const View3D<const double> &state,
                         View4D<double> &diameters, View3D<int> &size_status);

/// @brief Route nonnegative returned material rates from an input snapshot.
void redistributeBulkResuspension(IAerosolPackage &package,
                                  const View3D<const double> &mass_rates,
                                  std::span<const BulkResuspensionRoute> routes,
                                  View3D<double> &routed_rates);

/// @brief Rebin selected-species integrated deposited kg/m2 by dry diameter.
void rebinBulkDeposition(IAerosolPackage &package,
                         const View3D<const double> &deposited_budget,
                         std::size_t species_index,
                         const View1D<const double> &diameter_edges,
                         View3D<double> &bin_mass, View3D<double> &tail_mass);
} // namespace exaero
