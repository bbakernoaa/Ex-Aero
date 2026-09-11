#pragma once
// Public attribute-query descriptor surface for EX-aero.
//
// MUST NOT include <Kokkos_Core.hpp> (research R4/R8, contract §1, SC-002):
// this header is part of the zero-dependency public interface.
//
// Curve-mapping model (ADR-003 / research R9): a table is a field over
// dynamic-length axes (radius, rh, lambda[, pol, moment]). Axis lengths are
// DATA, never code: no bin count, RH-grid length, or band count is fixed by
// anything declared here.
#include <cstddef>
#include <string>
#include <vector>

namespace exaero {

/// @brief Category of a surfaced GEOSmie MIE table attribute.
enum class AttributeCategory : int {
  Microphysical = 0,   ///< No wavelength axis; indexed by (radius, rh).
  SpectralOptical = 1, ///< Indexed by (radius, rh, lambda).
  PolarizedMoment = 2, ///< Indexed by (radius, rh, lambda, pol, moment).
  NUM_CATEGORIES = 3
};

/// @brief Availability outcome of a single attribute query (FR-008: never a
/// silent 0).
enum class AttributeStatus : int {
  Available = 0,     ///< Exact grid point, baked-in data.
  AvailableFile = 1, ///< Exact grid point, supplied by a runtime file (FR-012).
  Interpolated = 2,  ///< Between grid points; see the provenance `interpolated`
                     ///< flag (FR-006).
  NotActivated = 3,  ///< Species/category deselected by the operator (FR-010).
  NotInSource =
      4, ///< Species lacks this category in the pinned tables (FR-008).
  AvailableConfig = 5 ///< Exact point on a config override curve (ADR-003 R9).
};

/// @brief Where a surfaced value came from, in precedence order (FR-012, R9).
enum class DeliverySource : int {
  BakedIn = 0,     ///< Generated at build time from the pinned snapshot.
  RuntimeFile = 1, ///< Supplied by an optional runtime data file.
  Config = 2 ///< Supplied by the runtime configuration (overrides file and
             ///< baked-in).
};

/// @brief Bitmask helper: a category's bit within an activation mask (FR-010).
/// @param category Attribute category to encode.
/// @return Single-bit mask for that category.
constexpr int attribute_category_bit(AttributeCategory category) noexcept {
  return 1 << static_cast<int>(category);
}

/// @brief Index codes for the Microphysical category (no wavelength axis).
/// Units are normative (data-model.md); a query result is meaningless without
/// them.
namespace microphysical_indices {
constexpr int WET_PARTICLE_DENSITY = 0; ///< [kg m^-3]
constexpr int GROWTH_FACTOR =
    1; ///< fraction (wet/dry radius), >= 1, monotone up in rh
constexpr int EFFECTIVE_RADIUS = 2; ///< [m]
constexpr int MASS_MEAN_RADIUS = 3; ///< [m]
constexpr int BIN_LOWER_RADIUS = 4; ///< [m]
constexpr int BIN_UPPER_RADIUS = 5; ///< [m]
constexpr int VOLUME_PER_MASS = 6;  ///< [m^3 kg^-1]
constexpr int AREA_PER_MASS = 7;    ///< [m^2 kg^-1]
constexpr int PARTICLE_MASS =
    8; ///< [kg] (source `rMass`; data-model parity, added at implement)
constexpr int NUM_ATTRIBUTES = 9;
} // namespace microphysical_indices

/// @brief Index codes for the SpectralOptical category (indexed by radius, rh,
/// lambda).
namespace spectral_optical_indices {
constexpr int EXTINCTION_EFFICIENCY = 0;    ///< dimensionless, >= 0
constexpr int SCATTERING_EFFICIENCY = 1;    ///< dimensionless, [0, qext]
constexpr int ABSORPTION_EFFICIENCY = 2;    ///< dimensionless, >= 0
constexpr int MASS_EXTINCTION = 3;          ///< [m^2 (kg dry mass)^-1]
constexpr int MASS_SCATTERING = 4;          ///< [m^2 (kg dry mass)^-1]
constexpr int MASS_BACKSCATTER = 5;         ///< [m^2 (kg dry mass)^-1 sr^-1]
constexpr int LIDAR_RATIO = 6;              ///< [sr], guarded division (FR-007)
constexpr int ASYMMETRY_FACTOR = 7;         ///< fraction, [-1, 1] (SC-005)
constexpr int SINGLE_SCATTERING_ALBEDO = 8; ///< fraction, [0, 1] (SC-005)
constexpr int REFRACTIVE_INDEX_REAL = 9;    ///< dimensionless, > 0 (wet)
constexpr int REFRACTIVE_INDEX_IMAG = 10;   ///< dimensionless, >= 0 (wet)
constexpr int NUM_ATTRIBUTES = 11;
} // namespace spectral_optical_indices

/// @brief Index codes for the PolarizedMoment category.
/// Element ordering is normative and fixed by the source tables; the moment
/// count is per-species DATA (301/751/2001 in the pinned release) and never
/// hardcoded.
namespace polarized_moment_indices {
constexpr int ELEMENT_STRIDE = 6; ///< Ordering: P11, P12, P33, P34, P22, P44.
constexpr int PHASE_FUNCTION_MOMENT =
    0; ///< Single attribute, addressed by (element, moment).
constexpr int NUM_ATTRIBUTES = 1;
} // namespace polarized_moment_indices

/// @brief Provenance attached to every surfaced value (FR-004, FR-016).
/// POD and Kokkos-free so it crosses the C/Fortran boundary by value.
struct ProvenanceInfo {
  char species[16] = {
      0};              ///< GOCART species label (DU, SS, SU, BC, OC, BR, NI).
  char unit[32] = {0}; ///< Canonical unit string; "1" for dimensionless.
  char source_version[96] = {0}; ///< Pinned table version (FR-016).
  char citation[192] = {
      0}; ///< Literature reference (Colarco/Dubovik/Chin/Bian).
  int delivery_source = 0; ///< DeliverySource value (FR-012, R9).
  int interpolated = 0;    ///< 1 when the value is not a grid point (FR-006).
  int status = 3;          ///< AttributeStatus value for this result (FR-008).
  int num_radius =
      0;          ///< Resolved radius-axis length for this species (dynamic).
  int num_rh = 0; ///< Resolved RH-axis length for this species (dynamic).
  int num_lambda = 0; ///< Resolved wavelength-axis length (0 when N/A).
  int num_pol =
      0; ///< Polarized element count (0 unless moments exist, FR-003).
  int num_moment =
      0; ///< Polarized moment count (0 unless moments exist, FR-003).
};

/// @brief Runtime curve binding for one configuration species (research R9,
/// ADR-003). All arrays are non-owning views over caller memory with explicit
/// counts: no fixed-size shapes cross the public boundary, so a config may
/// declare any bin structure.
struct SpeciesCurveConfig {
  std::string species_name; ///< Config-side species name (arbitrary, e.g.
                            ///< "dust_mode_1").
  std::string source_label; ///< Curve to bind (e.g. "DU"); empty => species has
                            ///< no curve.
  std::vector<double> radius_nodes; ///< Config-owned bin centres [m]; empty =>
                                    ///< keep source axis.
  std::string interpolate =
      "linear"; ///< Named strategy behind the seam: "linear" (seam: "tension").
  int solver_radius_node =
      0; ///< Node used by the device hot path (in [0, size)).

  /// @brief One per-attribute override row at the config's own radius nodes.
  struct Override {
    AttributeCategory category = AttributeCategory::Microphysical;
    int attribute_index =
        0; ///< Index code from the category's *_indices namespace.
    std::vector<double> values; ///< Length == radius_nodes.size().
  };
  std::vector<Override>
      overrides; ///< Config values that become this species' curve (R9).
};

} // namespace exaero
