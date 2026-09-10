#pragma once
// Internal provenance record for the runtime table store (exaero_impl, no Kokkos).
//
// Mirrors the public POD exaero::ProvenanceInfo but uses std::string for host-side
// composition; it is copied into the POD only at the public query boundary.
// See specs/001-geosmie-lut-attributes/data-model.md "Provenance Record".
#include <exaero/AttributeQuery.hpp>

#include <string>

namespace exaero {

/// @brief Full provenance attached to one surfaced attribute value (FR-004).
struct Provenance {
    std::string species;          ///< GOCART species label.
    std::string attribute;        ///< Attribute name (e.g. "wet_particle_density").
    std::string unit;             ///< Canonical unit string; "1" when dimensionless.
    std::string source_version;   ///< Pinned GEOSmie table version (FR-016); must be non-empty.
    std::string citation;         ///< Literature reference for the source table.
    DeliverySource delivery = DeliverySource::BakedIn; ///< Winning layer (FR-012, R9).
    bool interpolated = false;    ///< True when the value is not a grid point (FR-006).
    AttributeStatus status = AttributeStatus::NotInSource; ///< Availability (FR-008).

    /// @brief Copy into the public POD form for delivery across the query boundary.
    /// @param num_radius Resolved radius-axis length for the queried species.
    /// @param num_rh Resolved RH-axis length for the queried species.
    /// @param num_lambda Resolved wavelength-axis length (0 for microphysical).
    /// @param num_pol Resolved polarized element count (0 unless moments exist, FR-003).
    /// @param num_moment Resolved polarized moment count (0 unless moments exist, FR-003).
    /// @return Populated public provenance record.
    ProvenanceInfo to_public(int num_radius, int num_rh, int num_lambda,
                             int num_pol = 0, int num_moment = 0) const;
};

/// @brief Copy a string into a fixed public POD field with guaranteed NUL termination.
/// @param dst Destination buffer.
/// @param dst_size Destination capacity in bytes.
/// @param src Source string (truncated to fit, never overflowed).
void copy_provenance_field(char* dst, std::size_t dst_size, const std::string& src);

} // namespace exaero
