#pragma once
// Native portable MIE curve file loader (exaero_impl — no NetCDF, no Kokkos).
//
// Reuses the SPEC-CRTM-004 ASCII stream-parse pattern: '#'/'!' comment lines,
// explicit dimension declarations, fail-fast on
// truncation/malformed/unit-mismatch, tolerance of Fortran exponent forms
// (1.0d-05). The runtime file EXTENDS or OVERRIDES the baked-in curves
// (FR-012); any validation failure aborts initialization with a "FATAL ERROR:"
// diagnostic and no silent fallback (FR-009).
//
// Schema (whitespace-separated float64 records):
//   # @format exaero_mie_v1
//   # @source_version <non-empty string>          (FR-016)
//   # @citation <string>
//   # @species <LABEL>                             (begins one species block)
//   # @axis <name> <units> <n> [n values...] (radius|rh|lambda|pol|moment) #
//   @field <name> <dims> <units>                 (dims = comma-joined axis
//   names) <value rows: product(dims lengths) doubles, may span lines>
#include <loader/MieTableStore.hpp>

#include <string>
#include <vector>

namespace exaero {

/// @brief Loader for the native portable MIE curve format.
class MieFileLoader {
public:
  /// @brief Parse and fully validate one portable file into curves.
  /// @param path File path; sanitized here (no traversal, bounded length) and
  /// again at
  ///             the C boundary (Security baseline).
  /// @return Parsed species curves with delivery = RuntimeFile.
  /// @throws std::runtime_error "FATAL ERROR: ..." on any
  /// schema/truncation/unit failure.
  static std::vector<SpeciesCurve> load(const std::string &path);

  /// @brief Path validation shared with the C boundary.
  /// @return true when the path is acceptable (non-empty, no traversal,
  /// bounded).
  static bool path_is_safe(const std::string &path);
};

} // namespace exaero
