#pragma once
// Host-side MIE table store (exaero_impl — NO Kokkos, NO mdspan on the solver
// path).
//
// SpheroidDatabase-style singleton holding, per configured species, a
// SpeciesCurve: a float64 field over dynamic-length axes (radius, rh, lambda[,
// pol, moment]) merged from three layers with precedence config > runtime-file
// > baked-in. Derived attributes are computed exactly
// at grid points from the raw baked arrays.
//
// Nothing here hardcodes a bin count, RH-grid length, or band count: axis
// lengths are data read from the curve. Validation aborts with a
// "FATAL ERROR:" diagnostic and no silent fallback ; absent data yields
// an explicit status, never a silent 0.
#include <exaero/AttributeQuery.hpp>
#include <loader/Provenance.hpp>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace exaero {

/// @brief One attribute field of a curve, stored flat in (radius, rh[,
/// lambda][, pol, moment]) order.
struct CurveField {
  std::vector<double>
      values; ///< float64, layout: radius-major (bin, rh[, lambda]).
  AttributeCategory category =
      AttributeCategory::Microphysical; ///<\< Optical, microphysical, or
                                        ///< mixing-state class (drives status
                                        ///< semantics).
  std::string unit; ///< canonical unit string.
  int rank =
      1; ///< 1 = (radius); 2 = (radius,rh); 3 = (radius,rh,lambda); 5 = pmom.
  DeliverySource delivery =
      DeliverySource::BakedIn; ///< winning layer for this field.
  std::vector<std::string>
      dims; ///< axis names in order (loader communicates shape;).
};

/// @brief A species' merged curve: dynamic axes + attribute fields.
struct SpeciesCurve {
  std::string source_label;   ///< bound source (e.g. "DU"); empty => no curve.
  std::string species_name;   ///< config-side name.
  std::string source_version; ///< provenance version string ; non-empty
                              ///< required.
  std::string citation;       ///< provenance citation string.
  std::vector<double> radius; ///< [n_radius] representative radii [m] (may be
                              ///< non-strict: modes).
  std::vector<double> rh;     ///< [n_rh] RH fractions [0,0.99].
  std::vector<double> lambda; ///< [n_lambda] band indices or wavelengths; empty
                              ///< for microphysical-only.
  std::vector<double>
      pol; ///< [n_pol] polarized element indices (P11..P44); empty unless pmom.
  std::vector<double>
      moment; ///< [n_moment] Legendre/GSF moment indices; empty unless pmom.
  std::map<std::string, CurveField> fields; ///< attribute name -> field.
  bool radius_resamplable = false; ///< true iff radius strictly increasing.
  bool config_modified =
      false; ///< radius resampled or overridden by config (status 5).
  int solver_radius_node = 0; ///< Default size-bin the device hot path queries.

  int n_radius() const {
    return static_cast<int>(radius.size());
  } ///<\< Size-axis extent (number of representative radii).
  int n_rh() const {
    return static_cast<int>(rh.size());
  } ///<\< Hygroscopicity-axis extent (number of RH bins).
  int n_lambda() const {
    return static_cast<int>(lambda.size());
  } ///<\< Spectral-axis extent; 0 for microphysical-only curves.
  int n_pol() const { return static_cast<int>(pol.size()); } ///< 0 unless pmom.
  int n_moment() const {
    return static_cast<int>(moment.size());
  } ///< 0 unless pmom.
};

/// @brief Serialized flat device-uploadable view of all activated curves.
/// A single contiguous pool + per-species descriptor; the
/// solver consumes raw pointers.
struct CurvePoolDescriptor {
  int species_index = -1; ///< State-axis index of this species (-1 unset).
  std::string
      source_label; ///< Bound Mie source (e.g. "DU"); empty => no curve.
  int curve_offset =
      0; ///< index into the double pool where this species' fields begin.
  int n_radius = 0; ///< Size-axis extent of this species' curve.
  int n_rh = 0;     ///< RH-axis extent of this species' curve.
  int n_lambda = 0; ///< Spectral-axis extent; 0 for microphysical-only.
  // Field sub-offsets within this species' block, per canonical serialized
  // field order.
  std::vector<std::int64_t> field_offsets; ///< -1 where a field is absent.
  std::vector<std::string> field_names;    ///< Canonical serialized field order
                                        ///< (index-aligned with field_offsets).
};

/// @brief Host-side store of merged GEOSmie MIE curves.
class MieTableStore {
public:
  /// @brief Process-wide singleton (SpheroidDatabase-style).
  static MieTableStore &instance();

  MieTableStore(const MieTableStore &) = delete;
  MieTableStore &operator=(const MieTableStore &) = delete;

  /// @brief Load baked-in curves from the generated header into the store
  /// (idempotent). Aborts with FATAL ERROR on any validation failure in the
  /// baked data.
  void ensure_baked_in_loaded();

  /// @brief Reset all merged curves to the baked-in baseline (called at package
  /// init).
  void reset_to_baked_in();

  /// @brief Merge an optional runtime file's curves (extends/overrides baked).
  /// @param path Runtime portable data-file path (already sanitized by the
  /// caller).
  /// @throws std::runtime_error prefixed "FATAL ERROR:" on validation failure.
  void apply_runtime_file(const std::string &path);

  /// @brief Resolve config curve bindings: resample radius onto config nodes,
  /// apply
  /// overrides, and register each config species (precedence config >
  /// file > baked).
  /// @param curves Per-species bindings from the configuration.
  /// @throws std::runtime_error "FATAL ERROR:" on unknown source/strategy or
  /// bad nodes.
  void apply_curve_config(const std::vector<SpeciesCurveConfig> &curves);

  /// @brief Activate a category set + species list; deactivate others.
  /// @param species Species names to activate (empty => all currently bound).
  /// @param categories_mask Bitmask of (1 << AttributeCategory).
  void set_activation(const std::vector<std::string> &species,
                      int categories_mask);

  /// @brief Is a (species, category) pair activated?
  bool is_activated(const std::string &species,
                    AttributeCategory category) const;

  /// @brief Query one attribute at (rh, wavelength) with grid/interpolated
  /// resolution.
  /// @param species Config species name.
  /// @param category Attribute category.
  /// @param attribute_index Index code within the category.
  /// @param rh Relative humidity fraction (clamped to [0,0.99]).
  /// @param wavelength_m Wavelength/band (ignored for Microphysical).
  /// @param radius_index Size-bin/mode index (0-based; <0 selects the first
  /// bin).
  /// @param value_out Receives the value when available/interpolated.
  /// @param prov_out Receives provenance when non-null.
  /// @return Availability status; never a silent 0.
  AttributeStatus query(const std::string &species, AttributeCategory category,
                        int attribute_index, double rh, double wavelength_m,
                        int radius_index, double *value_out,
                        Provenance *prov_out) const;

  /// @brief Serialize all activated curves into one flat pool + descriptors.
  /// @param pool_out Contiguous float64 pool (concatenation of activated field
  /// values).
  /// @param desc_out Per-species descriptors (offsets into pool_out).
  /// @param species_order Which config species to serialize, in solver-index
  /// order.
  void serialize_activated(const std::vector<std::string> &species_order,
                           std::vector<double> &pool_out,
                           std::vector<CurvePoolDescriptor> &desc_out) const;

  /// @brief Canonical unit string for a field name (unit-mismatch
  /// check).
  static const char *unit_for_field(const std::string &name);

  /// @brief Number of species currently bound to a curve.
  int num_bound_species() const { return static_cast<int>(curves_.size()); }

  /// @brief Access a resolved curve (for tests / wiring); nullptr if absent.
  const SpeciesCurve *find_curve(const std::string &species) const;

  /// @brief Canonical serialized field-name order (shared by serialize + solver
  /// reader).
  static const std::vector<std::string> &serialized_field_order();

  ~MieTableStore();

private:
  MieTableStore() = default;

  // Build the derived microphysical + spectral fields for a curve from its raw
  // arrays.
  static void derive_fields(SpeciesCurve &curve);

  // Resample every rank-2/rank-3 field along radius onto new nodes (
  // radius-only).
  static void resample_radius(SpeciesCurve &curve,
                              const std::vector<double> &nodes);

  // Validate one merged curve (monotone axes where resamplable, bounds, units,
  // version).
  void validate_curve(const SpeciesCurve &curve) const;

  bool baked_loaded_ = false;
  int activation_mask_ = 0;                    // category bitmask
  std::vector<std::string> activated_species_; // empty => all bound
  std::map<std::string, SpeciesCurve> curves_; // keyed by config species_name
};

} // namespace exaero
