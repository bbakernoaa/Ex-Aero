#pragma once
/// @file IAerosolPackage.hpp
/// @brief Abstract aerosol microphysics/optics interface for EX-aero.
///
/// This is the single zero-dependency public contract every aerosol
/// package (GOCART, the MAM4xx wrapper, test doubles) implements. It
/// defines the non-owning multi-dimensional views that cross the
/// C++/Fortran/CCPP boundary and the lifecycle + compute steps a host
/// model drives each physics timestep.
///
/// @note Kokkos/mdspan-define-free mandate (ADR-001): this header may
/// pull in the <experimental/mdspan> backport only under the forced
/// @c exaero_mdspan namespace so it never collides with a system Kokkos
/// build. No heavy dependency is exposed to host applications.
///
/// @note Memory layout: all multi-dimensional views are column-major
/// (@c layout_left) so a Fortran host array is consumed zero-copy with
/// no transpose. The normative index order is written explicitly on
/// every routine, e.g. (cell, level, species).
#include <exaero/AerosolIndices.hpp>
#include <exaero/AttributeQuery.hpp>
#include <exaero/GocartConfig.hpp>
#include <stdexcept>
#include <string>
#include <vector>

// Force backport to compile under exaero_mdspan namespace to prevent
// redefinition conflicts with system Kokkos
/// @brief Namespace override for the mdspan backport: all std::mdspan
/// symbols are emitted under @c exaero_mdspan so the bundled backport can
/// coexist with a system Kokkos that vendors its own std::mdspan (ADR-001
/// zero-dependency public-header rule).
#define MDSPAN_IMPL_STANDARD_NAMESPACE exaero_mdspan
#include <experimental/mdspan>

namespace exaero {

/// @brief 1D non-owning view over contiguous data (rank-1 needs no
/// layout tag). Used for per-band wavelength and supersaturation axes.
/// @tparam T Element type; use @c const T for read-only inputs.
template <typename T>
using View1D =
    exaero_mdspan::mdspan<T,
                          exaero_mdspan::extents<size_t, std::dynamic_extent>>;

/// @brief 2D column-major (layout_left) non-owning view.
///
/// The 2D/3D/4D layouts are column-major for zero-copy Fortran/CCPP
/// memory alignment. Normative index order for environmental fields is
/// (cell, level).
/// @tparam T Element type; use @c const T for read-only inputs.
template <typename T>
using View2D = exaero_mdspan::mdspan<
    T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent>,
    exaero_mdspan::layout_left>;

/// @brief 3D column-major non-owning view.
///
/// Normative index order is (cell, level, species) for aerosol state
/// and (cell, level, diagnostic_index) for derived output.
/// @tparam T Element type; use @c const T for read-only inputs.
template <typename T>
using View3D = exaero_mdspan::mdspan<
    T,
    exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent,
                           std::dynamic_extent>,
    exaero_mdspan::layout_left>;

/// @brief 4D column-major non-owning view.
///
/// Used for the optical/CCN spectra: normative order is
/// (cell, level, band_index, optics_index) or (cell, level, S_index, ...).
/// @tparam T Element type; use @c const T for read-only inputs.
template <typename T>
using View4D = exaero_mdspan::mdspan<
    T,
    exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent,
                           std::dynamic_extent, std::dynamic_extent>,
    exaero_mdspan::layout_left>;

/// @brief Read-only environmental boundary conditions for one timestep.
///
/// Every 2D field is a (cell, level) column-major view supplied by the
/// host model. The solver never mutates these; they drive wet growth,
/// the RH-dependent optical tables, and vertical layer integration
/// (@f$ \Delta z @f$).
struct EnvironmentalStateView {
  View2D<const double> temperature; ///< [K] air temperature (cell, level).
  View2D<const double> pressure;    ///< [Pa] air pressure (cell, level).
  View2D<const double> air_density; ///< [kg m^-3] air density (cell, level).
  View2D<const double>
      relative_humidity;                ///< RH [fraction 0-1] (cell, level);
                                        ///< clamped to <= 0.99 in the solvers.
  View2D<const double> layer_thickness; ///< Grid layer vertical thickness
                                        ///< @f$ \Delta z @f$ [m] (cell, level).
};

/// @brief Physical interpretation of the incoming emissions flux array.
///
/// Selects how computeEmissions() converts a raw flux into per-species
/// mass/number emission rates: either a concentration tendency applied
/// directly, or an area flux divided by the layer thickness.
enum class FluxType {
  MASS_CONCENTRATION_RATE, ///< Tendency [kg m^-3 s^-1]: applied directly.
  AREA_FLUX ///< Vertical area flux [kg m^-2 s^-1]: divided by @f$ \Delta z @f$.
};

/// @brief Raw emissions input: the flux field plus its unit convention.
struct EmissionsInputView {
  View3D<const double> flux; ///< Raw flux (cell, level, raw_species_index).
  FluxType flux_type;        ///< How to interpret @c flux.
};

/// @brief Abstract base class for an EX-aero aerosol package.
///
/// Defines the full physics contract a host model (CCPP component,
/// Fortran driver, or C++ test) drives each timestep: initialization,
/// passive microphysics, derived diagnostics, optics, cloud activation,
/// and emissions mapping, plus the GEOSmie MIE attribute surface.
///
/// @note Additive/backward-compatibility pattern: virtuals added after
/// the initial release carry a default implementation (throw or no-op)
/// so out-of-tree packages (the MAM4xx wrapper, test doubles) compile
/// and behave unchanged until they opt in.
class IAerosolPackage {
public:
  /// @brief Virtual destructor: packages are owned polymorphically.
  virtual ~IAerosolPackage() = default;

  /// @brief Initialize the package from a YAML configuration string.
  ///
  /// Parses species, activation, and emissions-mapping blocks and pins
  /// the MIE table data. Value precedence: config curves > runtime file
  /// > baked-in tables.
  /// @param config_yaml YAML document text (not a file path).
  /// @throws std::exception on a malformed or missing required input,
  ///         surfaced through the C/CCPP layer as a @c "FATAL ERROR:".
  virtual void initialize(const std::string &config_yaml) = 0;

  /// @brief Structured (no-YAML) initialization from an in-memory config.
  ///
  /// Same semantics and precedence as the YAML string path — config
  /// curves > runtime file > baked-in — but the caller supplies a
  /// GocartConfig directly, avoiding text parsing on the hot path.
  /// @param config Fully-populated package configuration.
  /// @note Default throws so packages without structured support
  /// (MAM4xx wrapper, test doubles) compile and behave unchanged.
  virtual void initialize(const GocartConfig &config) {
    (void)config;
    throw std::logic_error("EX-aero Error: this package does not support "
                           "structured configuration");
  }

  /// @brief Passive microphysics transport step (tendency / advection).
  ///
  /// Advances the prognostic aerosol mass state by @p delta_time_sec
  /// under the supplied environment. Mass-conserving: no thermodynamic
  /// state is mutated.
  /// @param env Read-only environmental boundary conditions.
  /// @param[in,out] state Aerosol state, (cell, level, species) [kg m^-3].
  /// @param delta_time_sec Integration timestep length [s].
  virtual void
  executeMicrophysics(const EnvironmentalStateView &env,
                      View3D<double> &state, // (cell, level, species)
                      double delta_time_sec) = 0;

  /// @brief Derived-diagnostics step (PM2.5, PM10, number, column mass, ALW).
  ///
  /// Fills the diagnostic slab indexed by exaero::diagnostic_indices.
  /// Wet sizing uses the RH-clamped growth factor; number concentration
  /// is derived from mass via the particle volume/mass relation.
  /// @param env Read-only environmental boundary conditions.
  /// @param state Aerosol mass state, (cell, level, species) [kg m^-3].
  /// @param[out] diagnostics_out (cell, level, diagnostic_index).
  virtual void computeDerivedDiagnostics(
      const EnvironmentalStateView &env,
      const View3D<const double> &state,    // (cell, level, species)
      View3D<double> &diagnostics_out) = 0; // (cell, level, diagnostic_index)

  /// @brief Optical-properties step (extinction, scattering, lidar, AOT).
  ///
  /// Evaluates the optical index set (exaero::optical_indices) for each
  /// requested band via RH-dependent, spectrally-interpolated MIE
  /// tables. Column AOT follows from layer integration of the
  /// extinction coefficient over @f$ \Delta z @f$ (Beer–Lambert).
  /// @param env Read-only environmental boundary conditions.
  /// @param state Aerosol mass state, (cell, level, species) [kg m^-3].
  /// @param wavelengths Queried free-space wavelengths [m], (band_index).
  /// @param[out] optics_out (cell, level, band_index, optics_index).
  virtual void
  computeOptics(const EnvironmentalStateView &env,
                const View3D<const double> &state, // (cell, level, species)
                const View1D<const double>
                    &wavelengths, // [m] queried wavelengths (band_index)
                View4D<double>
                    &optics_out) = 0; // (cell, level, band_index, optics_index)

  /// @brief Cloud-microphysics CCN activation spectra.
  ///
  /// Computes liquid cloud-condensation-number activation over the
  /// requested supersaturations (Köhler theory: dry number, wet growth,
  /// and critical supersaturation per species).
  /// @param env Read-only environmental boundary conditions.
  /// @param state Aerosol mass state, (cell, level, species) [kg m^-3].
  /// @param supersaturations Queried supersaturations [fraction], (S_index).
  /// @param[out] ccn_out (cell, level, S_index,
  ///             activated_ccn_number_concentration) [particles m^-3].
  virtual void computeCCN(
      const EnvironmentalStateView &env, const View3D<const double> &state,
      const View1D<const double> &
          supersaturations, // [fraction 0-1] queried supersaturations (S_index)
      View4D<double> &ccn_out) = 0; // (cell, level, S_index,
                                    // activated_ccn_number_concentration)

  /// @brief Map raw emissions to package-specific mass/number rates.
  ///
  /// Converts each raw (CECE) species flux into the package's own
  /// species via the configured split fractions and modal/number
  /// parameters. Flux units are given by EmissionsInputView::flux_type.
  /// @param env Read-only environmental boundary conditions (thickness).
  /// @param emissions_in Raw flux field plus its FluxType interpretation.
  /// @param[out] emissions_out (cell, level, target_species_index).
  virtual void computeEmissions(const EnvironmentalStateView &env,
                                const EmissionsInputView &emissions_in,
                                View3D<double> &emissions_out) = 0;

  /// @brief Resolve a species name to its packed state index.
  ///
  /// Dynamic name-to-index mapping so the host never hardcodes species
  /// ordering (which depends on the active configuration).
  /// @param name Species label (e.g. "DU", "SS").
  /// @return Zero-based index into the species axis, or -1 when absent.
  virtual int getSpeciesIndex(const std::string &name) const = 0;

  /// @brief Resolve a packed state index back to its species name.
  /// @param index Zero-based species index from getSpeciesIndex().
  /// @return Species label; undefined if @p index is out of range.
  virtual std::string getSpeciesName(int index) const = 0;

  // --- GEOSmie MIE table attribute surface (additive, backward compatible) ---
  // Defaults return NotActivated / do nothing so existing packages (MAM4xx
  // wrapper, tests) compile and behave unchanged.

  /// @brief Hot-path, device-resident bulk attribute query (optics-like).
  ///
  /// Fills a per-gridpoint attribute slab for one species and category,
  /// mirroring the computeOptics() execution model so it stays on-device.
  /// @param env Read-only environmental boundary conditions.
  /// @param state Aerosol mass state, (cell, level, species) [kg m^-3].
  /// @param species_index Target species (from getSpeciesIndex()).
  /// @param category Attribute category (see AttributeCategory).
  /// @param wavelengths Queried wavelengths [m]; ignored for Microphysical.
  /// @param[out] attributes_out (cell, level, attribute_index) for one
  ///             species + category.
  /// @param[out] status_out Optional; receives the AttributeStatus code per
  ///             (cell, level, attribute_index). May be null.
  virtual void computeAttributes(const EnvironmentalStateView &env,
                                 const View3D<const double> &state,
                                 int species_index, AttributeCategory category,
                                 const View1D<const double> &wavelengths,
                                 View3D<double> &attributes_out,
                                 View3D<int> *status_out = nullptr) {
    (void)env;
    (void)state;
    (void)species_index;
    (void)category;
    (void)wavelengths;
    (void)attributes_out;
    (void)status_out;
  }

  /// @brief Scalar attribute query for init-time / retrieval consumers.
  ///
  /// Returns the availability status; @p value_out and @p provenance_out
  /// are written ONLY when the status is available or interpolated —
  /// never a silent 0 on a not-activated / not-in-source result.
  /// @param species_index Target species (from getSpeciesIndex()).
  /// @param category Attribute category (see AttributeCategory).
  /// @param attribute_index Index code from the category's *_indices.
  /// @param rh Relative humidity [fraction 0-1] at the query point.
  /// @param wavelength_m Wavelength [m]; ignored (pass NaN) for Microphysical.
  /// @param[out] value_out Receives the surfaced value on an ok status.
  /// @param[out] provenance_out Optional unit/version/citation provenance.
  /// @return AttributeStatus describing availability and data source.
  virtual AttributeStatus
  queryAttribute(int species_index, AttributeCategory category,
                 int attribute_index, double rh, double wavelength_m,
                 double *value_out,
                 ProvenanceInfo *provenance_out = nullptr) const {
    (void)species_index;
    (void)category;
    (void)attribute_index;
    (void)rh;
    (void)wavelength_m;
    (void)value_out;
    (void)provenance_out;
    return AttributeStatus::NotActivated;
  }

  /// @brief Enable attribute categories/species with an optional override file.
  ///
  /// Controls which (species, category) pairs are surfaced. A runtime
  /// extension/override file that fails validation aborts initialization
  /// with a @c "FATAL ERROR:" diagnostic — no silent fallback.
  /// @param species Species labels to activate; empty => all species.
  /// @param categories_mask Bitmask of categories (see attribute_category_bit).
  /// @param runtime_file_path Optional runtime MIE file; empty => none.
  virtual void
  setAttributeActivation(const std::vector<std::string> &species,
                         int categories_mask,
                         const std::string &runtime_file_path = "") {
    (void)species;
    (void)categories_mask;
    (void)runtime_file_path;
  }

  /// @brief Bind config species to source MIE curves (CATChem-style).
  ///
  /// Declares how runtime-configurable species map onto source curves
  /// with optional per-attribute overrides on the config's own radius
  /// nodes. Resolved during initialize(); a bad config aborts with a
  /// @c "FATAL ERROR:" diagnostic.
  /// @param curves One SpeciesCurveConfig per runtime species.
  virtual void
  setSpeciesCurveConfig(const std::vector<SpeciesCurveConfig> &curves) {
    (void)curves;
  }

  /// @brief Polarized-moment (element, moment) slot counts for a species.
  ///
  /// Sizes the output buffer for a PolarizedMoment computeAttributes()
  /// call: the number of (element, moment) slots the species' curve
  /// carries, derived from data (never a fixed literal). Both outputs
  /// are 0 when the species ships no moments or the category is not
  /// activated.
  /// @param species_index Target species (from getSpeciesIndex()).
  /// @param[out] num_pol_out Polarized element count (0 if none).
  /// @param[out] num_moment_out Moment count per element (0 if none).
  /// @note Defaults to zeros so packages without moments compile unchanged.
  virtual void momentCounts(int species_index, int *num_pol_out,
                            int *num_moment_out) const {
    (void)species_index;
    if (num_pol_out)
      *num_pol_out = 0;
    if (num_moment_out)
      *num_moment_out = 0;
  }
};

} // namespace exaero
