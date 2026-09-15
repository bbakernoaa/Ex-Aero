#ifndef EXAERO_AEROSOL_PACKAGE_C_H
#define EXAERO_AEROSOL_PACKAGE_C_H
/// @file IAerosolPackage_C.h
/// @brief extern "C" surface of IAerosolPackage for Fortran/CCPP hosts.
///
/// Every routine follows the CCPP error convention: a caller-sized
/// @p errmsg buffer (null-terminated; empty string on success) plus an
/// @p errflg code (0 = success, 1 = fatal). Multidimensional arrays are
/// passed as flat pointers with explicit dimension counts in Fortran
/// column-major (layout_left) order — (cell, level, species) — so the
/// binding layer forwards host arrays zero-copy. No fixed-size shapes
/// cross this boundary; axis lengths are always data.

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Opaque handle to a C++ aerosol package instance (Fortran-safe).
///
/// Created by exaero_create_gocart_package(); never dereference on the
/// Fortran side — pass back to the exaero_* routines unchanged.
typedef void *exaero_package_t;

/// @brief Allocate a GocartPackage and return its opaque handle.
/// @return Handle for use in all other calls, or NULL on failure.
exaero_package_t exaero_create_gocart_package();

/// @brief Release a package handle (safe to call with NULL).
/// @param pkg Handle from exaero_create_gocart_package().
void exaero_free_package(exaero_package_t pkg);

/// @brief Initialize a package from a YAML configuration string.
/// @param pkg Opaque package handle.
/// @param config_yaml NUL-terminated YAML document text.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_initialize_package(exaero_package_t pkg, const char *config_yaml,
                               char *errmsg, int *errflg);

/// @brief Derived-diagnostics wrapper (PM2.5/PM10/number/column mass/ALW).
/// @param pkg Opaque package handle.
/// @param num_cells Horizontal cell count (fastest-varying Fortran dim).
/// @param num_levels Vertical level count.
/// @param num_species Species axis length of @p state_ptr.
/// @param temp_ptr [K] temperature, column-major (cell, level).
/// @param pres_ptr [Pa] pressure, column-major (cell, level).
/// @param dens_ptr [kg m^-3] air density, column-major (cell, level).
/// @param rh_ptr [fraction 0-1] relative humidity, column-major (cell, level).
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param state_ptr [kg m^-3] aerosol mass, column-major (cell, level,
/// species).
/// @param[out] diags_ptr Column-major (cell, level, NUM_DIAGNOSTICS); slot
///             ordering per exaero::diagnostic_indices.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_compute_diagnostics(exaero_package_t pkg, int num_cells,
                                int num_levels, int num_species,
                                const double *temp_ptr, const double *pres_ptr,
                                const double *dens_ptr, const double *rh_ptr,
                                const double *thick_ptr,
                                const double *state_ptr, double *diags_ptr,
                                char *errmsg, int *errflg);

/// @brief Emissions-mapping wrapper (raw CECE species -> package species).
/// @param pkg Opaque package handle.
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_raw_species Raw emissions species axis length.
/// @param num_target_species Package species axis length of the output.
/// @param flux_type_code 0 = mass-concentration rate [kg m^-3 s^-1],
///             1 = area flux [kg m^-2 s^-1] (see exaero::FluxType).
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param raw_emissions_ptr Column-major (cell, level, num_raw_species).
/// @param[out] target_emissions_out_ptr Column-major
///             (cell, level, num_target_species).
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_compute_emissions(exaero_package_t pkg, int num_cells,
                              int num_levels, int num_raw_species,
                              int num_target_species, int flux_type_code,
                              const double *thick_ptr,
                              const double *raw_emissions_ptr,
                              double *target_emissions_out_ptr, char *errmsg,
                              int *errflg);

/// @brief Optics wrapper (extinction/scattering/lidar/AOT across bands).
/// @param pkg Opaque package handle.
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_bands Wavelength band count.
/// @param num_species Species axis length of @p state_ptr.
/// @param wavelengths_ptr [m] queried free-space wavelengths, length num_bands.
/// @param temp_ptr [K] temperature, column-major (cell, level).
/// @param pres_ptr [Pa] pressure, column-major (cell, level).
/// @param dens_ptr [kg m^-3] air density, column-major (cell, level).
/// @param rh_ptr [fraction 0-1] relative humidity, column-major (cell, level).
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param state_ptr [kg m^-3] aerosol mass, column-major (cell, level,
/// species).
/// @param[out] optics_ptr Column-major (cell, level, num_bands, NUM_OPTICS);
///             slot ordering per exaero::optical_indices.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_compute_optics(exaero_package_t pkg, int num_cells, int num_levels,
                           int num_bands, int num_species,
                           const double *wavelengths_ptr,
                           const double *temp_ptr, const double *pres_ptr,
                           const double *dens_ptr, const double *rh_ptr,
                           const double *thick_ptr, const double *state_ptr,
                           double *optics_ptr, char *errmsg, int *errflg);

/// @brief Cloud CCN activation wrapper (Köhler liquid activation vs. S).
/// @param pkg Opaque package handle.
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param num_ss Supersaturation sample count.
/// @param num_species Species axis length of @p state_ptr.
/// @param ss_ptr [fraction] queried supersaturations, length num_ss.
/// @param temp_ptr [K] temperature, column-major (cell, level).
/// @param rh_ptr [fraction 0-1] relative humidity, column-major (cell, level).
/// @param state_ptr [kg m^-3] aerosol mass, column-major (cell, level,
/// species).
/// @param[out] ccn_ptr Column-major (cell, level, num_ss, 1) activated CCN
///             number concentration [particles m^-3].
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_compute_ccn(exaero_package_t pkg, int num_cells, int num_levels,
                        int num_ss, int num_species, const double *ss_ptr,
                        const double *temp_ptr, const double *rh_ptr,
                        const double *state_ptr, double *ccn_ptr, char *errmsg,
                        int *errflg);

/// @brief Resolve a species name to its packed state index.
/// @param pkg Opaque package handle.
/// @param name NUL-terminated species label (e.g. "DU").
/// @return Zero-based species index, or -1 on a null handle/name or unknown
///         species.
int exaero_get_species_index(exaero_package_t pkg, const char *name);

/// @brief Resolve a packed state index back to its species name.
/// @param pkg Opaque package handle.
/// @param index Zero-based species index.
/// @param[out] name_out Receives the NUL-terminated name (truncated to fit).
/// @param max_len Size of @p name_out in bytes.
void exaero_get_species_name(exaero_package_t pkg, int index, char *name_out,
                             int max_len);

// --- GEOSmie MIE attribute surface (CCPP errmsg/errflg) ---
/// @brief Activate attribute categories/species with an optional data file.
///
/// A runtime file that fails validation aborts with a @c "FATAL ERROR:"
/// diagnostic and errflg=1 — no silent fallback.
/// @param pkg Opaque package handle.
/// @param species_names Array of num_species NUL-terminated labels; NULL or
///             length 0 activates all species.
/// @param num_species Length of @p species_names.
/// @param categories_mask Bitmask of (1 << AttributeCategory).
/// @param runtime_file_path Optional runtime MIE file path; NULL/"" => none.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_set_attribute_activation(exaero_package_t pkg,
                                     const char *const *species_names,
                                     int num_species, int categories_mask,
                                     const char *runtime_file_path,
                                     char *errmsg, int *errflg);

/// @brief Install config-species curve bindings (jagged-array flat form).
///
/// Flat arrays + explicit counts only: no fixed shapes cross the boundary.
/// species_names/source_labels: length num_curves. radius_nodes[i]:
/// length num_radius_nodes[i] (0 => keep source axis). interpolate[i]:
/// "linear"/"tension". Overrides are delivered through the same call by
/// pointing override_values[i] at a flat (num_attributes *
/// num_radius_nodes[i]) row-major buffer with matching
/// override_attribute_ids[i] (category<<8 | index) and
/// num_override_attributes[i]. A bad config aborts with @c "FATAL ERROR:"
/// and errflg=1.
/// @param pkg Opaque package handle.
/// @param species_names Config-side species labels, length num_curves.
/// @param source_labels Source curve labels (e.g. "DU"), length num_curves.
/// @param radius_nodes Per-curve bin-centre arrays [m], length num_curves.
/// @param num_radius_nodes Length of each radius_nodes[i], length num_curves.
/// @param override_values Per-curve flat override rows, length num_curves.
/// @param num_override_attributes Override count per curve, length num_curves.
/// @param override_attribute_ids Per-curve (category<<8 | index) codes.
/// @param interpolate Per-curve strategy name ("linear"/"tension").
/// @param num_curves Number of curve bindings to install.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_set_species_curve_config(
    exaero_package_t pkg, const char *const *species_names,
    const char *const *source_labels, const double *const *radius_nodes,
    const int *num_radius_nodes, const double *const *override_values,
    const int *num_override_attributes,
    const int *const *override_attribute_ids, const char *const *interpolate,
    int num_curves, char *errmsg, int *errflg);

/// @brief Hot-path bulk attribute query mirroring exaero_compute_optics.
///
/// attributes_out / status_out are layout_left (num_cells, num_levels,
/// num_slots). num_slots = num_attributes for Microphysical (wavelength
/// ignored); for SpectralOptical it is num_bands * num_attributes in
/// band-major order (slot = band * num_attributes + attribute), mirroring
/// computeOptics. For PolarizedMoment the slots enumerate the addressable
/// (element, moment) pairs carried by THIS species' curve in dense
/// row-major order (element fastest), so num_slots = num_pol * num_moment
/// (data); call exaero_get_moment_counts to size the output, and a species
/// without moments yields 0 slots. num_attributes comes from the public
/// *_indices NUM_ATTRIBUTES constants (never a fixed literal).
/// @param pkg Opaque package handle.
/// @param num_cells Horizontal cell count.
/// @param num_levels Vertical level count.
/// @param temp_ptr [K] temperature, column-major (cell, level).
/// @param pres_ptr [Pa] pressure, column-major (cell, level).
/// @param dens_ptr [kg m^-3] air density, column-major (cell, level).
/// @param rh_ptr [fraction 0-1] relative humidity, column-major (cell, level).
/// @param thick_ptr [m] layer thickness, column-major (cell, level).
/// @param state_ptr [kg m^-3] the target species' mass, column-major
///             (cell, level, 1).
/// @param species_index Target species (from exaero_get_species_index()).
/// @param category AttributeCategory value as int.
/// @param num_bands Wavelength count (SpectralOptical only; else ignored).
/// @param wavelengths_ptr [m] queried wavelengths, length num_bands.
/// @param[out] attributes_out Column-major (cell, level, num_slots).
/// @param[out] status_out Optional AttributeStatus code per slot; NULL to skip.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_compute_attributes(
    exaero_package_t pkg, int num_cells, int num_levels, const double *temp_ptr,
    const double *pres_ptr, const double *dens_ptr, const double *rh_ptr,
    const double *thick_ptr, const double *state_ptr, int species_index,
    int category, int num_bands, const double *wavelengths_ptr,
    double *attributes_out, int *status_out, char *errmsg, int *errflg);

/// @brief Scalar attribute query for init-time / retrieval consumers.
///
/// value_out is written only on an available or interpolated result —
/// never a silent 0 for not-activated / not-in-source (FR-008).
/// @param pkg Opaque package handle.
/// @param species_index Target species (from exaero_get_species_index()).
/// @param category AttributeCategory value as int.
/// @param attribute_index Index code from the category's *_indices constants.
/// @param rh [fraction 0-1] relative humidity at the query point.
/// @param wavelength_m [m] query wavelength; may be NaN for Microphysical.
/// @param[out] value_out Receives the value on an ok status only.
/// @param[out] unit_out Canonical unit string; caller-sized (unit_max).
/// @param unit_max Size of @p unit_out in bytes.
/// @param[out] version_out Pinned table source version; caller-sized.
/// @param version_max Size of @p version_out in bytes.
/// @param[out] status_out Receives the AttributeStatus code.
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
void exaero_query_attribute(exaero_package_t pkg, int species_index,
                            int category, int attribute_index, double rh,
                            double wavelength_m, double *value_out,
                            char *unit_out, int unit_max, char *version_out,
                            int version_max, int *status_out, char *errmsg,
                            int *errflg);

/// @brief Query a species' polarized-moment (element, moment) slot counts.
///
/// Sized from the curve data (never a literal); outputs are 0 for a
/// species with no moments. Use to size a PolarizedMoment
/// exaero_compute_attributes() output buffer.
/// @param pkg Opaque package handle.
/// @param species_index Target species (from exaero_get_species_index()).
/// @param[out] num_pol_out Polarized element count (0 if none).
/// @param[out] num_moment_out Moment count per element (0 if none).
/// @param[out] errmsg Caller-sized error buffer; "" on success.
/// @param[out] errflg CCPP flag: 0 = success, 1 = fatal error.
/// @return 0 on success, 1 on error (bad handle / species index).
int exaero_get_moment_counts(exaero_package_t pkg, int species_index,
                             int *num_pol_out, int *num_moment_out,
                             char *errmsg, int *errflg);

/// @brief Initialize the Kokkos runtime (call once before any compute).
void exaero_init_environment();

/// @brief Finalize the Kokkos runtime (call once at teardown).
void exaero_finalize_environment();

#ifdef __cplusplus
}
#endif

#endif // EXAERO_AEROSOL_PACKAGE_C_H
