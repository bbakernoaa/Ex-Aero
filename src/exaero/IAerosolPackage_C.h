#ifndef EXAERO_AEROSOL_PACKAGE_C_H
#define EXAERO_AEROSOL_PACKAGE_C_H

#ifdef __cplusplus
extern "C" {
#endif

    // Opaque handle representing the C++ class instance to Fortran
    typedef void* exaero_package_t;

    // Lifecycle helper to create and free GocartPackage
    exaero_package_t exaero_create_gocart_package();
    void exaero_free_package(exaero_package_t pkg);

    // Initializer with CCPP-standard error flags
    void exaero_initialize_package(exaero_package_t pkg, const char* config_yaml, char* errmsg, int* errflg);

    // Diagnostics calculation wrapper with CCPP-standard error flags
    void exaero_compute_diagnostics(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_species,
        const double* temp_ptr, const double* pres_ptr, const double* dens_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* diags_ptr,
        char* errmsg, int* errflg
    );

    // Emissions calculation wrapper with CCPP-standard error flags
    void exaero_compute_emissions(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_raw_species, int num_target_species,
        int flux_type_code,
        const double* thick_ptr,
        const double* raw_emissions_ptr,
        double* target_emissions_out_ptr,
        char* errmsg, int* errflg
    );

    // Optics calculation wrapper with CCPP-standard error flags
    void exaero_compute_optics(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_bands, int num_species,
        const double* wavelengths_ptr,
        const double* temp_ptr, const double* pres_ptr, const double* dens_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* optics_ptr,
        char* errmsg, int* errflg
    );

    // Cloud CCN Activation calculation wrapper with CCPP-standard error flags
    void exaero_compute_ccn(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_ss, int num_species,
        const double* ss_ptr,
        const double* temp_ptr, const double* rh_ptr, const double* state_ptr, double* ccn_ptr,
        char* errmsg, int* errflg
    );

    // Dynamic Species-to-Index mapping query wrappers (Hole 4)
    int exaero_get_species_index(exaero_package_t pkg, const char* name);
    void exaero_get_species_name(exaero_package_t pkg, int index, char* name_out, int max_len);

    // --- GEOSmie MIE attribute surface (contract §3, CCPP errmsg/errflg) ---
    // Activation control (FR-010) with optional runtime extension/override file (FR-012).
    // categories_mask is a bitmask of (1 << AttributeCategory). A bad file aborts with a
    // "FATAL ERROR:" diagnostic and errflg=1, no fallback (FR-009).
    void exaero_set_attribute_activation(
        exaero_package_t pkg,
        const char* const* species_names, int num_species,
        int categories_mask, const char* runtime_file_path,
        char* errmsg, int* errflg);

    // Curve mapping (R9): flat arrays + explicit counts, no fixed shapes across the boundary.
    // species_names/source_labels: length num_curves. radius_nodes[i]: length
    // num_radius_nodes[i] (0 => keep source axis). interpolate[i]: "linear"/"tension".
    // Overrides are delivered through the same call by pointing override_values[i] at a
    // flat (num_attributes * num_radius_nodes[i]) row-major buffer with matching
    // override_attribute_ids[i] (category<<8 | index) and num_override_attributes[i].
    void exaero_set_species_curve_config(
        exaero_package_t pkg,
        const char* const* species_names, const char* const* source_labels,
        const double* const* radius_nodes, const int* num_radius_nodes,
        const double* const* override_values, const int* num_override_attributes,
        const int* const* override_attribute_ids, const char* const* interpolate,
        int num_curves, char* errmsg, int* errflg);

    // Hot-path bulk query mirroring computeOptics (FR-015). attributes_out / status_out are
    // (num_cells, num_levels, num_attributes) layout_left; num_attributes is derived from
    // category by the caller via the public *_indices NUM_ATTRIBUTES constants.
    void exaero_compute_attributes(
        exaero_package_t pkg,
        int num_cells, int num_levels,
        const double* temp_ptr, const double* pres_ptr, const double* dens_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr,
        int species_index, int category,
        int num_bands, const double* wavelengths_ptr,
        double* attributes_out, int* status_out,
        char* errmsg, int* errflg);

    // Scalar query (FR-001..FR-003). wavelength_m may be NaN for Microphysical. unit_out /
    // version_out are caller-sized. status_out receives the AttributeStatus code.
    void exaero_query_attribute(
        exaero_package_t pkg,
        int species_index, int category, int attribute_index,
        double rh, double wavelength_m, double* value_out,
        char* unit_out, int unit_max, char* version_out, int version_max,
        int* status_out, char* errmsg, int* errflg);

    // Public environment helpers
    void exaero_init_environment();
    void exaero_finalize_environment();

#ifdef __cplusplus
}
#endif

#endif // EXAERO_AEROSOL_PACKAGE_C_H
