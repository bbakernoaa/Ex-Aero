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

    // Public environment helpers
    void exaero_init_environment();
    void exaero_finalize_environment();

#ifdef __cplusplus
}
#endif

#endif // EXAERO_AEROSOL_PACKAGE_C_H
