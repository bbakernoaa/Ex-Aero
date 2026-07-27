#include <exaero/IAerosolPackage_C.h>
#include <gocart/GocartPackage.hpp>
#include <exaero/Environment.hpp>
#include <exaero/AerosolIndices.hpp>

extern "C" {

    exaero_package_t exaero_create_gocart_package() {
        return static_cast<exaero_package_t>(new exaero::GocartPackage());
    }

    void exaero_free_package(exaero_package_t pkg) {
        if (pkg) {
            delete static_cast<exaero::GocartPackage*>(pkg);
        }
    }

    void exaero_initialize_package(exaero_package_t pkg, const char* config_yaml) {
        if (pkg && config_yaml) {
            static_cast<exaero::GocartPackage*>(pkg)->initialize(std::string(config_yaml));
        }
    }

    void exaero_compute_diagnostics(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_species,
        const double* temp_ptr, const double* pres_ptr, const double* dens_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* diags_ptr) {

        if (!pkg) return;

        auto* package = static_cast<exaero::GocartPackage*>(pkg);

        // Map incoming flat Fortran/C raw pointers to standard layout C++20 mdspan Views on-the-fly (LayoutLeft)
        exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
        exaero::View2D<const double> pressure(pres_ptr, num_cells, num_levels);
        exaero::View2D<const double> air_density(dens_ptr, num_cells, num_levels);
        exaero::View2D<const double> relative_humidity(rh_ptr, num_cells, num_levels);
        exaero::View2D<const double> layer_thickness(thick_ptr, num_cells, num_levels);

        exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

        exaero::View3D<const double> state(state_ptr, num_cells, num_levels, num_species);
        exaero::View3D<double> diagnostics_out(diags_ptr, num_cells, num_levels, exaero::diagnostic_indices::NUM_DIAGNOSTICS);

        // Execute dynamic diagnostics solver
        package->computeDerivedDiagnostics(env, state, diagnostics_out);
    }

    void exaero_compute_optics(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_bands, int num_species,
        const double* wavelengths_ptr,
        const double* temp_ptr, const double* pres_ptr, const double* dens_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* optics_ptr) {

        if (!pkg) return;

        auto* package = static_cast<exaero::GocartPackage*>(pkg);

        exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
        exaero::View2D<const double> pressure(pres_ptr, num_cells, num_levels);
        exaero::View2D<const double> air_density(dens_ptr, num_cells, num_levels);
        exaero::View2D<const double> relative_humidity(rh_ptr, num_cells, num_levels);
        exaero::View2D<const double> layer_thickness(thick_ptr, num_cells, num_levels);

        exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

        exaero::View1D<const double> wavelengths(wavelengths_ptr, num_bands);
        exaero::View3D<const double> state(state_ptr, num_cells, num_levels, num_species);
        exaero::View4D<double> optics_out(optics_ptr, num_cells, num_levels, num_bands, exaero::optical_indices::NUM_OPTICS);

        // Execute dynamic optics solver
        package->computeOptics(env, state, wavelengths, optics_out);
    }

    void exaero_compute_ccn(
        exaero_package_t pkg,
        int num_cells, int num_levels, int num_ss, int num_species,
        const double* ss_ptr,
        const double* temp_ptr, const double* pres_ptr, const double* dens_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* ccn_ptr) {

        if (!pkg) return;

        auto* package = static_cast<exaero::GocartPackage*>(pkg);

        exaero::View2D<const double> temperature(temp_ptr, num_cells, num_levels);
        exaero::View2D<const double> pressure(pres_ptr, num_cells, num_levels);
        exaero::View2D<const double> air_density(dens_ptr, num_cells, num_levels);
        exaero::View2D<const double> relative_humidity(rh_ptr, num_cells, num_levels);
        exaero::View2D<const double> layer_thickness(thick_ptr, num_cells, num_levels);

        exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

        exaero::View1D<const double> supersaturations(ss_ptr, num_ss);
        exaero::View3D<const double> state(state_ptr, num_cells, num_levels, num_species);
        exaero::View4D<double> ccn_out(ccn_ptr, num_cells, num_levels, num_ss, 1);

        // Execute dynamic cloud activation solver
        package->computeCCN(env, state, supersaturations, ccn_out);
    }

    // --- Dynamic Species-to-Index Queries (Hole 4) ---
    int exaero_get_species_index(exaero_package_t pkg, const char* name) {
        if (!pkg || !name) return -1;
        return static_cast<exaero::GocartPackage*>(pkg)->getSpeciesIndex(std::string(name));
    }

    void exaero_get_species_name(exaero_package_t pkg, int index, char* name_out, int max_len) {
        if (!pkg || !name_out || max_len <= 0) return;
        try {
            std::string name = static_cast<exaero::GocartPackage*>(pkg)->getSpeciesName(index);
            int len = static_cast<int>(name.length());
            if (len >= max_len) {
                len = max_len - 1;
            }
            for (int i = 0; i < len; ++i) {
                name_out[i] = name[i];
            }
            name_out[len] = '\0'; // ensure null-termination
        } catch (...) {
            name_out[0] = '\0';
        }
    }

    void exaero_init_environment() {
        exaero::initialize_environment();
    }

    void exaero_finalize_environment() {
        exaero::finalize_environment();
    }

} // extern "C"
