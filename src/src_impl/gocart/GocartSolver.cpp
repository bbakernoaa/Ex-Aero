#include <gocart/GocartSpeciesParams.hpp>
#include <Kokkos_Core.hpp>

namespace exaero {

    // Implement public environment lifecycle routines (declared in exaero/Environment.hpp)
    void initialize_environment() {
        if (!Kokkos::is_initialized()) {
            Kokkos::initialize();
        }
    }

    void finalize_environment() {
        if (Kokkos::is_initialized()) {
            Kokkos::finalize();
        }
    }

    // Definition of our private solver state managing GPU-allocated memory views
    struct GocartSolverState {
        int num_species;
        Kokkos::View<GocartSpeciesParams*, Kokkos::DefaultExecutionSpace> d_species_params;
    };

    // Extern C/C++ helper functions declared in GocartPackage.cpp
    GocartSolverState* create_solver_state(int num_species, const GocartSpeciesParams* params) {
        auto* state = new GocartSolverState();
        state->num_species = num_species;

        // Allocate Device View
        state->d_species_params = Kokkos::View<GocartSpeciesParams*, Kokkos::DefaultExecutionSpace>(
            "d_species_params", num_species
        );

        // Allocate Host View Mirror
        auto h_view = Kokkos::create_mirror_view(state->d_species_params);
        for (int i = 0; i < num_species; ++i) {
            h_view(i) = params[i];
        }

        // Deep copy values to GPU Default Execution Space
        Kokkos::deep_copy(state->d_species_params, h_view);

        return state;
    }

    void free_solver_state(GocartSolverState* state) {
        if (state) {
            delete state;
        }
    }

    void run_gocart_diagnostics(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_species,
        const double* rh_ptr, const double* thick_ptr, const double* state_ptr, double* diags_ptr) {

        using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

        // Wrap incoming raw pointers into unmanaged Default Execution Space memory Views (zero copy!)
        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_diags = Kokkos::View<double***, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            diags_ptr, num_cells, num_levels, diagnostic_indices::NUM_DIAGNOSTICS
        );

        auto d_species_params = state->d_species_params;

        // Execute parallel diagnostics calculation on the GPU
        Kokkos::parallel_for("GocartDiagnosticsKernel", 
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_cells, num_levels}),
            KOKKOS_LAMBDA(int i_cell, int i_level) {
                double total_pm2_5_3d = 0.0;
                double total_pm10_3d = 0.0;
                double total_number_3d = 0.0;
                double total_sad_3d = 0.0;

                double current_rh = d_rh(i_cell, i_level);

                for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                    const auto& params = d_species_params(i_spec);
                    double mass_conc = d_state(i_cell, i_level, i_spec); // [kg/m³]
                    
                    if (mass_conc <= 0.0) continue;

                    // 1. Derived Number Concentration from Bulk Mass (lognormal population)
                    double ln_sig = Kokkos::log(params.lognormal_sigma);
                    double vol_factor = (M_PI / 6.0) * params.dry_density * 
                                        Kokkos::pow(params.lognormal_dg, 3) * 
                                        Kokkos::exp(4.5 * ln_sig * ln_sig);
                    double num_conc = mass_conc / vol_factor; // [particles/m³]
                    total_number_3d += num_conc;

                    // 2. Surface Area Density (SAD) [m²/m³]
                    double sad_factor = M_PI * Kokkos::pow(params.lognormal_dg, 2) * Kokkos::exp(2.0 * ln_sig * ln_sig);
                    total_sad_3d += num_conc * sad_factor;

                    // 3. 3D PM2.5 and PM10 size cuts
                    if (params.dry_particle_diameter <= 2.5e-6) {
                        total_pm2_5_3d += mass_conc;
                    }
                    if (params.dry_particle_diameter <= 10.0e-6) {
                        total_pm10_3d += mass_conc;
                    }
                }

                // Write 3D grid cell diagnostics
                d_diags(i_cell, i_level, diagnostic_indices::MASS_CONCENTRATION) = d_state(i_cell, i_level, 0); // principal species mass
                d_diags(i_cell, i_level, diagnostic_indices::PM2_5_CONCENTRATION) = total_pm2_5_3d;
                d_diags(i_cell, i_level, diagnostic_indices::PM10_CONCENTRATION) = total_pm10_3d;
                d_diags(i_cell, i_level, diagnostic_indices::NUMBER_CONCENTRATION) = total_number_3d;
                d_diags(i_cell, i_level, diagnostic_indices::SURFACE_AREA_DENSITY) = total_sad_3d;
            }
        );
        Kokkos::fence();

        // 2D Column Integration and Surface extractions (k=0 surface layer, integrated sum vertically)
        Kokkos::parallel_for("GocartColumnIntegratorKernel", 
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, num_cells),
            KOKKOS_LAMBDA(int i_cell) {
                double total_col_mass = 0.0;
                double total_col_pm25 = 0.0;

                for (int i_level = 0; i_level < num_levels; ++i_level) {
                    double dz = d_thick(i_cell, i_level);
                    
                    // Sum vertical column mass
                    for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                        double mass_conc = d_state(i_cell, i_level, i_spec);
                        const auto& params = d_species_params(i_spec);

                        if (mass_conc > 0.0) {
                            total_col_mass += mass_conc * dz;
                            if (params.dry_particle_diameter <= 2.5e-6) {
                                total_col_pm25 += mass_conc * dz;
                            }
                        }
                    }
                }

                // Surface layer is bottom layer (i_level = 0)
                double surface_mass = 0.0;
                double surface_pm25 = 0.0;
                for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                    double mass_conc = d_state(i_cell, 0, i_spec);
                    const auto& params = d_species_params(i_spec);
                    if (mass_conc > 0.0) {
                        surface_mass += mass_conc;
                        if (params.dry_particle_diameter <= 2.5e-6) {
                            surface_pm25 += mass_conc;
                        }
                    }
                }

                // Store 2D column-integrated values in the bottom-most level slots of output
                d_diags(i_cell, 0, diagnostic_indices::SURFACE_MASS) = surface_mass;
                d_diags(i_cell, 0, diagnostic_indices::COLUMN_MASS) = total_col_mass;
                d_diags(i_cell, 0, diagnostic_indices::SURFACE_PM2_5_MASS) = surface_pm25;
                d_diags(i_cell, 0, diagnostic_indices::COLUMN_PM2_5_MASS) = total_col_pm25;
            }
        );
        Kokkos::fence();
    }

    void run_gocart_optics(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_species,
        const double* rh_ptr, const double* thick_ptr, const double* state_ptr, double* optics_ptr) {
        
        // Passive placeholder for Task 3 & 4 (will be fully implemented in Task 6)
    }

} // namespace exaero
