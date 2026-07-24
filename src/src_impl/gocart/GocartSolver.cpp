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
        int num_cells, int num_levels, int num_bands, int num_species,
        const double* wavelengths_ptr, const double* rh_ptr, const double* thick_ptr,
        const double* state_ptr, double* optics_ptr) {

        using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

        // Wrap raw pointers directly to unmanaged default space Kokkos views (zero copy!)
        auto d_wavelengths = Kokkos::View<const double*, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            wavelengths_ptr, num_bands
        );
        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_optics = Kokkos::View<double****, Kokkos::LayoutRight, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            optics_ptr, num_cells, num_levels, num_bands, optical_indices::NUM_OPTICS
        );

        auto d_species_params = state->d_species_params;

        // 1. Calculate 3D Optical Coefficients (Extinction, Scattering, Asymmetry) across all cells, levels, and bands simultaneously
        Kokkos::parallel_for("GocartOptics3D_Kernel", 
            Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {num_cells, num_levels, num_bands}),
            KOKKOS_LAMBDA(int i_cell, int i_level, int i_band) {
                double total_ext_coeff = 0.0;
                double total_sca_coeff = 0.0;
                double weighted_asymmetry = 0.0;

                double current_rh = Kokkos::min(d_rh(i_cell, i_level), 0.99); // defensively clamp relative humidity to prevent infinity
                double wavelength = d_wavelengths(i_band);                   // dynamic queried wavelength in meters

                for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                    const auto& params = d_species_params(i_spec);
                    double mass_conc = d_state(i_cell, i_level, i_spec);

                    if (mass_conc <= 0.0) continue;

                    double ext_coeff_spec = 0.0;
                    double sca_coeff_spec = 0.0;
                    double asm_spec = 0.0;

                    if (params.has_optics_lookup) {
                        // --- MODE B: RH Lookup Table 1D Linear Interpolation ---
                        // Find RH bin
                        int i_bin = 0;
                        while (i_bin < 6 && current_rh > params.rh_bins[i_bin+1]) {
                            i_bin++;
                        }
                        
                        double rh_lower = params.rh_bins[i_bin];
                        double rh_upper = params.rh_bins[i_bin+1];
                        double weight = (current_rh - rh_lower) / (rh_upper - rh_lower + 1e-15);

                        // Linear Interpolation
                        double mee = params.ext_lookup[i_bin] + weight * (params.ext_lookup[i_bin+1] - params.ext_lookup[i_bin]);
                        double ssa = params.ssa_lookup[i_bin] + weight * (params.ssa_lookup[i_bin+1] - params.ssa_lookup[i_bin]);
                        asm_spec = params.asm_lookup[i_bin] + weight * (params.asm_lookup[i_bin+1] - params.asm_lookup[i_bin]);

                        // Extinction and Scattering Coefficients (MEE is in m²/g, mass_conc is in kg/m³)
                        // Conversion factor: MEE [m²/g] * mass_conc [kg/m³] * 10^3 [g/kg] = [m⁻¹]
                        ext_coeff_spec = mass_conc * mee * 1000.0;
                        sca_coeff_spec = ext_coeff_spec * ssa;

                    } else {
                        // --- MODE A: Anomalous Diffraction Theory (ADT) Analytical Solver ---
                        // 1. Wet size growth
                        double wet_diameter = params.dry_particle_diameter * 
                                              Kokkos::pow(1.0 + params.hygroscopicity * (current_rh / (1.0 - current_rh)), 1.0/3.0);
                        
                        // 2. Size parameter x based on dynamic wavelength
                        double x = (M_PI * wet_diameter) / wavelength;

                        // 3. ADT Phase shift parameter rho
                        double rho = 2.0 * x * (params.refractive_index_real - 1.0);
                        rho = Kokkos::max(rho, 1e-12); // prevent division by zero

                        // 4. Extinction efficiency Q_ext (Anomalous Diffraction Theory)
                        double q_ext = 2.0 - (4.0 / rho) * Kokkos::sin(rho) + (4.0 / (rho * rho)) * (1.0 - Kokkos::cos(rho));

                        // 5. Scattering efficiency Q_sca (scaled based on imaginary index absorption)
                        double q_sca = q_ext * Kokkos::exp(-2.0 * x * params.refractive_index_imag);

                        // 6. Number Concentration
                        double ln_sig = Kokkos::log(params.lognormal_sigma);
                        double vol_factor = (M_PI / 6.0) * params.dry_density * 
                                            Kokkos::pow(params.lognormal_dg, 3) * 
                                            Kokkos::exp(4.5 * ln_sig * ln_sig);
                        double num_conc = mass_conc / vol_factor;

                        // 7. Coefficients calculations
                        double cross_section = (M_PI / 4.0) * wet_diameter * wet_diameter;
                        ext_coeff_spec = num_conc * cross_section * q_ext;
                        sca_coeff_spec = num_conc * cross_section * q_sca;
                        
                        // Parameterize asymmetry as a smooth function of wet size parameter
                        asm_spec = 0.7 * (x / (x + 1.0)); // standard asymptotic growth curve
                    }

                    total_ext_coeff += ext_coeff_spec;
                    total_sca_coeff += sca_coeff_spec;
                    weighted_asymmetry += sca_coeff_spec * asm_spec;
                }

                d_optics(i_cell, i_level, i_band, optical_indices::EXTINCTION_COEFF) = total_ext_coeff;
                d_optics(i_cell, i_level, i_band, optical_indices::SCATTERING_COEFF) = total_sca_coeff;
                d_optics(i_cell, i_level, i_band, optical_indices::BACKSCATTER_COEFF) = total_ext_coeff * 0.05; // standard backscatter ratio
                d_optics(i_cell, i_level, i_band, optical_indices::ASYMMETRY_FACTOR) = 
                    total_sca_coeff > 0.0 ? (weighted_asymmetry / total_sca_coeff) : 0.0;
            }
        );
        Kokkos::fence();

        // 2. Calculate 2D Column-Integrated Optical AOTs
        Kokkos::parallel_for("GocartOpticsAOT_Kernel", 
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_cells, num_bands}),
            KOKKOS_LAMBDA(int i_cell, int i_band) {
                double total_ext_aot = 0.0;
                double total_sca_aot = 0.0;
                double finemode_ext_aot = 0.0;
                double finemode_sca_aot = 0.0;
                double pm25_ext_aot = 0.0;
                double pm25_sca_aot = 0.0;

                double wavelength = d_wavelengths(i_band);

                for (int i_level = 0; i_level < num_levels; ++i_level) {
                    double dz = d_thick(i_cell, i_level);
                    double current_rh = Kokkos::min(d_rh(i_cell, i_level), 0.99);

                    for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                        double mass_conc = d_state(i_cell, i_level, i_spec);
                        if (mass_conc <= 0.0) continue;

                        const auto& params = d_species_params(i_spec);
                        double ext_coeff_spec = 0.0;
                        double sca_coeff_spec = 0.0;

                        if (params.has_optics_lookup) {
                            int i_bin = 0;
                            while (i_bin < 6 && current_rh > params.rh_bins[i_bin+1]) {
                                i_bin++;
                            }
                            double weight = (current_rh - params.rh_bins[i_bin]) / (params.rh_bins[i_bin+1] - params.rh_bins[i_bin] + 1e-15);
                            double mee = params.ext_lookup[i_bin] + weight * (params.ext_lookup[i_bin+1] - params.ext_lookup[i_bin]);
                            double ssa = params.ssa_lookup[i_bin] + weight * (params.ssa_lookup[i_bin+1] - params.ssa_lookup[i_bin]);
                            ext_coeff_spec = mass_conc * mee * 1000.0;
                            sca_coeff_spec = ext_coeff_spec * ssa;
                        } else {
                            // ADT values
                            double wet_diameter = params.dry_particle_diameter * 
                                                  Kokkos::pow(1.0 + params.hygroscopicity * (current_rh / (1.0 - current_rh)), 1.0/3.0);
                            double x = (M_PI * wet_diameter) / wavelength;
                            double rho = Kokkos::max(2.0 * x * (params.refractive_index_real - 1.0), 1e-12);
                            double q_ext = 2.0 - (4.0 / rho) * Kokkos::sin(rho) + (4.0 / (rho * rho)) * (1.0 - Kokkos::cos(rho));
                            double q_sca = q_ext * Kokkos::exp(-2.0 * x * params.refractive_index_imag);
                            double ln_sig = Kokkos::log(params.lognormal_sigma);
                            double vol_factor = (M_PI / 6.0) * params.dry_density * Kokkos::pow(params.lognormal_dg, 3) * Kokkos::exp(4.5 * ln_sig * ln_sig);
                            double num_conc = mass_conc / vol_factor;
                            double cross_section = (M_PI / 4.0) * wet_diameter * wet_diameter;
                            ext_coeff_spec = num_conc * cross_section * q_ext;
                            sca_coeff_spec = num_conc * cross_section * q_sca;
                        }

                        total_ext_aot += ext_coeff_spec * dz;
                        total_sca_aot += sca_coeff_spec * dz;

                        // Fine mode (sub-micron dry diameter)
                        if (params.dry_particle_diameter <= 1.0e-6) {
                            finemode_ext_aot += ext_coeff_spec * dz;
                            finemode_sca_aot += sca_coeff_spec * dz;
                        }

                        // PM2.5 (dry diameter <= 2.5 um)
                        if (params.dry_particle_diameter <= 2.5e-6) {
                            pm25_ext_aot += ext_coeff_spec * dz;
                            pm25_sca_aot += sca_coeff_spec * dz;
                        }
                    }
                }

                // Store 2D column AOT values in the bottom-most level slots of output
                d_optics(i_cell, 0, i_band, optical_indices::EXTINCTION_AOT) = total_ext_aot;
                d_optics(i_cell, 0, i_band, optical_indices::SCATTERING_AOT) = total_sca_aot;
                d_optics(i_cell, 0, i_band, optical_indices::FINE_MODE_EXTINCTION_AOT) = finemode_ext_aot;
                d_optics(i_cell, 0, i_band, optical_indices::FINE_MODE_SCATTERING_AOT) = finemode_sca_aot;
                d_optics(i_cell, 0, i_band, optical_indices::PM2_5_EXTINCTION_AOT) = pm25_ext_aot;
                d_optics(i_cell, 0, i_band, optical_indices::PM2_5_SCATTERING_AOT) = pm25_sca_aot;
            }
        );
        Kokkos::fence();
    }

} // namespace exaero
