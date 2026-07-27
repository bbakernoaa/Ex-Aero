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

        // Wrap incoming raw pointers into unmanaged Default Execution Space column-major (LayoutLeft) memory Views (zero copy!)
        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_diags = Kokkos::View<double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
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
                double total_alw_3d = 0.0;
                double weighted_vg_numerator = 0.0;
                double total_mass_vg_weight = 0.0;

                double current_rh = Kokkos::min(d_rh(i_cell, i_level), 0.99);

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

                    // 2. Wet size calculation (Kohler hygroscopic growth approximation)
                    double wet_diameter = params.dry_particle_diameter * 
                                          Kokkos::pow(1.0 + params.hygroscopicity * (current_rh / (1.0 - current_rh)), 1.0/3.0);

                    // 3. Surface Area Density (SAD) [m²/m³]
                    double sad_factor = M_PI * Kokkos::pow(params.lognormal_dg, 2) * Kokkos::exp(2.0 * ln_sig * ln_sig);
                    total_sad_3d += num_conc * sad_factor;

                    // 4. Aerosol Liquid Water (ALW) Content [kg/m³] (density of liquid water = 1000 kg/m³)
                    double dry_vol = (M_PI / 6.0) * Kokkos::pow(params.dry_particle_diameter, 3);
                    double wet_vol = (M_PI / 6.0) * Kokkos::pow(wet_diameter, 3);
                    double alw_mass_spec = num_conc * (wet_vol - dry_vol) * 1000.0; // [kg/m³]
                    total_alw_3d += alw_mass_spec;

                    // 5. Gravitational Settling Fall Velocity (vg) [m/s] (Stokes Settling)
                    // Wet density calculation
                    double dry_mass = dry_vol * params.dry_density;
                    double water_mass = (wet_vol - dry_vol) * 1000.0;
                    double wet_density = (dry_mass + water_mass) / wet_vol;

                    double g_acc = 9.80665;       // [m/s²]
                    double dyn_visc = 1.825e-5;   // [kg/m-s] air dynamic viscosity
                    double vg_spec = (wet_density * wet_diameter * wet_diameter * g_acc) / (18.0 * dyn_visc);
                    
                    weighted_vg_numerator += mass_conc * vg_spec;
                    total_mass_vg_weight += mass_conc;

                    // 6. 3D PM2.5 and PM10 size cuts
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
                d_diags(i_cell, i_level, diagnostic_indices::AEROSOL_LIQUID_WATER) = total_alw_3d;
                d_diags(i_cell, i_level, diagnostic_indices::GRAVITATIONAL_SETTLING_VELOCITY) = 
                    total_mass_vg_weight > 0.0 ? (weighted_vg_numerator / total_mass_vg_weight) : 0.0;
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

        // Wrap raw pointers directly into unmanaged column-major (LayoutLeft) default space Kokkos views (zero copy!)
        auto d_wavelengths = Kokkos::View<const double*, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            wavelengths_ptr, num_bands
        );
        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_optics = Kokkos::View<double****, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            optics_ptr, num_cells, num_levels, num_bands, optical_indices::NUM_OPTICS
        );

        auto d_species_params = state->d_species_params;

        // 1. Calculate 3D Optical Coefficients (Extinction, Scattering, Lidar Backscatter, Asymmetry)
        Kokkos::parallel_for("GocartOptics3D_Kernel", 
            Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {num_cells, num_levels, num_bands}),
            KOKKOS_LAMBDA(int i_cell, int i_level, int i_band) {
                double total_ext_coeff = 0.0;
                double total_sca_coeff = 0.0;
                double weighted_asymmetry = 0.0;
                double total_lidar_backscatter = 0.0;

                double current_rh = Kokkos::min(d_rh(i_cell, i_level), 0.99); // clamp Relative Humidity
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
                        double q_ext = 0.0;
                        if (rho < 0.01) {
                            // --- stable Taylor series expansion to eliminate small-particle floating-point cancellation (Hole 3) ---
                            q_ext = 0.5 * rho * rho - (4.0 / 45.0) * Kokkos::pow(rho, 4) + (1.0 / 72.0) * Kokkos::pow(rho, 6);
                        } else {
                            // --- Standard ADT formula ---
                            q_ext = 2.0 - (4.0 / rho) * Kokkos::sin(rho) + (4.0 / (rho * rho)) * (1.0 - Kokkos::cos(rho));
                        }

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

                    // Henyey-Greenstein (HG) backscatter phase function evaluated at 180 degrees (theta = pi, cos_theta = -1):
                    // P_HG(pi, g) = (1 - g^2) / (1 + g^2 + 2g)^1.5 = (1 - g) / (1 + g)^2
                    double phase_hg_backscatter = (1.0 - asm_spec) / (4.0 * M_PI * Kokkos::pow(1.0 + asm_spec, 2));
                    double backscatter_spec = sca_coeff_spec * phase_hg_backscatter; // [m⁻¹ sr⁻¹]

                    total_ext_coeff += ext_coeff_spec;
                    total_sca_coeff += sca_coeff_spec;
                    weighted_asymmetry += sca_coeff_spec * asm_spec;
                    total_lidar_backscatter += backscatter_spec;
                }

                d_optics(i_cell, i_level, i_band, optical_indices::EXTINCTION_COEFF) = total_ext_coeff;
                d_optics(i_cell, i_level, i_band, optical_indices::SCATTERING_COEFF) = total_sca_coeff;
                d_optics(i_cell, i_level, i_band, optical_indices::BACKSCATTER_COEFF) = total_ext_coeff * 0.05; // standard backscatter ratio
                d_optics(i_cell, i_level, i_band, optical_indices::ASYMMETRY_FACTOR) = 
                    total_sca_coeff > 0.0 ? (weighted_asymmetry / total_sca_coeff) : 0.0;
                d_optics(i_cell, i_level, i_band, optical_indices::LIDAR_BACKSCATTER) = total_lidar_backscatter;
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
                            
                            double q_ext = 0.0;
                            if (rho < 0.01) {
                                // --- stable Taylor series expansion to eliminate small-particle floating-point cancellation (Hole 3) ---
                                q_ext = 0.5 * rho * rho - (4.0 / 45.0) * Kokkos::pow(rho, 4) + (1.0 / 72.0) * Kokkos::pow(rho, 6);
                            } else {
                                // --- Standard ADT formula ---
                                q_ext = 2.0 - (4.0 / rho) * Kokkos::sin(rho) + (4.0 / (rho * rho)) * (1.0 - Kokkos::cos(rho));
                            }

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

                        // Fine mode
                        if (params.dry_particle_diameter <= 1.0e-6) {
                            finemode_ext_aot += ext_coeff_spec * dz;
                            finemode_sca_aot += sca_coeff_spec * dz;
                        }

                        // PM2.5
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

    void run_gocart_ccn(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_ss, int num_species,
        const double* ss_ptr, const double* temp_ptr, const double* rh_ptr,
        const double* state_ptr, double* ccn_ptr) {

        using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

        auto d_ss = Kokkos::View<const double*, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            ss_ptr, num_ss
        );
        auto d_temp = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            temp_ptr, num_cells, num_levels
        );
        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_ccn = Kokkos::View<double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            ccn_ptr, num_cells, num_levels, num_ss
        );

        auto d_species_params = state->d_species_params;

        // Parallel Cloud CCN Activation Spectrum Solver on the GPU
        Kokkos::parallel_for("GocartCcnSpectrumKernel", 
            Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {num_cells, num_levels, num_ss}),
            KOKKOS_LAMBDA(int i_cell, int i_level, int i_ss) {
                double total_activated_ccn = 0.0;
                
                double temp = d_temp(i_cell, i_level);
                double query_ss = d_ss(i_ss); // supersaturation (fraction, e.g. 0.005 for 0.5% SS)

                // Kelvin parameter: A_kelvin = 2 * sigma * M_w / (rho_water * R * T) ≈ 1.2e-9 / T [m]
                double a_kelvin = 1.2e-9 / temp;

                for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                    const auto& params = d_species_params(i_spec);
                    double mass_conc = d_state(i_cell, i_level, i_spec);

                    if (mass_conc <= 0.0) continue;

                    // Calculate analytical critical supersaturation S_c for liquid cloud droplet activation
                    // S_c = sqrt( (4 * A^3) / (27 * kappa * D_dry^3) ) [fraction]
                    double s_crit = Kokkos::sqrt((4.0 * Kokkos::pow(a_kelvin, 3)) / 
                                                 (27.0 * params.hygroscopicity * Kokkos::pow(params.dry_particle_diameter, 3) + 1e-30));

                    // If queried supersaturation exceeds the critical activation threshold, the species activates!
                    if (query_ss >= s_crit) {
                        double ln_sig = Kokkos::log(params.lognormal_sigma);
                        double vol_factor = (M_PI / 6.0) * params.dry_density * 
                                            Kokkos::pow(params.lognormal_dg, 3) * 
                                            Kokkos::exp(4.5 * ln_sig * ln_sig);
                        double num_conc = mass_conc / vol_factor; // [particles/m³]
                        total_activated_ccn += num_conc;
                    }
                }

                d_ccn(i_cell, i_level, i_ss) = total_activated_ccn;
            }
        );
        Kokkos::fence();
    }

} // namespace exaero
