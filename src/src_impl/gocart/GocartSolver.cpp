#include <gocart/GocartSpeciesParams.hpp>
#include <Kokkos_Core.hpp>

namespace exaero {

namespace {

    // Device-safe bracket search on a strictly-increasing coordinate axis with
    // clamp at both edges (never extrapolates beyond declared points, FR-007).
    // log_space != 0 blends in log coordinates (reference spectral interpolation,
    // research R3); otherwise linear. Mirrors the host MieTableStore locate rules.
    KOKKOS_INLINE_FUNCTION void curve_locate(const double* axis, int n, double v,
                                             int log_space, int& i0, int& i1, double& w) {
        if (n <= 1) { i0 = 0; i1 = 0; w = 0.0; return; }
        if (v <= axis[0]) { i0 = 0; i1 = 0; w = 0.0; return; }
        if (v >= axis[n - 1]) { i0 = n - 1; i1 = n - 1; w = 0.0; return; }
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            const int mid = (lo + hi) / 2;
            if (axis[mid] <= v) lo = mid; else hi = mid;
        }
        i0 = lo; i1 = hi;
        const double v0 = axis[i0], v1 = axis[i1];
        if (log_space && v0 > 0.0 && v1 > 0.0 && v > 0.0) {
            const double l0 = Kokkos::log(v0);
            const double span = Kokkos::log(v1) - l0;
            w = (span > 0.0) ? (Kokkos::log(v) - l0) / span : 0.0;
        } else {
            const double span = v1 - v0;
            w = (span > 0.0) ? (v - v0) / span : 0.0;
        }
    }

    // GEOSmie table spectral read (T028): mass extinction [m^2/kg], single-scattering
    // albedo and asymmetry at (rh, band coordinate) for the species' solver radius node.
    // Block layout documented in GocartSpeciesParams; the extents are data, never literals.
    // Returns false when the requested band coordinate lies OUTSIDE the curve's declared
    // lambda domain (e.g. a physical wavelength in metres against a band-index axis):
    // the caller must then fall back to its analytical path rather than silently clamping
    // to an unrelated band (FR-007 fail-loud intent, no silent wrong physics).
    KOKKOS_INLINE_FUNCTION bool species_spectral_table(const GocartSpeciesParams& params,
                                                        const double* pool, double rh,
                                                        double band, double& ext_per_mass,
                                                        double& ssa_v, double& g_v) {
        const double* sb = pool + params.spec_offset;
        const int nH = params.n_spec_rh;
        const int nL = params.n_spec_lambda;
        const int nR = params.n_spec_radius;
        const double* srh = sb;
        const double* slam = sb + nH;
        if (!(band >= slam[0] && band <= slam[nL - 1])) return false;
        const long nfield = static_cast<long>(nR) * nH * nL;
        const double* f_ext = sb + nH + nL;
        const double* f_ssa = f_ext + nfield;
        const double* f_g   = f_ssa + nfield;
        int rnode = params.solver_radius_node;
        if (rnode < 0) rnode = 0;
        if (rnode > nR - 1) rnode = nR - 1;
        int h0, h1, l0, l1; double wh, wl;
        curve_locate(srh, nH, rh, 0, h0, h1, wh);
        curve_locate(slam, nL, band, 1, l0, l1, wl);
        const long base = static_cast<long>(rnode) * nH;
        const double e00 = f_ext[(base + h0) * nL + l0];
        const double e10 = f_ext[(base + h1) * nL + l0];
        const double e01 = f_ext[(base + h0) * nL + l1];
        const double e11 = f_ext[(base + h1) * nL + l1];
        const double s00 = f_ssa[(base + h0) * nL + l0];
        const double s10 = f_ssa[(base + h1) * nL + l0];
        const double s01 = f_ssa[(base + h0) * nL + l1];
        const double s11 = f_ssa[(base + h1) * nL + l1];
        const double g00 = f_g[(base + h0) * nL + l0];
        const double g10 = f_g[(base + h1) * nL + l0];
        const double g01 = f_g[(base + h0) * nL + l1];
        const double g11 = f_g[(base + h1) * nL + l1];
        const double e_l = e00 + wh * (e10 - e00);
        const double e_h = e01 + wh * (e11 - e01);
        ext_per_mass = e_l + wl * (e_h - e_l);
        const double s_l = s00 + wh * (s10 - s00);
        const double s_h = s01 + wh * (s11 - s01);
        ssa_v = s_l + wl * (s_h - s_l);
        const double g_l = g00 + wh * (g10 - g00);
        const double g_h = g01 + wh * (g11 - g01);
        g_v = g_l + wl * (g_h - g_l);
        return true;
    }

} // anonymous namespace

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
        // ONE flat device double pool holding every species' curve block (ADR-003 R10).
        // Uploaded once here; never reallocated or H2D-copied inside the timestep loop.
        Kokkos::View<double*, Kokkos::DefaultExecutionSpace> d_curve_pool;
    };

    // Extern C/C++ helper functions declared in GocartPackage.cpp
    GocartSolverState* create_solver_state(int num_species, const GocartSpeciesParams* params,
                                           const double* pool, int pool_size) {
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

        // Single H2D upload of the flat curve pool (zero transfers in the loop, FR-015).
        if (pool_size > 0) {
            state->d_curve_pool = Kokkos::View<double*, Kokkos::DefaultExecutionSpace>(
                "d_curve_pool", pool_size);
            auto h_pool = Kokkos::create_mirror_view(state->d_curve_pool);
            for (int i = 0; i < pool_size; ++i) h_pool(i) = pool[i];
            Kokkos::deep_copy(state->d_curve_pool, h_pool);
        }

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
        auto d_curve_pool = state->d_curve_pool;

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

                    // 0. GEOSmie table-backed microphysics (T021): when this species is
                    //    curve-bound, growth factor and wet particle density are read from
                    //    the flat device pool (RH-linear, extent-driven scan — no literals).
                    double gf_table = -1.0;   // wet/dry effective-radius growth factor
                    double wd_table = -1.0;   // wet particle density [kg/m^3]
                    if (params.micro_offset >= 0) {
                        const double* mb = d_curve_pool.data() + params.micro_offset;
                        const int nH = params.n_micro_rh;
                        const int nR = params.n_micro_radius;
                        const double* mrh = mb;                       // [nH]
                        const double* mgf = mb + nH;                  // [nR*nH]
                        const double* mwd = mb + nH + nR * nH;        // [nR*nH]
                        int rnode = Kokkos::min(Kokkos::max(params.solver_radius_node, 0), nR - 1);
                        int i_bin = 0;
                        while (i_bin < nH - 2 && current_rh > mrh[i_bin + 1]) {
                            i_bin++;
                        }
                        double w = (current_rh - mrh[i_bin]) / (mrh[i_bin + 1] - mrh[i_bin] + 1e-15);
                        w = Kokkos::min(Kokkos::max(w, 0.0), 1.0); // clamp: never extrapolate (FR-007)
                        gf_table = mgf[rnode * nH + i_bin] + w * (mgf[rnode * nH + i_bin + 1] - mgf[rnode * nH + i_bin]);
                        wd_table = mwd[rnode * nH + i_bin] + w * (mwd[rnode * nH + i_bin + 1] - mwd[rnode * nH + i_bin]);
                    }

                    // 1. Derived Number Concentration from Bulk Mass (lognormal population)
                    double ln_sig = Kokkos::log(params.lognormal_sigma);
                    double vol_factor = (M_PI / 6.0) * params.dry_density * 
                                        Kokkos::pow(params.lognormal_dg, 3) * 
                                        Kokkos::exp(4.5 * ln_sig * ln_sig);
                    double num_conc = mass_conc / vol_factor; // [particles/m³]
                    total_number_3d += num_conc;

                    // 2. Wet size calculation: GEOSmie growth factor (table-backed, T021)
                    //    when available, else the kappa-Kohler approximation.
                    double growth = (gf_table > 0.0)
                        ? gf_table
                        : Kokkos::pow(1.0 + params.hygroscopicity * (current_rh / (1.0 - current_rh)), 1.0/3.0);
                    double wet_diameter = params.dry_particle_diameter * growth;

                    // 3. Surface Area Density (SAD) [m²/m³]
                    double sad_factor = M_PI * Kokkos::pow(params.lognormal_dg, 2) * Kokkos::exp(2.0 * ln_sig * ln_sig);
                    total_sad_3d += num_conc * sad_factor;

                    // 4. Aerosol Liquid Water (ALW) Content [kg/m³] (density of liquid water = 1000 kg/m³)
                    double dry_vol = (M_PI / 6.0) * Kokkos::pow(params.dry_particle_diameter, 3);
                    double wet_vol = (M_PI / 6.0) * Kokkos::pow(wet_diameter, 3);
                    double alw_mass_spec = num_conc * (wet_vol - dry_vol) * 1000.0; // [kg/m³]
                    total_alw_3d += alw_mass_spec;

                    // 5. Gravitational Settling Fall Velocity (vg) [m/s] (Stokes Settling)
                    // Wet density: GEOSmie table value (T021) when curve-bound, else the
                    // dry-volume + condensate-mass mixture approximation.
                    double dry_mass = dry_vol * params.dry_density;
                    double water_mass = (wet_vol - dry_vol) * 1000.0;
                    double wet_density = (wd_table > 0.0) ? wd_table : (dry_mass + water_mass) / wet_vol;

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
        auto d_curve_pool = state->d_curve_pool;

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

                    double tbl_ext = 0.0, tbl_ssa = 0.0, tbl_g = 0.0;
                    const bool table_hit =
                        params.spec_offset >= 0 &&
                        species_spectral_table(params, d_curve_pool.data(), current_rh,
                                               wavelength, tbl_ext, tbl_ssa, tbl_g);
                    if (table_hit) {
                        // --- MODE C: GEOSmie table spectral read (T028) ---
                        // Mass extinction per dry mass [m^2/kg] interpolated in RH (linear)
                        // and band coordinate (linear-in-log), at the solver radius node.
                        ext_coeff_spec = mass_conc * tbl_ext;
                        sca_coeff_spec = ext_coeff_spec * tbl_ssa;
                        asm_spec = tbl_g;
                    } else if (params.has_optics_lookup) {
                        // --- MODE B: RH Lookup Table 1D Linear Interpolation ---
                        // Curve read from the flat device pool; loop bound from n_rh extent,
                        // never a literal (ADR-003 R10). Block layout: {rh,ext,ssa,asm} x n_rh.
                        const double* blk = d_curve_pool.data() + params.curve_offset;
                        const int n_rh = params.n_rh;
                        const double* rh_axis = blk + curve_slot::RH_AXIS * n_rh;
                        const double* ext     = blk + curve_slot::EXT * n_rh;
                        const double* ssa     = blk + curve_slot::SSA * n_rh;
                        const double* asm_l   = blk + curve_slot::ASM * n_rh;
                        int i_bin = 0;
                        while (i_bin < n_rh - 2 && current_rh > rh_axis[i_bin + 1]) {
                            i_bin++;
                        }

                        double rh_lower = rh_axis[i_bin];
                        double rh_upper = rh_axis[i_bin + 1];
                        double weight = (current_rh - rh_lower) / (rh_upper - rh_lower + 1e-15);

                        // Linear Interpolation
                        double mee = ext[i_bin] + weight * (ext[i_bin + 1] - ext[i_bin]);
                        double ssa_v = ssa[i_bin] + weight * (ssa[i_bin + 1] - ssa[i_bin]);
                        asm_spec = asm_l[i_bin] + weight * (asm_l[i_bin + 1] - asm_l[i_bin]);

                        ext_coeff_spec = mass_conc * mee * 1000.0;
                        sca_coeff_spec = ext_coeff_spec * ssa_v;

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

                        double tbl_ext = 0.0, tbl_ssa = 0.0, tbl_g = 0.0;
                        const bool table_hit =
                            params.spec_offset >= 0 &&
                            species_spectral_table(params, d_curve_pool.data(), current_rh,
                                                   wavelength, tbl_ext, tbl_ssa, tbl_g);
                        if (table_hit) {
                            // --- MODE C: GEOSmie table spectral read (T028) ---
                            ext_coeff_spec = mass_conc * tbl_ext;
                            sca_coeff_spec = ext_coeff_spec * tbl_ssa;
                        } else if (params.has_optics_lookup) {
                            const double* blk = d_curve_pool.data() + params.curve_offset;
                            const int n_rh = params.n_rh;
                            const double* rh_axis = blk + curve_slot::RH_AXIS * n_rh;
                            const double* ext     = blk + curve_slot::EXT * n_rh;
                            const double* ssa     = blk + curve_slot::SSA * n_rh;
                            int i_bin = 0;
                            while (i_bin < n_rh - 2 && current_rh > rh_axis[i_bin + 1]) {
                                i_bin++;
                            }
                            double weight = (current_rh - rh_axis[i_bin]) / (rh_axis[i_bin + 1] - rh_axis[i_bin] + 1e-15);
                            double mee = ext[i_bin] + weight * (ext[i_bin + 1] - ext[i_bin]);
                            double ssa_v = ssa[i_bin] + weight * (ssa[i_bin + 1] - ssa[i_bin]);
                            ext_coeff_spec = mass_conc * mee * 1000.0;
                            sca_coeff_spec = ext_coeff_spec * ssa_v;
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

    void run_gocart_emissions(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_raw_species, int num_target_species,
        int flux_type_code,
        const double* thick_ptr,
        const double* raw_emissions_ptr,
        double* target_emissions_out_ptr) {

        using MemSpace = typename Kokkos::DefaultExecutionSpace::memory_space;

        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_raw = Kokkos::View<const double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            raw_emissions_ptr, num_cells, num_levels, num_raw_species
        );
        auto d_out = Kokkos::View<double***, Kokkos::LayoutLeft, MemSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            target_emissions_out_ptr, num_cells, num_levels, num_target_species
        );

        auto d_species_params = state->d_species_params;

        Kokkos::parallel_for("GocartEmissions_Kernel",
            Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {num_cells, num_levels, num_target_species}),
            KOKKOS_LAMBDA(int i_cell, int i_level, int i_spec) {
                
                const auto& params = d_species_params(i_spec);
                const auto& map = params.emissions_mapping;
                
                if (!map.is_active || map.raw_cece_index < 0 || map.raw_cece_index >= num_raw_species) {
                    d_out(i_cell, i_level, i_spec) = 0.0;
                    return;
                }

                // 1. Fetch raw flux and clamp negative or NaN values
                double raw_flux = d_raw(i_cell, i_level, map.raw_cece_index);
                if (Kokkos::isnan(raw_flux) || !Kokkos::isfinite(raw_flux) || raw_flux < 0.0) {
                    raw_flux = 0.0;
                }

                // 2. Perform Unit Scaling (Area Flux -> Volumetric)
                double vol_mass_flux = 0.0;
                if (flux_type_code == 1) { // AREA_FLUX
                    double dz = d_thick(i_cell, i_level);
                    if (dz > 1e-12) {
                        vol_mass_flux = raw_flux / dz;
                    } else {
                        vol_mass_flux = 0.0; // Avoid division-by-zero on micro layers
                    }
                } else { // MASS_CONCENTRATION_RATE
                    vol_mass_flux = raw_flux;
                }

                // 3. Apply mass split fraction
                double final_mass_rate = vol_mass_flux * map.mass_split_fraction;

                // 4. Handle Modal Number Conversion if applicable
                if (map.is_modal_mode && map.emitted_particle_diameter > 1e-12) {
                    double ln_sig = Kokkos::log(map.lognormal_sigma);
                    double vol_factor = (M_PI / 6.0) * params.dry_density * 
                                        Kokkos::pow(map.emitted_particle_diameter, 3) * 
                                        Kokkos::exp(4.5 * ln_sig * ln_sig);
                    
                    if (vol_factor > 1e-30) {
                        d_out(i_cell, i_level, i_spec) = final_mass_rate / vol_factor; // Emitted Number concentration
                    } else {
                        d_out(i_cell, i_level, i_spec) = 0.0;
                    }
                } else {
                    d_out(i_cell, i_level, i_spec) = final_mass_rate; // Standard mass rate
                }
            }
        );
        Kokkos::fence();
    }

} // namespace exaero
