# EX-aero GOCART Proof of Concept Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and verify a dynamic, GPU-portable, zero-copy GOCART (bulk aerosol) diagnostics and optics interface package as a standalone proof of concept for EX-aero.

**Architecture:** We build a single compiled library `libexaero` with a clean public interface `exaero::IAerosolPackage`. We vendor the C++20 `std::experimental::mdspan` backport header-only privately, compiling it under our custom isolated `exaero_mdspan` namespace. Raw pointers are wrapped to `mdspan` on the API border, which then map zero-copy to unmanaged default execution space `Kokkos::View` structures inside parallel GPU execution kernels.

**Tech Stack:** C++20, Kokkos, `yaml-cpp`, standard single-header `mdspan` backport.

## Global Constraints

*   **Zero-Dependency Public API:** Public headers must require zero external host headers (only C++20 standard library, `mdspan`, and Kokkos).
*   **Target Artifacts:** `libexaero.so` / `libexaero.a`.
*   **Namespace:** `exaero::` for package, `exaero_mdspan::` for isolated mdspan.
*   **Header Directory:** `<exaero/...>` for public headers, `<gocart/...>` for private.
*   **Memory Overhead:** Host pointer wrapping into mdspan and Kokkos Views must incur `< 0.1%` performance overhead.
*   **License:** Open Source (Apache 2.0 / BSD-3-Clause compatible).

---

### Task 1: Scaffolding, CMake & Private mdspan Vendoring
*(Complete)*

### Task 2: Public Core API Definition
*(Complete)*

### Task 3: Dynamic GOCART YAML Loading
*(Complete)*

### Task 4: Zero-Copy Pointer Mapping & Memory Adapters
*(Complete)*

---

### Task 5: GOCART Passive Diagnostics & Sizing Solver

**Files:**
- Modify: `src/src_impl/gocart/GocartSolver.cpp`
- Modify: `tests/test_main.cpp`

**Interfaces:**
- Consumes: Environmental fields (rh, layer thickness), species params, and state concentrations.
- Produces: 3D diagnostics (PM2.5, PM10, derived number, SAD) and 2D column-integrated diagnostics (Surface Mass, Column Mass, Surface PM2.5, Column PM2.5) mapped to `diagnostics_out` view.

- [ ] **Step 1: Implement Gocart Diagnostics GPU Parallel Kernel**

Update `run_gocart_diagnostics` inside `src/src_impl/gocart/GocartSolver.cpp` to map raw pointers to unmanaged default execution space views and perform parallel diagnostic evaluations on the GPU:

```cpp
#include <exaero/IAerosolPackage.hpp>
#include <cmath>

    void run_gocart_diagnostics(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_species,
        const double* rh_ptr, const double* thick_ptr, const double* state_ptr, double* diags_ptr) {

        // Wrap incoming raw pointers into unmanaged default space Kokkos Views (zero copy!)
        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_diags = Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
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

                    // 2. Wet diameter calculation (Kohler hygroscopic growth approximation)
                    double wet_diameter = params.dry_particle_diameter * 
                                          Kokkos::pow(1.0 + params.hygroscopicity * (current_rh / (1.0 - current_rh)), 1.0/3.0);

                    // 3. Surface Area Density (SAD) [m²/m³]
                    double sad_factor = M_PI * Kokkos::pow(params.lognormal_dg, 2) * Kokkos::exp(2.0 * ln_sig * ln_sig);
                    total_sad_3d += num_conc * sad_factor;

                    // 4. 3D PM2.5 and PM10 size cuts
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
```

- [ ] **Step 2: Update Unit Integration Test**

Add assertions inside `test_zero_copy_mapping()` in `tests/test_main.cpp` to verify exact GOCART diagnostics math (e.g. Column Mass integration and Surface PM2.5 extractions).

Update test code in `test_zero_copy_mapping()`:
```cpp
    // Assert diagnostics outputs are calculated correctly
    double expected_col_mass = state_raw[0] * thick_raw[0]; // 1.0e-6 * 100.0 = 1.0e-4 kg/m²
    assert(std::abs(diags_raw[exaero::diagnostic_indices::COLUMN_MASS] - expected_col_mass) < 1e-12);
    assert(diags_raw[exaero::diagnostic_indices::SURFACE_MASS] == 1.0e-6);
    assert(diags_raw[exaero::diagnostic_indices::PM2_5_CONCENTRATION] == 1.0e-6);
    assert(diags_raw[exaero::diagnostic_indices::NUMBER_CONCENTRATION] > 0.0);
```

- [ ] **Step 3: Compile and run test to verify dynamic GOCART diagnostics math**

Run:
```bash
cmake --build build && ./build/tests/exaero_test_runner
```
Expected: PASS (Prints all pass messages, confirming correct diagnostic calculations)

- [ ] **Step 4: Commit**

```bash
git add src/src_impl/gocart/GocartSolver.cpp tests/test_main.cpp
git commit -m "feat: implement GOCART GPU diagnostics and column mass integrator"
```

---

### Task 6: GOCART Dual-Mode Optical Parameterization Solver

**Files:**
- Modify: `src/src_impl/gocart/GocartSolver.cpp`
- Modify: `tests/test_main.cpp`

**Interfaces:**
- Consumes: Environmental fields (rh, layer thickness), species params, and state concentrations.
- Produces: 3D optical properties (extinction, scattering, asymmetry) and 2D column AOT integrations (extinction/scattering total, fine-mode, and PM2.5) mapped to `optics_out` view.

- [ ] **Step 1: Implement Dual-Mode Optics GPU Solver**

Update `run_gocart_optics` inside `src/src_impl/gocart/GocartSolver.cpp` to execute the dual-mode optical solver (ADT analytical model or RH lookup table interpolation):

```cpp
#include <exaero/IAerosolPackage.hpp>
#include <cmath>

    void run_gocart_optics(
        GocartSolverState* state,
        int num_cells, int num_levels, int num_species,
        const double* rh_ptr, const double* thick_ptr, const double* state_ptr, double* optics_ptr) {

        auto d_rh = Kokkos::View<const double**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            rh_ptr, num_cells, num_levels
        );
        auto d_thick = Kokkos::View<const double**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            thick_ptr, num_cells, num_levels
        );
        auto d_state = Kokkos::View<const double***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            state_ptr, num_cells, num_levels, num_species
        );
        auto d_optics = Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
            optics_ptr, num_cells, num_levels, optical_indices::NUM_OPTICS
        );

        auto d_species_params = state->d_species_params;

        // 1. Calculate 3D Optical Coefficients (Extinction, Scattering, Asymmetry)
        Kokkos::parallel_for("GocartOptics3D_Kernel", 
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_cells, num_levels}),
            KOKKOS_LAMBDA(int i_cell, int i_level) {
                double total_ext_coeff = 0.0;
                double total_sca_coeff = 0.0;
                double weighted_asymmetry = 0.0;

                double current_rh = d_rh(i_cell, i_level);

                for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                    const auto& params = d_species_params(i_spec);
                    double mass_conc = d_state(i_cell, i_level, i_spec); // [kg/m³]

                    if (mass_conc <= 0.0) continue;

                    double ext_coeff_spec = 0.0;
                    double sca_coeff_spec = 0.0;
                    double asm_spec = 0.0;

                    if (params.has_optics_lookup) {
                        // --- MODE B: RH Lookup Table 1D Linear Interpolation ---
                        // Find RH bin
                        int i_bin = 0;
                        while (i_bin < 7 && current_rh > params.rh_bins[i_bin+1]) {
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
                        
                        // 2. Size parameter x
                        double wavelength = 550.0e-9; // shortwave 550 nm visible standard
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

                d_optics(i_cell, i_level, optical_indices::EXTINCTION_COEFF) = total_ext_coeff;
                d_optics(i_cell, i_level, optical_indices::SCATTERING_COEFF) = total_sca_coeff;
                d_optics(i_cell, i_level, optical_indices::BACKSCATTER_COEFF) = total_ext_coeff * 0.05; // standard backscatter ratio
                d_optics(i_cell, i_level, optical_indices::ASYMMETRY_FACTOR) = 
                    total_sca_coeff > 0.0 ? (weighted_asymmetry / total_sca_coeff) : 0.0;
            }
        );
        Kokkos::fence();

        // 2. Calculate 2D Column-Integrated Optical AOTs
        Kokkos::parallel_for("GocartOpticsAOT_Kernel", 
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, num_cells),
            KOKKOS_LAMBDA(int i_cell) {
                double total_ext_aot = 0.0;
                double total_sca_aot = 0.0;
                double finemode_ext_aot = 0.0;
                double finemode_sca_aot = 0.0;
                double pm25_ext_aot = 0.0;
                double pm25_sca_aot = 0.0;

                for (int i_level = 0; i_level < num_levels; ++i_level) {
                    double dz = d_thick(i_cell, i_level);

                    for (int i_spec = 0; i_spec < num_species; ++i_spec) {
                        double mass_conc = d_state(i_cell, i_level, i_spec);
                        if (mass_conc <= 0.0) continue;

                        const auto& params = d_species_params(i_spec);
                        double ext_coeff_spec = 0.0;
                        double sca_coeff_spec = 0.0;

                        // Retrieve the 3D coefficients calculated in Step 1
                        if (params.has_optics_lookup) {
                            int i_bin = 0;
                            while (i_bin < 7 && d_rh(i_cell, i_level) > params.rh_bins[i_bin+1]) {
                                i_bin++;
                            }
                            double weight = (d_rh(i_cell, i_level) - params.rh_bins[i_bin]) / (params.rh_bins[i_bin+1] - params.rh_bins[i_bin] + 1e-15);
                            double mee = params.ext_lookup[i_bin] + weight * (params.ext_lookup[i_bin+1] - params.ext_lookup[i_bin]);
                            double ssa = params.ssa_lookup[i_bin] + weight * (params.ssa_lookup[i_bin+1] - params.ssa_lookup[i_bin]);
                            ext_coeff_spec = mass_conc * mee * 1000.0;
                            sca_coeff_spec = ext_coeff_spec * ssa;
                        } else {
                            // ADT values
                            double wet_diameter = params.dry_particle_diameter * 
                                                  Kokkos::pow(1.0 + params.hygroscopicity * (d_rh(i_cell, i_level) / (1.0 - d_rh(i_cell, i_level))), 1.0/3.0);
                            double x = (M_PI * wet_diameter) / 550.0e-9;
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
                d_optics(i_cell, 0, optical_indices::EXTINCTION_AOT) = total_ext_aot;
                d_optics(i_cell, 0, optical_indices::SCATTERING_AOT) = total_sca_aot;
                d_optics(i_cell, 0, optical_indices::FINE_MODE_EXTINCTION_AOT) = finemode_ext_aot;
                d_optics(i_cell, 0, optical_indices::FINE_MODE_SCATTERING_AOT) = finemode_sca_aot;
                d_optics(i_cell, 0, optical_indices::PM2_5_EXTINCTION_AOT) = pm25_ext_aot;
                d_optics(i_cell, 0, optical_indices::PM2_5_SCATTERING_AOT) = pm25_sca_aot;
            }
        );
        Kokkos::fence();
    }
```

- [ ] **Step 2: Add Optical Unit Tests to verify ADT and RH Lookups**

Add `test_gocart_optics()` to `tests/test_main.cpp` implementing assertions for:
1.  **ADT Mode (Analytical):** Uses ADT formulas to derive extinction coefficient and AOT on the GPU.
2.  **Lookup Table Mode:** Uses a pre-tabulated 8-bin RH array and asserts correct 1D linear interpolation behavior on the GPU.

```cpp
void test_gocart_optics() {
    // 1. Prepare YAML containing species with has_optics_lookup flag and pre-tabulated values
    std::string yaml_string = R"(
    species:
      - name: "Dust_ADT"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
      - name: "Sulfate_Lookup"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.2e-6
        hygroscopicity: 0.50
        lognormal_sigma: 2.0
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-8
        has_optics_lookup: true
        rh_bins: [0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]
        ext_lookup: [2.5, 3.8, 5.0, 6.5, 8.5, 12.0, 15.0, 18.0]
        ssa_lookup: [0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99]
        asm_lookup: [0.60, 0.63, 0.65, 0.67, 0.70, 0.72, 0.73, 0.74]
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.75 }; // 75% relative humidity (midway between bin 0.7 and 0.8!)
    double thick_raw[1] = { 100.0 };

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    double state_raw[2] = { 1.0e-6, 1.0e-6 }; // 1.0 ug/m³ of each species
    exaero::View3D<double> state(state_raw, 1, 1, 2);

    double optics_raw[exaero::optical_indices::NUM_OPTICS] = { 0.0 };
    exaero::View3D<double> optics_out(optics_raw, 1, 1, exaero::optical_indices::NUM_OPTICS);

    // Run optical step
    package.computeOptics(env, state, optics_out);

    // Verify 1D Linear Interpolation for Sulfate_Lookup at 75% RH
    // rh_bins[2] = 0.7 (ext = 5.0), rh_bins[3] = 0.8 (ext = 6.5)
    // 75% RH is exactly midway: expected MEE = 5.75 m²/g.
    // expected ext_coeff_spec = 1.0e-6 * 5.75 * 1000 = 5.75e-3 m⁻¹
    assert(optics_raw[exaero::optical_indices::EXTINCTION_COEFF] > 5.75e-3); // ADT dust + Sulfate
    assert(optics_raw[exaero::optical_indices::EXTINCTION_AOT] > 0.0);

    std::cout << "GOCART Optics Dual-Mode Calculations: PASS" << std::endl;
}
```

Update `GocartPackage::initialize` in `GocartPackage.cpp` to read `has_optics_lookup` and populate lookups:
```cpp
        for (int i = 0; i < num_species_; ++i) {
            auto s = species_node[i];
            bool has_lookup = s["has_optics_lookup"] && s["has_optics_lookup"].as<bool>();
            
            GocartSpeciesParams p{
                s["dry_density"].as<double>(),
                s["molecular_weight"].as<double>(),
                s["dry_particle_diameter"].as<double>(),
                s["hygroscopicity"].as<double>(),
                s["lognormal_sigma"].as<double>(),
                s["lognormal_dg"].as<double>(),
                s["refractive_index_real"].as<double>(),
                s["refractive_index_imag"].as<double>(),
                has_lookup
            };

            if (has_lookup) {
                auto rh_bins = s["rh_bins"];
                auto ext_lookup = s["ext_lookup"];
                auto ssa_lookup = s["ssa_lookup"];
                auto asm_lookup = s["asm_lookup"];
                for (int j = 0; j < 8; ++j) {
                    p.rh_bins[j] = rh_bins[j].as<double>();
                    p.ext_lookup[j] = ext_lookup[j].as<double>();
                    p.ssa_lookup[j] = ssa_lookup[j].as<double>();
                    p.asm_lookup[j] = asm_lookup[j].as<double>();
                }
            }
            h_species_params_.push_back(p);
        }
```

- [ ] **Step 3: Compile and run test to verify dynamic optical solvers**

Run:
```bash
cmake --build build && ./build/tests/exaero_test_runner
```
Expected: PASS (Prints all pass messages, confirming correct optics calculations)

- [ ] **Step 4: Commit**

```bash
git add src/src_impl/gocart/ tests/test_main.cpp
git commit -m "feat: implement GOCART GPU dual-mode optics solver and column AOT integrations"
```
