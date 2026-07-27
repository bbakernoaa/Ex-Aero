#include <gocart/GocartPackage.hpp>
#include <exaero/Environment.hpp>
#include <cassert>
#include <iostream>
#include <cmath>

void test_gocart_yaml_parsing() {
    std::string yaml_string = R"(
    species:
      - name: "Dust"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    assert(package.get_num_species() == 1);
    auto p = package.get_species_params(0);
    assert(p.dry_density == 2600.0);
    assert(p.molecular_weight == 100.0);
    assert(p.dry_particle_diameter == 2.0e-6);
    assert(p.hygroscopicity == 0.1);
    assert(p.lognormal_sigma == 1.5);
    assert(p.lognormal_dg == 1.0e-6);
    assert(p.refractive_index_real == 1.55);
    assert(p.refractive_index_imag == 0.002);

    std::cout << "YAML Parsing Unit Test: PASS" << std::endl;
}

void test_zero_copy_mapping() {
    std::string yaml_string = R"(
    species:
      - name: "Dust"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.5 };
    double thick_raw[1] = { 100.0 }; // 100 meters vertical layer thickness

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    double state_raw[1] = { 1.0e-6 };
    exaero::View3D<double> state(state_raw, 1, 1, 1);

    double diags_raw[exaero::diagnostic_indices::NUM_DIAGNOSTICS] = { 0.0 };
    exaero::View3D<double> diagnostics_out(diags_raw, 1, 1, exaero::diagnostic_indices::NUM_DIAGNOSTICS);

    // Run diagnostic calculation step (maps raw pointers to unmanaged views zero-copy)
    package.computeDerivedDiagnostics(env, state, diagnostics_out);

    // Verify 3D diagnostics using dynamic mdspan bracket accessors
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::MASS_CONCENTRATION) == 1.0e-6);
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::PM2_5_CONCENTRATION) == 1.0e-6); // 2.0e-6 dry size is <= 2.5e-6
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::PM10_CONCENTRATION) == 1.0e-6);
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::NUMBER_CONCENTRATION) > 0.0);
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::SURFACE_AREA_DENSITY) > 0.0);
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::AEROSOL_LIQUID_WATER) > 0.0); // Absorbed ALW
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::GRAVITATIONAL_SETTLING_VELOCITY) > 0.0); // vg fall velocity

    // Verify 2D Column Mass and Surface Mass diagnostics
    double expected_col_mass = state_raw[0] * thick_raw[0]; // 1.0e-6 kg/m³ * 100 m = 1.0e-4 kg/m²
    assert(std::abs(diagnostics_out(0, 0, exaero::diagnostic_indices::COLUMN_MASS) - expected_col_mass) < 1e-12);
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::SURFACE_MASS) == 1.0e-6);
    assert(diagnostics_out(0, 0, exaero::diagnostic_indices::SURFACE_PM2_5_MASS) == 1.0e-6);
    assert(std::abs(diagnostics_out(0, 0, exaero::diagnostic_indices::COLUMN_PM2_5_MASS) - expected_col_mass) < 1e-12);

    std::cout << "Memory Mapping Zero-Copy Integration Test: PASS" << std::endl;
}

void test_gocart_optics() {
    // YAML containing two species: ADT Dust (scaled to small size) and Lookup-Table Sulfate
    std::string yaml_string = R"(
    species:
      - name: "Dust_ADT"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 0.15e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 0.1e-6
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
    double rh_raw[1] = { 0.75 }; // 75% relative humidity (exactly midway between bin 2: 0.7 and bin 3: 0.8!)
    double thick_raw[1] = { 100.0 }; // 100 meters grid height

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    // 1.0 ug/m³ of ADT Dust and 1.0 ug/m³ of Lookup Sulfate
    double state_raw[2] = { 1.0e-6, 1.0e-6 };
    exaero::View3D<double> state(state_raw, 1, 1, 2);

    // Query two wavelength bands simultaneously: 550nm and 870nm
    double wavelengths_raw[2] = { 550.0e-9, 870.0e-9 };
    exaero::View1D<const double> wavelengths(wavelengths_raw, 2);

    // 4D output array of size (cells, levels, bands, optical_indices::NUM_OPTICS)
    double optics_raw[1 * 1 * 2 * exaero::optical_indices::NUM_OPTICS] = { 0.0 };
    exaero::View4D<double> optics_out(optics_raw, 1, 1, 2, exaero::optical_indices::NUM_OPTICS);

    // Run multi-band optics calculations
    package.computeOptics(env, state, wavelengths, optics_out);

    // ---- Mode B: RH Lookup Table Verification ----
    // rh_bins[2] = 0.7 (ext = 5.0), rh_bins[3] = 0.8 (ext = 6.5)
    // At RH = 75% (exactly midway), expected Sulfate MEE = 5.75 m²/g.
    // Conversion: 1.0e-6 kg/m³ * 5.75 m²/g * 1000 g/kg = 5.75e-3 m⁻¹ (Extinction Coefficient)
    double expected_sulfate_ext = 5.75e-3;

    // ---- Mode A: ADT Size Parameter Wavelength Scaling Verification ----
    // Longer wavelength (870 nm) has a smaller size parameter x = pi*D_wet/lambda
    // and therefore a smaller extinction coefficient than shorter wavelength (550 nm).
    double ext_550 = optics_out(0, 0, 0, exaero::optical_indices::EXTINCTION_COEFF);
    double ext_870 = optics_out(0, 0, 1, exaero::optical_indices::EXTINCTION_COEFF);

    assert(ext_550 > expected_sulfate_ext); // Total visible includes ADT Dust + Lookup Sulfate
    assert(ext_870 > 0.0);
    assert(ext_550 > ext_870); // Shorter wavelength has physically larger extinction!
    assert(optics_out(0, 0, 0, exaero::optical_indices::LIDAR_BACKSCATTER) > 0.0); // Lidar backscatter populated!

    // Verify 2D Column AOT calculations
    double aot_550 = optics_out(0, 0, 0, exaero::optical_indices::EXTINCTION_AOT);
    double aot_870 = optics_out(0, 0, 1, exaero::optical_indices::EXTINCTION_AOT);
    assert(aot_550 > 0.0);
    assert(aot_870 > 0.0);
    assert(aot_550 > aot_870); // Total visible AOT is physically larger than NIR AOT

    std::cout << "GOCART Optics Dual-Mode Calculations: PASS" << std::endl;
}

void test_optical_precision() {
    // High-precision physical validation of our GPU ADT Mie solver against standard analytical results.
    // For n = 1.5, x = 10.0, we have:
    // phase shift rho = 2 * x * (n - 1) = 10.0
    // sin(10.0) ≈ -0.54402111, cos(10.0) ≈ -0.83907153
    // Analytical ADT Q_ext = 2 - (4/10)*sin(10) + (4/100)*(1 - cos(10)) = 2.29117135
    std::string yaml_string = R"(
    species:
      - name: "Dust_ADT_Precision"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.0            # No wet size growth (D_wet = D_dry = 2.0e-6)
        lognormal_sigma: 1.0001        # Monodisperse limit (ln(sigma) ≈ 0)
        lognormal_dg: 2.0e-6           # dg = dry_particle_diameter = 2.0e-6
        refractive_index_real: 1.50    # n = 1.5
        refractive_index_imag: 0.0     # Non-absorbing (k = 0)
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.0 };       // Dry (0% RH)
    double thick_raw[1] = { 1.0 };

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    // Calculate state mass concentration to yield exactly 1.0 particle/m³
    // N = C / ( (pi/6)*rho*D^3 * exp(4.5*ln(sigma)^2) )
    // For N = 1.0, C = (pi/6)*rho*D^3 ≈ (3.14159265/6)*2600*(2.0e-6)^3 ≈ 1.08908549e-14 [kg/m³]
    double pi_val = 3.141592653589793;
    double expected_vol = (pi_val / 6.0) * 2600.0 * std::pow(2.0e-6, 3);
    double state_raw[1] = { expected_vol };
    exaero::View3D<double> state(state_raw, 1, 1, 1);

    // Query wavelength that yields exactly size parameter x = 10.0
    // x = pi * D_wet / lambda -> lambda = pi * D_wet / 10.0 = pi * 2.0e-6 / 10.0 ≈ 6.2831853e-7 [m]
    double lambda_query = pi_val * 2.0e-6 / 10.0;
    double wavelengths_raw[1] = { lambda_query };
    exaero::View1D<const double> wavelengths(wavelengths_raw, 1);

    double optics_raw[exaero::optical_indices::NUM_OPTICS] = { 0.0 };
    exaero::View4D<double> optics_out(optics_raw, 1, 1, 1, exaero::optical_indices::NUM_OPTICS);

    // Run multi-band optics calculations
    package.computeOptics(env, state, wavelengths, optics_out);

    double total_ext_coeff = optics_out(0, 0, 0, exaero::optical_indices::EXTINCTION_COEFF);

    // Under N = 1.0 particle/m³, extinction coefficient is:
    // b_ext = N * cross_section * Q_ext
    // where cross_section = (pi/4)*D^2 = (3.14159265/4)*(2.0e-6)^2 ≈ 3.14159265e-12 [m²]
    // So Q_ext = b_ext / (N * cross_section)
    double cross_section = (pi_val / 4.0) * std::pow(2.0e-6, 2);
    double calculated_q_ext = total_ext_coeff / (1.0 * cross_section);

    double expected_q_ext = 2.29117135; // Van de Hulst Analytical Limit for rho = 10.0

    // Assert that our C++/Kokkos GPU ADT kernel calculates and returns exactly the expected Mie efficiency
    // with high floating-point numerical precision (< 1e-6 tolerance!)
    double numerical_tolerance = 1.0e-6;
    double numerical_diff = std::abs(calculated_q_ext - expected_q_ext);

    assert(numerical_diff < numerical_tolerance);

    std::cout << "GOCART Optics High-Precision Physical Validation: PASS" << std::endl;
    std::cout << "  - Expected Q_ext: " << expected_q_ext << std::endl;
    std::cout << "  - Calculated Q_ext: " << calculated_q_ext << " (Diff: " << numerical_diff << ")" << std::endl;
}

void test_gocart_ccn() {
    // YAML containing two species: ADT Dust (scaled to small size) and Lookup-Table Sulfate
    std::string yaml_string = R"(
    species:
      - name: "Dust_ADT"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 0.15e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 0.1e-6
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
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    double temp_raw[1] = { 298.0 }; // T = 298 K
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.50 };
    double thick_raw[1] = { 100.0 };

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    // 1.0 ug/m³ of ADT Dust and 1.0 ug/m³ of Sulfate
    double state_raw[2] = { 1.0e-6, 1.0e-6 };
    exaero::View3D<double> state(state_raw, 1, 1, 2);

    // Query three supersaturation levels simultaneously: 0.05% (0.0005), 0.1% (0.001), 0.5% (0.005)
    double ss_raw[3] = { 0.0005, 0.001, 0.005 };
    exaero::View1D<const double> supersaturations(ss_raw, 3);

    // 4D output view of size (cells, levels, ss, 1) -> we can use ccn_out view directly
    double ccn_raw[1 * 1 * 3] = { 0.0 };
    exaero::View4D<double> ccn_out(ccn_raw, 1, 1, 3, 1);

    // Run multi-supersaturation CCN activation solver
    package.computeCCN(env, state, supersaturations, ccn_out);

    double ccn_05ss = ccn_out(0, 0, 0, 0); // Activated particles/m3 at 0.05% SS
    double ccn_10ss = ccn_out(0, 0, 1, 0); // Activated particles/m3 at 0.1% SS
    double ccn_50ss = ccn_out(0, 0, 2, 0); // Activated particles/m3 at 0.5% SS

    // Verify physical activation spectrum monotonicity:
    // Higher supersaturation must activate larger or equal number concentrations!
    assert(ccn_50ss >= ccn_10ss);
    assert(ccn_10ss >= ccn_05ss);

    std::cout << "GOCART Cloud CCN Activation Spectra: PASS" << std::endl;
}

void test_gocart_yaml_emissions_parsing() {
    std::string yaml_string = R"(
    species:
      - name: "Dust_1"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 0.15e-6
        hygroscopicity: 0.14
        lognormal_sigma: 1.5
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.53
        refractive_index_imag: 0.003
    emissions_mapping:
      - raw_name: "CECE_Dust"
        mappings:
          - target_species: "Dust_1"
            mass_split_fraction: 0.35
            is_modal_mode: true
            emitted_particle_diameter: 0.25e-6
            lognormal_sigma: 1.8
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    assert(package.get_num_species() == 1);
    auto p = package.get_species_params(0);
    assert(p.emissions_mapping.is_active == true);
    assert(p.emissions_mapping.raw_cece_index == 0);
    assert(p.emissions_mapping.mass_split_fraction == 0.35);
    assert(p.emissions_mapping.is_modal_mode == true);
    assert(p.emissions_mapping.emitted_particle_diameter == 0.25e-6);
    assert(p.emissions_mapping.lognormal_sigma == 1.8);

    std::cout << "YAML Emissions Parsing Unit Test: PASS" << std::endl;
}

int main() {
    exaero::initialize_environment();
    
    test_gocart_yaml_parsing();
    test_gocart_yaml_emissions_parsing();
    test_zero_copy_mapping();
    test_gocart_optics();
    test_optical_precision();
    test_gocart_ccn();
    
    exaero::finalize_environment();
    return 0;
}
