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

    // Verify 3D diagnostics
    assert(diags_raw[exaero::diagnostic_indices::MASS_CONCENTRATION] == 1.0e-6);
    assert(diags_raw[exaero::diagnostic_indices::PM2_5_CONCENTRATION] == 1.0e-6); // 2.0e-6 dry size is <= 2.5e-6
    assert(diags_raw[exaero::diagnostic_indices::PM10_CONCENTRATION] == 1.0e-6);
    assert(diags_raw[exaero::diagnostic_indices::NUMBER_CONCENTRATION] > 0.0);
    assert(diags_raw[exaero::diagnostic_indices::SURFACE_AREA_DENSITY] > 0.0);

    // Verify 2D Column Mass and Surface Mass diagnostics
    double expected_col_mass = state_raw[0] * thick_raw[0]; // 1.0e-6 kg/m³ * 100 m = 1.0e-4 kg/m²
    assert(std::abs(diags_raw[exaero::diagnostic_indices::COLUMN_MASS] - expected_col_mass) < 1e-12);
    assert(diags_raw[exaero::diagnostic_indices::SURFACE_MASS] == 1.0e-6);
    assert(diags_raw[exaero::diagnostic_indices::SURFACE_PM2_5_MASS] == 1.0e-6);
    assert(std::abs(diags_raw[exaero::diagnostic_indices::COLUMN_PM2_5_MASS] - expected_col_mass) < 1e-12);

    std::cout << "Memory Mapping Zero-Copy Integration Test: PASS" << std::endl;
}

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

int main() {
    exaero::initialize_environment();
    
    test_gocart_yaml_parsing();
    test_zero_copy_mapping();
    test_gocart_optics();
    
    exaero::finalize_environment();
    return 0;
}
