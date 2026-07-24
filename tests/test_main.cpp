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
    assert(diags_raw[exaero::diagnostic_indices::AEROSOL_LIQUID_WATER] > 0.0); // Absorbed ALW
    assert(diags_raw[exaero::diagnostic_indices::GRAVITATIONAL_SETTLING_VELOCITY] > 0.0); // vg fall velocity

    // Verify 2D Column Mass and Surface Mass diagnostics
    double expected_col_mass = state_raw[0] * thick_raw[0]; // 1.0e-6 kg/m³ * 100 m = 1.0e-4 kg/m²
    assert(std::abs(diags_raw[exaero::diagnostic_indices::COLUMN_MASS] - expected_col_mass) < 1e-12);
    assert(diags_raw[exaero::diagnostic_indices::SURFACE_MASS] == 1.0e-6);
    assert(diags_raw[exaero::diagnostic_indices::SURFACE_PM2_5_MASS] == 1.0e-6);
    assert(std::abs(diags_raw[exaero::diagnostic_indices::COLUMN_PM2_5_MASS] - expected_col_mass) < 1e-12);

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

    // Helper lambda to index 4D optics outputs: (cell, level, band, optics)
    auto get_optics_val = [&](int i_cell, int i_level, int i_band, int i_opt) {
        int idx = i_cell * (1 * 2 * exaero::optical_indices::NUM_OPTICS) +
                  i_level * (2 * exaero::optical_indices::NUM_OPTICS) +
                  i_band * (exaero::optical_indices::NUM_OPTICS) +
                  i_opt;
        return optics_raw[idx];
    };

    // ---- Mode B: RH Lookup Table Verification ----
    // rh_bins[2] = 0.7 (ext = 5.0), rh_bins[3] = 0.8 (ext = 6.5)
    // At RH = 75% (exactly midway), expected Sulfate MEE = 5.75 m²/g.
    // Conversion: 1.0e-6 kg/m³ * 5.75 m²/g * 1000 g/kg = 5.75e-3 m⁻¹ (Extinction Coefficient)
    double expected_sulfate_ext = 5.75e-3;

    // ---- Mode A: ADT Size Parameter Wavelength Scaling Verification ----
    // Longer wavelength (870 nm) has a smaller size parameter x = pi*D_wet/lambda
    // and therefore a smaller extinction coefficient than shorter wavelength (550 nm).
    double ext_550 = get_optics_val(0, 0, 0, exaero::optical_indices::EXTINCTION_COEFF);
    double ext_870 = get_optics_val(0, 0, 1, exaero::optical_indices::EXTINCTION_COEFF);

    assert(ext_550 > expected_sulfate_ext); // Total visible includes ADT Dust + Lookup Sulfate
    assert(ext_870 > 0.0);
    assert(ext_550 > ext_870); // Shorter wavelength has physically larger extinction!
    assert(get_optics_val(0, 0, 0, exaero::optical_indices::LIDAR_BACKSCATTER) > 0.0); // Lidar backscatter populated!

    // Verify 2D Column AOT calculations
    double aot_550 = get_optics_val(0, 0, 0, exaero::optical_indices::EXTINCTION_AOT);
    double aot_870 = get_optics_val(0, 0, 1, exaero::optical_indices::EXTINCTION_AOT);
    assert(aot_550 > 0.0);
    assert(aot_870 > 0.0);
    assert(aot_550 > aot_870); // Total visible AOT is physically larger than NIR AOT

    std::cout << "GOCART Optics Dual-Mode Calculations: PASS" << std::endl;
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

    double ccn_05ss = ccn_raw[0]; // Activated particles/m3 at 0.05% SS
    double ccn_10ss = ccn_raw[1]; // Activated particles/m3 at 0.1% SS
    double ccn_50ss = ccn_raw[2]; // Activated particles/m3 at 0.5% SS

    // Verify physical activation spectrum monotonicity:
    // Higher supersaturation must activate larger or equal number concentrations!
    assert(ccn_50ss >= ccn_10ss);
    assert(ccn_10ss >= ccn_05ss);

    std::cout << "GOCART Cloud CCN Activation Spectra: PASS" << std::endl;
}

int main() {
    exaero::initialize_environment();
    
    test_gocart_yaml_parsing();
    test_zero_copy_mapping();
    test_gocart_optics();
    test_gocart_ccn();
    
    exaero::finalize_environment();
    return 0;
}
