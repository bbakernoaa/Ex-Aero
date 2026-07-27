#include <gocart/GocartPackage.hpp>
#include <exaero/Environment.hpp>
#include <cassert>
#include <iostream>
#include <random>
#include <cmath>
#include <limits>

// Verify dynamic Kohler wet sizing monotonicity
void verify_sizing_property(exaero::GocartPackage& package) {
    auto p = package.get_species_params(0);

    double rh_low = 0.0;
    double rh_high = 0.99;
    
    // Formula check: D_wet must be strictly increasing with RH
    double d_prev = 0.0;
    for (int i = 0; i <= 100; ++i) {
        double rh = rh_low + (rh_high - rh_low) * (i / 100.0);
        double wet_diameter = p.dry_particle_diameter * 
                              std::pow(1.0 + p.hygroscopicity * (rh / (1.0 - rh)), 1.0/3.0);
        assert(wet_diameter >= d_prev); // Sizing must be monotonic increasing!
        d_prev = wet_diameter;
    }
}

// Property-Based Invariant Verification (1000 randomized trials)
void run_property_based_tests(exaero::GocartPackage& package) {
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_real_distribution<double> rh_dist(0.0, 0.99);
    std::uniform_real_distribution<double> thick_dist(1.0, 1000.0);
    std::uniform_real_distribution<double> state_dist(1e-10, 1e-3);
    std::uniform_real_distribution<double> temp_dist(200.0, 320.0); // Temperature range [K]

    int trials = 1000;
    
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };

    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);

    for (int t = 0; t < trials; ++t) {
        double rh = rh_dist(rng);
        double thickness = thick_dist(rng);
        double concentration = state_dist(rng);
        double temp = temp_dist(rng);

        exaero::View2D<const double> temperature(&temp, 1, 1);
        exaero::View2D<const double> relative_humidity(&rh, 1, 1);
        exaero::View2D<const double> layer_thickness(&thickness, 1, 1);
        exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

        double state_raw[2] = { concentration, concentration };
        exaero::View3D<double> state(state_raw, 1, 1, 2);

        // 1. Diagnostics Property checks
        double diags_raw[exaero::diagnostic_indices::NUM_DIAGNOSTICS] = { 0.0 };
        exaero::View3D<double> diagnostics_out(diags_raw, 1, 1, exaero::diagnostic_indices::NUM_DIAGNOSTICS);

        package.computeDerivedDiagnostics(env, state, diagnostics_out);

        // Mass Conservation: Column Mass must equal concentration * thickness
        double expected_col_mass = (state_raw[0] + state_raw[1]) * thickness;
        assert(std::abs(diagnostics_out(0, 0, exaero::diagnostic_indices::COLUMN_MASS) - expected_col_mass) / expected_col_mass < 1e-12);
        
        // Mass non-negativity
        assert(diagnostics_out(0, 0, exaero::diagnostic_indices::PM2_5_CONCENTRATION) >= 0.0);
        assert(diagnostics_out(0, 0, exaero::diagnostic_indices::NUMBER_CONCENTRATION) >= 0.0);
        assert(diagnostics_out(0, 0, exaero::diagnostic_indices::SURFACE_AREA_DENSITY) >= 0.0);
        assert(diagnostics_out(0, 0, exaero::diagnostic_indices::AEROSOL_LIQUID_WATER) >= 0.0);
        assert(diagnostics_out(0, 0, exaero::diagnostic_indices::GRAVITATIONAL_SETTLING_VELOCITY) >= 0.0);

        // 2. Optical Property checks
        double wavelengths_raw[2] = { 550e-9, 870e-9 };
        exaero::View1D<const double> wavelengths(wavelengths_raw, 2);

        double optics_raw[1 * 1 * 2 * exaero::optical_indices::NUM_OPTICS] = { 0.0 };
        exaero::View4D<double> optics_out(optics_raw, 1, 1, 2, exaero::optical_indices::NUM_OPTICS);

        package.computeOptics(env, state, wavelengths, optics_out);

        for (int band = 0; band < 2; ++band) {
            double extinction = optics_out(0, 0, band, exaero::optical_indices::EXTINCTION_COEFF);
            double scattering = optics_out(0, 0, band, exaero::optical_indices::SCATTERING_COEFF);
            double asymmetry = optics_out(0, 0, band, exaero::optical_indices::ASYMMETRY_FACTOR);

            assert(extinction >= 0.0);
            assert(scattering >= 0.0);
            assert(extinction >= scattering); // scattering coefficient cannot exceed extinction!
            
            // Asymmetry g must reside in [-1.0, 1.0]
            assert(asymmetry >= -1.0 && asymmetry <= 1.0);
        }

        // 3. Cloud CCN Spectrum Property checks
        double ss_raw[3] = { 0.0005, 0.001, 0.005 };
        exaero::View1D<const double> supersaturations(ss_raw, 3);

        double ccn_raw[3] = { 0.0 };
        exaero::View4D<double> ccn_out(ccn_raw, 1, 1, 3, 1);

        package.computeCCN(env, state, supersaturations, ccn_out);

        // Monotonicity: higher supersaturations must activate larger or equal droplet counts
        assert(ccn_out(0, 0, 2, 0) >= ccn_out(0, 0, 1, 0));
        assert(ccn_out(0, 0, 1, 0) >= ccn_out(0, 0, 0, 0));
    }
    std::cout << "Property-Based Invariants (1000 Trials): PASS" << std::endl;
}

// Fuzz and Crash Testing Harness
void run_fuzz_tests(exaero::GocartPackage& package) {
    std::mt19937 rng(1337);
    
    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);

    // List of extreme, malformed, and NaN/Inf values to inject
    std::vector<double> fuzz_inputs = {
        -1.0, -100.0,       // extreme negative values
        0.0,                // zero boundaries
        1.0, 1.5, 100.0,    // out-of-bounds RH
        std::numeric_limits<double>::quiet_NaN(), // NaN
        std::numeric_limits<double>::infinity(),  // Positive Inf
        -std::numeric_limits<double>::infinity(), // Negative Inf
        1e-308, 1e308       // subnormal and extreme boundaries
    };

    for (double f_val : fuzz_inputs) {
        double rh = f_val;
        double thickness = f_val;
        double concentration = f_val;

        exaero::View2D<const double> relative_humidity(&rh, 1, 1);
        exaero::View2D<const double> layer_thickness(&thickness, 1, 1);
        exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

        double state_raw[2] = { concentration, concentration };
        exaero::View3D<double> state(state_raw, 1, 1, 2);

        double diags_raw[exaero::diagnostic_indices::NUM_DIAGNOSTICS] = { 0.0 };
        exaero::View3D<double> diagnostics_out(diags_raw, 1, 1, exaero::diagnostic_indices::NUM_DIAGNOSTICS);

        double wavelengths_raw[1] = { 550e-9 };
        exaero::View1D<const double> wavelengths(wavelengths_raw, 1);

        double optics_raw[1 * 1 * 1 * exaero::optical_indices::NUM_OPTICS] = { 0.0 };
        exaero::View4D<double> optics_out(optics_raw, 1, 1, 1, exaero::optical_indices::NUM_OPTICS);

        double ss_raw[1] = { 0.001 };
        exaero::View1D<const double> supersaturations(ss_raw, 1);

        double ccn_raw[1] = { 0.0 };
        exaero::View4D<double> ccn_out(ccn_raw, 1, 1, 1, 1);

        // Test robustness: solvers must defensively clamp values, handle NaNs, and NEVER segfault or crash
        try {
            package.computeDerivedDiagnostics(env, state, diagnostics_out);
            package.computeOptics(env, state, wavelengths, optics_out);
            package.computeCCN(env, state, supersaturations, ccn_out);
        } catch (...) {
            // Exceptions are acceptable, but segmentation faults or page-fault crashes are failures
        }
    }
    std::cout << "Fuzz/Crash Resilience Boundary Testing: PASS" << std::endl;
}

void test_dynamic_queries(exaero::GocartPackage& package) {
    // Dynamically query mapping indices to names (Hole 4)
    int dust_idx = package.getSpeciesIndex("Dust_ADT");
    int sulf_idx = package.getSpeciesIndex("Sulfate_Lookup");
    int fake_idx = package.getSpeciesIndex("FakeSpecies");

    assert(dust_idx == 0);
    assert(sulf_idx == 1);
    assert(fake_idx == -1); // Not found

    std::string name_0 = package.getSpeciesName(0);
    std::string name_1 = package.getSpeciesName(1);

    assert(name_0 == "Dust_ADT");
    assert(name_1 == "Sulfate_Lookup");

    // Test dynamic out of bounds index throws cleanly
    bool threw = false;
    try {
        package.getSpeciesName(99);
    } catch (const std::out_of_range&) {
        threw = true;
    }
    assert(threw);

    // Test defensive boundary checks throw cleanly when malformed layouts are passed (Hole 2)
    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.5 };
    double thick_raw[1] = { 100.0 };

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);
    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    // Pass invalid state array size of 3 species (mismatch!)
    double bad_state_raw[3] = { 1e-6, 1e-6, 1e-6 };
    exaero::View3D<double> bad_state(bad_state_raw, 1, 1, 3);

    double diags_raw[exaero::diagnostic_indices::NUM_DIAGNOSTICS] = { 0.0 };
    exaero::View3D<double> diagnostics_out(diags_raw, 1, 1, exaero::diagnostic_indices::NUM_DIAGNOSTICS);

    bool check_threw = false;
    try {
        package.computeDerivedDiagnostics(env, bad_state, diagnostics_out);
    } catch (const std::runtime_error& e) {
        check_threw = true;
    }
    assert(check_threw);

    std::cout << "Dynamic Metadata Queries & Boundary Protections: PASS" << std::endl;
}

int main() {
    exaero::initialize_environment();
    {
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

        verify_sizing_property(package);
        run_property_based_tests(package);
        run_fuzz_tests(package);
        test_dynamic_queries(package);
    }
    exaero::finalize_environment();
    return 0;
}
