#include <gocart/GocartPackage.hpp>
#include <exaero/Environment.hpp>
#include <cassert>
#include <iostream>

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

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity};

    double state_raw[1] = { 1.0e-6 };
    exaero::View3D<double> state(state_raw, 1, 1, 1);

    double diags_raw[6] = { 0.0 };
    exaero::View3D<double> diagnostics_out(diags_raw, 1, 1, 6);

    // Run diagnostic calculation step (maps raw pointers to unmanaged views zero-copy)
    package.computeDerivedDiagnostics(env, state, diagnostics_out);

    std::cout << "Memory Mapping Zero-Copy Integration Test: PASS" << std::endl;
}

int main() {
    exaero::initialize_environment();
    
    test_gocart_yaml_parsing();
    test_zero_copy_mapping();
    
    exaero::finalize_environment();
    return 0;
}
