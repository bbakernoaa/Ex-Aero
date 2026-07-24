#include <gocart/GocartPackage.hpp>
#include <exaero/Environment.hpp>
#include <cassert>
#include <iostream>

void test_gocart_yaml_parsing() {
    exaero::initialize_environment();
    {
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
    exaero::finalize_environment();
}

int main() {
    test_gocart_yaml_parsing();
    return 0;
}
