#pragma once

namespace exaero {

    struct GocartSpeciesParams {
        double dry_density;             // [kg/m³]
        double molecular_weight;        // [g/mol]
        double dry_particle_diameter;   // [m]
        double hygroscopicity;          // kappa value
        double lognormal_sigma;         // GSD
        double lognormal_dg;            // GMD [m]
        double refractive_index_real;   // n
        double refractive_index_imag;   // k
    };

} // namespace exaero
