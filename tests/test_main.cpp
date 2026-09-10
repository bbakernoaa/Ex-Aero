#include <gocart/GocartPackage.hpp>
#include <loader/MieTableStore.hpp>
#include <utils/LutGenerator.hpp>
#include <utils/DubovikSpheroidDatabase.hpp>
#include <exaero/Environment.hpp>
#include <exaero/GocartConfig.hpp>
#include <cassert>
#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <regex>
#include <algorithm>
#include <chrono>
#include <vector>

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

void test_gocart_arbitrary_lookup_length() {
    // ADR-003 R10 / invariant C12: the RH lookup path must honor ANY declared length.
    // The retired fixed-capacity arrays capped the bracket scan at 6 points; a prime
    // (7-point) and a short (3-point) table must now both interpolate correctly.
    std::string yaml_string = R"(
    species:
      - name: "SevenPoint"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.2e-6
        hygroscopicity: 0.50
        lognormal_sigma: 2.0
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-8
        has_optics_lookup: true
        rh_bins: [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0]
        ext_lookup: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        ssa_lookup: [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        asm_lookup: [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
      - name: "ThreePoint"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.2e-6
        hygroscopicity: 0.50
        lognormal_sigma: 2.0
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-8
        has_optics_lookup: true
        rh_bins: [0.0, 0.50, 0.99]
        ext_lookup: [2.0, 4.0, 8.0]
        ssa_lookup: [0.9, 0.9, 0.9]
        asm_lookup: [0.7, 0.7, 0.7]
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    // Extents are data, not capacity: 7 and 3 are both representable.
    assert(package.get_species_params(0).n_rh == 7);
    assert(package.get_species_params(1).n_rh == 3);

    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.5 };   // inside the 0.4-0.6 bracket of the 7-point table
    double thick_raw[1] = { 1.0 };
    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);
    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    double state_raw[2] = { 1.0e-6, 1.0e-6 };
    exaero::View3D<double> state(state_raw, 1, 1, 2);
    double wavelengths_raw[1] = { 550.0e-9 };
    exaero::View1D<const double> wavelengths(wavelengths_raw, 1);
    double optics_raw[1 * 1 * 1 * exaero::optical_indices::NUM_OPTICS] = { 0.0 };
    exaero::View4D<double> optics_out(optics_raw, 1, 1, 1, exaero::optical_indices::NUM_OPTICS);

    package.computeOptics(env, state, wavelengths, optics_out);
    double total_ext = optics_out(0, 0, 0, exaero::optical_indices::EXTINCTION_COEFF);
    // SevenPoint @ rh=0.5: MEE = 3.5 -> 3.5e-3 m^-1. ThreePoint @ rh=0.5: MEE = 4.0 -> 4.0e-3.
    double expected = 3.5e-3 + 4.0e-3;
    assert(std::abs(total_ext - expected) / expected < 1e-9);

    // Above the last declared edge the curve clamps to the final bracket (no extrapolation).
    rh_raw[0] = 0.99;
    package.computeOptics(env, state, wavelengths, optics_out);
    total_ext = optics_out(0, 0, 0, exaero::optical_indices::EXTINCTION_COEFF);
    // SevenPoint @ 0.99: bracket [0.9,1.0] -> 6.9e-3. ThreePoint @ 0.99: last point 8.0 -> 8.0e-3.
    expected = 6.9e-3 + 8.0e-3;
    assert(std::abs(total_ext - expected) / expected < 1e-9);

    std::cout << "GOCART Arbitrary Lookup Length (R10/C12): PASS" << std::endl;
}

void test_gocart_lookup_length_mismatch_fails() {
    // Fail-fast, no silent fallback (FR-009): mismatched lookup list lengths must abort.
    std::string yaml_string = R"(
    species:
      - name: "Broken"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.2e-6
        hygroscopicity: 0.5
        lognormal_sigma: 2.0
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-8
        has_optics_lookup: true
        rh_bins: [0.0, 0.5, 0.99]
        ext_lookup: [2.0, 4.0]
        ssa_lookup: [0.9, 0.9, 0.9]
        asm_lookup: [0.7, 0.7, 0.7]
    )";
    exaero::GocartPackage package;
    bool threw = false;
    try {
        package.initialize(yaml_string);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw); // must abort, never silently pad or truncate
    std::cout << "GOCART Lookup Length Mismatch Fail-Fast: PASS" << std::endl;
}

// --- GEOSmie MIE attribute surface: US1 microphysical golden tests (T015, T016) ---
// Golden values are the pinned snapshot's float64-widened grid values (FR-005, C1);
// provenance must carry unit + version + citation + delivery source on every value (C5).

// Relative tolerance for grid-point reproduction (FR-005).
static const double kGridTol = 1e-7;

// Six baked species configured by their source labels so the package resolves each
// name straight to its baked curve (no mie_table binding needed for the golden path).
static const char* const kMieSixYaml = R"YAML(
    species:
      - name: "DU"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
      - name: "SS"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.2e-6
        hygroscopicity: 0.5
        lognormal_sigma: 2.0
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-8
      - name: "SU"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.2e-6
        hygroscopicity: 0.5
        lognormal_sigma: 2.0
        lognormal_dg: 0.15e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-8
      - name: "BC"
        dry_density: 1800.0
        molecular_weight: 12.0
        dry_particle_diameter: 0.1e-6
        hygroscopicity: 0.0
        lognormal_sigma: 1.8
        lognormal_dg: 0.05e-6
        refractive_index_real: 1.85
        refractive_index_imag: 0.75
      - name: "OC"
        dry_density: 1300.0
        molecular_weight: 150.0
        dry_particle_diameter: 0.1e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.8
        lognormal_dg: 0.05e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.0
      - name: "NI"
        dry_density: 2150.0
        molecular_weight: 85.0
        dry_particle_diameter: 0.3e-6
        hygroscopicity: 0.6
        lognormal_sigma: 1.6
        lognormal_dg: 0.1e-6
        refractive_index_real: 1.52
        refractive_index_imag: 0.01
)YAML";

namespace mie_micro {
struct Golden {
    const char* species;
    int index;                 // package species index
    double reff00;             // rEff(bin0, rh=0)      [m]
    double mass00;             // rMass(bin0, rh=0)     [kg]
    double rlow0, rupp0;       // bin boundaries        [m]
    double growth_rh50;        // rEff(bin0, rh=0.5) / rEff(bin0, rh=0)
    double wetdens00;          // rMass/(4/3 pi rEff^3) at (bin0, rh=0) [kg m^-3]
    double volmass00, areamass00;
};
// Golden values read from the pinned snapshot (tools/geosmie_snapshot) at full float64.
static const Golden GOLDEN[6] = {
    {"DU", 0, 6.358845325848961e-07, 5.644337028201749e-15, 1.0000000116860974e-07,
     9.999999974752427e-07, 1.0, 5240.70296164897, 0.00019081409637560404, 225.0574828419758},
    {"SS", 1, 7.841114069151445e-08, 4.9905266329105316e-18, 2.999999892949745e-08,
     1.0000000116860974e-07, 1.2618481798599361, 2471.2942166388125, 0.0004046462753269791,
     3870.4284087538726},
    {"SU", 2, 1.5664450359054172e-07, 4.6720689626767e-17, 4.999999969612645e-09,
     3.000000106112566e-07, 1.3904448902042543, 2901.8512928604996, 0.00034460759669536687,
     1649.9506308699547},
    {"BC", 3, 3.921032032394578e-08, 9.94582633186942e-19, 4.999999969612645e-09,
     3.000000106112566e-07, 1.0, 3938.683571448429, 0.00025389193669910767, 4856.347791885845},
    {"OC", 4, 8.765797332443981e-08, 1.6408043966108013e-17, 4.999999969612645e-09,
     3.000000106112566e-07, 1.0, 5815.591933351836, 0.0001719515418998194, 1471.2142151352748},
    {"NI", 5, 1.560005102874129e-07, 1.234258558954136e-16, 4.999999969612645e-09,
     0.0, 1.4733953629039604, 7761.389545995164, 0.0001288429080996192, 619.4350319539391},
};
} // namespace mie_micro

static void check_close(double got, double want, const char* what) {
    const double tol = kGridTol * (std::abs(want) > 0.0 ? std::abs(want) : 1.0);
    assert(std::abs(got - want) <= tol);
    (void)what;
}

void test_mie_baked_in_microphysical() {
    using namespace exaero::microphysical_indices;
    exaero::GocartPackage package;
    package.initialize(kMieSixYaml);

    for (const auto& g : mie_micro::GOLDEN) {
        double v = 0.0;
        exaero::ProvenanceInfo p{};
        auto st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                         EFFECTIVE_RADIUS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.reff00, g.species);

        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    PARTICLE_MASS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.mass00, g.species);

        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    BIN_LOWER_RADIUS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.rlow0, g.species);

        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    BIN_UPPER_RADIUS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.rupp0, g.species);

        // rh = 0.5 is an exact source grid point (rh[10] == 0.5): no interpolation.
        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    GROWTH_FACTOR, 0.5, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        assert(p.interpolated == 0);
        check_close(v, g.growth_rh50, g.species);

        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    WET_PARTICLE_DENSITY, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.wetdens00, g.species);

        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    VOLUME_PER_MASS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.volmass00, g.species);

        st = package.queryAttribute(g.index, exaero::AttributeCategory::Microphysical,
                                    AREA_PER_MASS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, g.areamass00, g.species);
    }
    std::cout << "MIE Baked-in Microphysical Grid Values (FR-005/C1): PASS" << std::endl;
}

void test_mie_provenance() {
    using namespace exaero::microphysical_indices;
    exaero::GocartPackage package;
    package.initialize(kMieSixYaml);

    double v = 0.0;
    exaero::ProvenanceInfo p{};
    package.queryAttribute(0, exaero::AttributeCategory::Microphysical, EFFECTIVE_RADIUS,
                           0.0, 0.0, &v, &p);
    assert(std::string(p.unit) == "m");
    assert(std::string(p.species) == "DU");
    assert(std::string(p.source_version) == "ufs-regtests-input-data-20260617/GOCART/p8c_5d");
    assert(std::strlen(p.citation) > 0);
    assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::BakedIn));
    assert(p.interpolated == 0);
    assert(p.num_radius == 5 && p.num_rh == 36 && p.num_lambda == 30); // DU extents are data

    package.queryAttribute(0, exaero::AttributeCategory::Microphysical, WET_PARTICLE_DENSITY,
                           0.0, 0.0, &v, &p);
    assert(std::string(p.unit) == "kg m^-3");
    package.queryAttribute(0, exaero::AttributeCategory::Microphysical, GROWTH_FACTOR,
                           0.0, 0.0, &v, &p);
    assert(std::string(p.unit) == "1");
    package.queryAttribute(0, exaero::AttributeCategory::Microphysical, PARTICLE_MASS,
                           0.0, 0.0, &v, &p);
    assert(std::string(p.unit) == "kg");

    // mass_mean_radius has no source variable in the pinned band tables (data-model
    // implement-time reconciliation): explicit NotInSource, never a silent 0 (FR-008).
    auto st = package.queryAttribute(0, exaero::AttributeCategory::Microphysical,
                                     MASS_MEAN_RADIUS, 0.0, 0.0, &v, &p);
    assert(st == exaero::AttributeStatus::NotInSource);

    std::cout << "MIE Provenance (unit+version+citation+delivery, C5): PASS" << std::endl;
}

void test_mie_not_available() {
    using namespace exaero::microphysical_indices;
    // A species with no bound curve: explicit NotInSource, value untouched (FR-008/C4).
    std::string yaml = std::string(kMieSixYaml) +
        R"YAML(      - name: "MYSTERY"
        dry_density: 1000.0
        molecular_weight: 50.0
        dry_particle_diameter: 0.1e-6
        hygroscopicity: 0.0
        lognormal_sigma: 1.5
        lognormal_dg: 0.1e-6
        refractive_index_real: 1.4
        refractive_index_imag: 0.0
)YAML";
    exaero::GocartPackage package;
    package.initialize(yaml);

    double v = -777.0; // sentinel: must NOT be overwritten with a silent 0
    exaero::ProvenanceInfo p{};
    auto st = package.queryAttribute(6, exaero::AttributeCategory::Microphysical,
                                     EFFECTIVE_RADIUS, 0.0, 0.0, &v, &p);
    assert(st == exaero::AttributeStatus::NotInSource);
    assert(v == -777.0);

    // Deactivate microphysical for DU only: explicit NotActivated (FR-010/C4).
    package.setAttributeActivation({"DU"},
        exaero::attribute_category_bit(exaero::AttributeCategory::SpectralOptical), "");
    st = package.queryAttribute(0, exaero::AttributeCategory::Microphysical,
                                EFFECTIVE_RADIUS, 0.0, 0.0, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotActivated);
    // Spectral stays active for DU.
    st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                exaero::spectral_optical_indices::EXTINCTION_EFFICIENCY,
                                0.0, 1.0, &v, nullptr);
    assert(st == exaero::AttributeStatus::Available);

    // Out-of-range species index: explicit status, no crash, no silent 0.
    st = package.queryAttribute(999, exaero::AttributeCategory::Microphysical,
                                EFFECTIVE_RADIUS, 0.0, 0.0, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);

    // Restore full activation so later tests see the default surface.
    package.setAttributeActivation({},
        exaero::attribute_category_bit(exaero::AttributeCategory::Microphysical) |
        exaero::attribute_category_bit(exaero::AttributeCategory::SpectralOptical), "");

    std::cout << "MIE Explicit Not-Available Statuses (FR-008/C4): PASS" << std::endl;
}

// Spectral RRTMG-band goldens at (radius bin 0, rh=0.5 grid point, band index 3.0),
// read from the pinned snapshot at full float64 (ravel order b*nH*nL + h*nL + l).
namespace mie_spectral {
struct SpectralGolden {
    const char* species;
    int index;
    double qext, qsca, qabs, bext, bsca, bbck, lidar, g, ssa, n, k;
};
static const SpectralGolden GOLDEN[6] = {
    {"DU", 0, 1.9339340925216675, 1.821953296661377, 0.11198079586029053,
     829.6609497070312, 779.1709594726562, 10.848725318908691, 76.47543147405354,
     0.5862439870834351, 0.9420968913608229, 1.5232499837875366, -0.011869008652865887},
    {"SS", 1, 0.00227008992806077, 0.0009156085434369743, 0.0013544813846237957,
     15.748342514038086, 6.351870536804199, 0.723405122756958, 21.769741488724634,
     0.017153162509202957, 0.4033358027446671, 1.3590805530548096, -0.002306509530171752},
    {"SU", 2, 0.029475990682840347, 0.026380717754364014, 0.0030952729284763336,
     160.50079345703125, 143.6468963623047, 10.787318229675293, 14.878655662118392,
     0.15355892479419708, 0.8949900289431741, 1.3059253692626953, -0.002059066668152809},
    {"BC", 3, 0.08555283397436142, 0.001218374352902174, 0.08433445962145925,
     1636.421630859375, 23.304594039916992, 2.423806667327881, 675.1452799094094,
     0.05321965739130974, 0.01424119221190613, 1.813500165939331, -0.5034999847412109},
    {"OC", 4, 0.00932903029024601, 0.0036768026184290648, 0.005652227671816945,
     44.34389877319336, 17.477035522460938, 1.7159570455551147, 25.84207972341641,
     0.07141537964344025, 0.39412484513780116, 1.419995665550232, -0.010700000450015068},
    {"NI", 5, 0.08747505396604538, 0.08508098870515823, 0.002394065260887146,
     529.2645874023438, 514.7792358398438, 14.933114051818848, 35.44234548572804,
     0.5163609981536865, 0.9726314514556752, 1.3417788743972778, -0.001385296112857759},
};
} // namespace mie_spectral

void test_mie_rrtmg_spectral() {
    using namespace exaero::spectral_optical_indices;
    exaero::GocartPackage package;
    package.initialize(kMieSixYaml);

    // rh=0.5 and band index 3.0 are both exact source grid points -> Available, no interp.
    const double rh = 0.5, band = 3.0;
    for (const auto& g : mie_spectral::GOLDEN) {
        double v = 0.0;
        exaero::ProvenanceInfo p{};
        auto q = [&](int idx, double want) {
            auto st = package.queryAttribute(g.index, exaero::AttributeCategory::SpectralOptical,
                                             idx, rh, band, &v, &p);
            assert(st == exaero::AttributeStatus::Available);
            assert(p.interpolated == 0);
            check_close(v, want, g.species);
        };
        q(EXTINCTION_EFFICIENCY, g.qext);
        q(SCATTERING_EFFICIENCY, g.qsca);
        q(ABSORPTION_EFFICIENCY, g.qabs);
        q(MASS_EXTINCTION, g.bext);
        q(MASS_SCATTERING, g.bsca);
        q(MASS_BACKSCATTER, g.bbck);
        q(LIDAR_RATIO, g.lidar);
        q(ASYMMETRY_FACTOR, g.g);
        q(SINGLE_SCATTERING_ALBEDO, g.ssa);
        q(REFRACTIVE_INDEX_REAL, g.n);
        q(REFRACTIVE_INDEX_IMAG, g.k);

        // Physical bounds (SC-005/C3): ssa in [0,1], g in [-1,1], qext/qsca/bext >= 0.
        double ssa_v = 0.0, g_v = 0.0, qe = 0.0;
        package.queryAttribute(g.index, exaero::AttributeCategory::SpectralOptical,
                               SINGLE_SCATTERING_ALBEDO, rh, band, &ssa_v, nullptr);
        package.queryAttribute(g.index, exaero::AttributeCategory::SpectralOptical,
                               ASYMMETRY_FACTOR, rh, band, &g_v, nullptr);
        package.queryAttribute(g.index, exaero::AttributeCategory::SpectralOptical,
                               EXTINCTION_EFFICIENCY, rh, band, &qe, nullptr);
        assert(ssa_v >= -1e-6 && ssa_v <= 1.0 + 1e-6);
        assert(g_v >= -1.0 - 1e-6 && g_v <= 1.0 + 1e-6);
        assert(qe >= 0.0);
    }
    std::cout << "MIE RRTMG Spectral Grid Values (FR-002/FR-005/C1): PASS" << std::endl;
}

void test_mie_interpolation() {
    using namespace exaero::spectral_optical_indices;
    exaero::GocartPackage package;
    package.initialize(kMieSixYaml);

    // SS radius bin 1 (solver default is 0, so address via the store directly through a
    // config curve bound to bin 1). SS is hygroscopic: qext varies in RH and band index.
    // Corners at (b=1, rh in {0.5,0.550000011920929}, band in {3.0,4.0}):
    //   v00=0.1515267938375473 v10=0.15872181951999664
    //   v01=0.22853879630565643 v11=0.24023324251174927
    // Reference (research R3): linear in RH, LINEAR-IN-LOG band index.
    auto& store = exaero::MieTableStore::instance();
    double v = 0.0;
    exaero::Provenance prov{};

    struct Case { const char* what; double rh, band, want; };
    const Case cases[] = {
        {"rh-only (0.52, band 3.0)", 0.52, 3.0, 0.15440480342435609},
        {"band-only log (0.5, band 3.5)", 0.5, 3.5, 0.192792669163537},
        {"both (0.52, band 3.5)", 0.52, 3.5, 0.1966350608006735},
    };
    for (const auto& cse : cases) {
        auto st = store.query("SS", exaero::AttributeCategory::SpectralOptical,
                              EXTINCTION_EFFICIENCY, cse.rh, cse.band, /*radius_index=*/1,
                              &v, &prov);
        assert(st == exaero::AttributeStatus::Interpolated);
        assert(prov.interpolated == 1);
        // Tolerance 1e-5 RELATIVE (FR-006/C2).
        const double tol = 1e-5 * std::abs(cse.want);
        assert(std::abs(v - cse.want) <= tol);
    }
    std::cout << "MIE Off-Grid Interpolation (FR-006/C2, log-wavelength): PASS" << std::endl;
}

void test_mie_compute_attributes_multiband() {
    using namespace exaero::spectral_optical_indices;
    exaero::GocartPackage package;
    package.initialize(kMieSixYaml);

    // Three bands: two exact grid points (band 3, band 4) and one off-grid (3.5).
    // Output is band-major: slot = band * NUM_ATTRIBUTES + attribute (T027).
    const int num_cells = 1, num_levels = 1;
    double temp[1] = {298.0}, pres[1] = {101325.0}, dens[1] = {1.2};
    double rh_raw[1] = {0.5}, thick[1] = {100.0}, state_raw[1] = {1.0e-6};
    exaero::View2D<const double> temperature(temp, num_cells, num_levels);
    exaero::View2D<const double> pressure(pres, num_cells, num_levels);
    exaero::View2D<const double> air_density(dens, num_cells, num_levels);
    exaero::View2D<const double> relative_humidity(rh_raw, num_cells, num_levels);
    exaero::View2D<const double> layer_thickness(thick, num_cells, num_levels);
    exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                       relative_humidity, layer_thickness};
    exaero::View3D<const double> state(state_raw, num_cells, num_levels, 1);

    const int num_bands = 3;
    double wavelengths_raw[num_bands] = {3.0, 4.0, 3.5};
    exaero::View1D<const double> wavelengths(wavelengths_raw, num_bands);

    const int slots = num_bands * NUM_ATTRIBUTES;
    double attrs[3 * 11];
    int statuses[3 * 11];
    for (double& x : attrs) x = -777.0;
    for (int& s : statuses) s = -1;
    exaero::View3D<double> attributes_out(attrs, num_cells, num_levels, slots);
    exaero::View3D<int> status_out(statuses, num_cells, num_levels, slots);

    package.computeAttributes(env, state, /*species_index=*/1 /*SS*/,
                              exaero::AttributeCategory::SpectralOptical,
                              wavelengths, attributes_out, &status_out);

    // Band 3 (grid): SS qext at (bin0, rh=0.5, band 3) golden.
    check_close(attrs[0 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY], 0.00227008992806077, "b3 qext");
    assert(statuses[0 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY] ==
           static_cast<int>(exaero::AttributeStatus::Available));
    // Band 4 (grid): must differ from band 3 -> proves each band is resolved, not just band 0.
    const double qe_b4 = attrs[1 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY];
    assert(qe_b4 != attrs[0 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY]);
    assert(std::isfinite(qe_b4));
    assert(statuses[1 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY] ==
           static_cast<int>(exaero::AttributeStatus::Available));
    // Band 3.5 (off-grid): interpolated flag, between the two grid values.
    const double qe_b35 = attrs[2 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY];
    assert(statuses[2 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY] ==
           static_cast<int>(exaero::AttributeStatus::Interpolated));
    // log-space blend weight between bands 3 and 4 is 0.5358..., so the value sits just
    // above the arithmetic midpoint of the two grid values.
    const double qe_b3 = attrs[0 * NUM_ATTRIBUTES + EXTINCTION_EFFICIENCY];
    assert(qe_b35 > qe_b3 && qe_b35 > (qe_b3 + qe_b4) * 0.5);

    // Microphysical category: wavelength ignored, single block, no band dimension.
    double micro_attrs[9];
    for (double& x : micro_attrs) x = -777.0;
    exaero::View3D<double> micro_out(micro_attrs, num_cells, num_levels,
                                     exaero::microphysical_indices::NUM_ATTRIBUTES);
    package.computeAttributes(env, state, 1, exaero::AttributeCategory::Microphysical,
                              wavelengths, micro_out, nullptr);
    // rh=0.5 grid point: SS bin-0 effective radius is the grown value (not the dry golden).
    check_close(micro_attrs[exaero::microphysical_indices::EFFECTIVE_RADIUS],
                9.894295516232887e-08, "SS reff @ rh=0.5");

    std::cout << "MIE Multi-Band computeAttributes (FR-002/T027): PASS" << std::endl;
}

void test_mie_monochromatic_file() {
    using namespace exaero::spectral_optical_indices;
    // Runtime-file path (T029/FR-012/C7): a portable file replaces the 30-band RRTMG
    // lambda axis with a 3-wavelength monochromatic set and adds a file-only category
    // (polarized moments). The merged curve must surface the FILE values (delivery=file)
    // at the new wavelengths and report the old band coordinate as NotInSource (the axis
    // no longer contains it) - proving file values EXTEND/override baked, never silently.
    const std::string file = std::string(EXAERO_TEST_DATA_DIR) + "/mie_dust_monochromatic.txt";
    std::string yaml = R"YAML(species:
      - name: "dust_mono"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 0.15e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 0.1e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
        mie_table: {source: DU}
activation:
  data_file: )YAML" + file + "\n";

    exaero::GocartPackage package;
    package.initialize(yaml);

    double v = 0.0;
    exaero::ProvenanceInfo p{};
    // Monochromatic 0.55 um = the reused RRTMG band-3 value (DU qext[0,10,2]).
    auto st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                     EXTINCTION_EFFICIENCY, 0.5, 0.55e-6, &v, &p);
    assert(st == exaero::AttributeStatus::AvailableFile);
    check_close(v, 1.9339340925216675, "mono 0.55um qext");
    assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::RuntimeFile));
    assert(p.num_lambda == 3); // the file's monochromatic axis length is data (R10)

    // 0.355 um = band-3 * 1.30 (shorter wavelength -> larger extinction, monotone).
    st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                EXTINCTION_EFFICIENCY, 0.5, 0.355e-6, &v, &p);
    assert(st == exaero::AttributeStatus::AvailableFile);
    check_close(v, 2.5141143202781677, "mono 0.355um qext");

    // 1.33 um = band-3 * 0.60.
    st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                EXTINCTION_EFFICIENCY, 0.5, 1.33e-6, &v, &p);
    assert(st == exaero::AttributeStatus::AvailableFile);
    check_close(v, 1.1603604555130005, "mono 1.33um qext");

    // The old RRTMG band coordinate 3.0 is far above the new axis max (1.33e-6): the
    // documented above-max clamp returns the last point's value, never extrapolated
    // garbage and never a silent 0 (FR-007). A NaN wavelength is rejected outright.
    st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                EXTINCTION_EFFICIENCY, 0.5, 3.0, &v, nullptr);
    assert(st == exaero::AttributeStatus::AvailableFile);
    check_close(v, 1.1603604555130005, "above-max clamps to 1.33um edge");
    double sentinel = -777.0;
    st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                EXTINCTION_EFFICIENCY, 0.5, std::nan(""), &sentinel, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);
    assert(sentinel == -777.0); // NaN wavelength rejected, value untouched (FR-008)

    // Microphysical fields survive the spectral-axis replacement (radius/rh inherited).
    st = package.queryAttribute(0, exaero::AttributeCategory::Microphysical,
                                exaero::microphysical_indices::EFFECTIVE_RADIUS,
                                0.0, 0.0, &v, &p);
    assert(st == exaero::AttributeStatus::Available);
    check_close(v, 6.358845325848961e-07, "mono DU reff survives");

    std::cout << "MIE Monochromatic Runtime File (FR-012/C7/T029): PASS" << std::endl;
}

void test_mie_moments_not_available() {
    // Polarized phase-function moments (FR-003/FR-008/FR-013, T030): the DU runtime
    // file carries a rank-5 pmom array over (radius, rh, lambda, pol, moment). The
    // attribute_index encodes the (element, moment) pair as idx = moment*ELEMENT_STRIDE
    // + element with the documented ordering P11,P12,P33,P34,P22,P44 (data-model).
    // The fixture value for element e, moment m is qe00 * 0.5^m * (1 + 0.1*e) with
    // qe00 = DU qext[bin 0, rh 0.5, band 3] - all float64-exact.
    using namespace exaero::polarized_moment_indices;
    const std::string file = std::string(EXAERO_TEST_DATA_DIR) + "/mie_dust_monochromatic.txt";
    std::string yaml = R"YAML(species:
      - name: "dust_mono"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 0.15e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 0.1e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
        mie_table: {source: DU}
      - name: "SS"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 0.15e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 0.1e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
activation:
  categories: [microphysical, spectral, polarized]
  data_file: )YAML" + file + "\n";

    exaero::GocartPackage package;
    package.initialize(yaml);

    constexpr double qe00 = 1.9339340925216675; // DU qext at (bin 0, rh 0.5, band 3)
    constexpr double band = 5.5e-7;             // the file's central wavelength (grid point)
    double v;
    exaero::ProvenanceInfo p{};

    // (1) Count + ordering + values: all 6 elements x 3 moments at the grid point.
    for (int m = 0; m < 3; ++m) {
        for (int e = 0; e < ELEMENT_STRIDE; ++e) {
            v = -777.0;
            const int idx = m * ELEMENT_STRIDE + e;
            auto st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                             idx, 0.5, band, &v, &p);
            assert(st == exaero::AttributeStatus::AvailableFile);
            check_close(v, qe00 * std::pow(0.5, m) * (1.0 + 0.1 * e), "pmom value");
            assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::RuntimeFile));
            assert(p.interpolated == 0);
            // The moment-0 sequence must follow the documented element ordering:
            // values strictly increase with e (the fixture encodes e in the magnitude).
            if (m == 0 && e > 0) assert(v > 0.0);
        }
    }

    // (2) Out-of-range (element, moment) decomposition: moment index 3 (idx 18) does
    // not exist in this file (M = 3) -> explicit NotInSource, value untouched (FR-008).
    v = -777.0;
    auto st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                     3 * ELEMENT_STRIDE + 0, 0.5, band, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);
    assert(v == -777.0);
    st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                -1, 0.5, band, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);

    // (3) Off-grid RH interpolates between moment tables (FR-006). DU is RH-invariant
    // in qext, so the value is unchanged but the provenance must flag interpolation.
    v = -777.0;
    st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                PHASE_FUNCTION_MOMENT, 0.52, band, &v, &p);
    assert(st == exaero::AttributeStatus::Interpolated);
    assert(p.interpolated == 1);
    check_close(v, qe00, "pmom RH-interpolated (DU is RH-invariant)");

    // (4) NaN RH rejected, value untouched (FR-007 never silent 0).
    v = -777.0;
    st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                PHASE_FUNCTION_MOMENT, std::nan(""), band, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);
    assert(v == -777.0);

    // (5) A spherical species with no moments: explicit NotInSource (FR-008, scenario 2).
    v = -777.0;
    st = package.queryAttribute(1, exaero::AttributeCategory::PolarizedMoment,
                                PHASE_FUNCTION_MOMENT, 0.5, band, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);
    assert(v == -777.0);

    // (6) Category deselected -> explicit NotActivated even though data exists (FR-010).
    package.setAttributeActivation({},
        exaero::attribute_category_bit(exaero::AttributeCategory::Microphysical) |
        exaero::attribute_category_bit(exaero::AttributeCategory::SpectralOptical), "");
    st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                PHASE_FUNCTION_MOMENT, 0.5, band, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotActivated);

    std::cout << "MIE Polarized Moments (FR-003/FR-008/T030): PASS" << std::endl;
}

// Build a DU-bound package whose activation loads `file` as the runtime data file.
static std::string mie_file_yaml(const std::string& file) {
    return std::string(R"YAML(species:
      - name: "DU"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
        mie_table: {source: DU}
activation:
  data_file: )YAML") + file + "\n";
}

void test_mie_file_override_and_failfast() {
    using namespace exaero::spectral_optical_indices;
    const std::string dir(EXAERO_TEST_DATA_DIR);

    // (1) Valid override (C7/FR-012): the DU override file redefines bext on the SAME
    //     30-band axis with a distinct constant. The query must report the FILE value
    //     with delivery=file, proving file overrides baked for the same key.
    {
        exaero::GocartPackage package;
        package.initialize(mie_file_yaml(dir + "/mie_dust_override.txt"));
        double v = 0.0;
        exaero::ProvenanceInfo p{};
        auto st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                         MASS_EXTINCTION, 0.5, 3.0, &v, &p);
        assert(st == exaero::AttributeStatus::AvailableFile);
        check_close(v, 1234.5, "override bext");
        assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::RuntimeFile));
    }

    // (2) Fail-fast (C6/FR-009): each corrupt file must abort initialize() with a
    //     "FATAL ERROR:" diagnostic and NO silent fallback to baked-in data.
    struct CorruptCase { const char* file; const char* why; };
    const CorruptCase cases[] = {
        {"mie_corrupt_schema.txt",  "unknown directive"},
        {"mie_corrupt_unit.txt",    "unit mismatch"},
        {"mie_corrupt_version.txt", "missing source_version"},
        {"mie_truncated.txt",       "truncated field"},
    };
    for (const auto& c : cases) {
        exaero::GocartPackage package;
        bool threw = false;
        std::string msg;
        try {
            package.initialize(mie_file_yaml(dir + "/" + c.file));
        } catch (const std::exception& e) {
            threw = true;
            msg = e.what();
        }
        assert(threw); // must abort, never silently fall back
        assert(msg.find("FATAL ERROR:") != std::string::npos); // clear diagnostic (FR-009)
    }

    std::cout << "MIE File Override + Fail-Fast (FR-009/FR-012/C6/C7/T036): PASS" << std::endl;
}

void test_mie_default_set() {
    using namespace exaero::microphysical_indices;
    using namespace exaero::spectral_optical_indices;
    // Invariant C9 / FR-017: with NO runtime file, the default activation is exactly the
    // pinned baked-in set -- every bound species microphysical + RRTMG-band spectral --
    // while the monochromatic spectral axis and polarized moments are file-only. The
    // pinned release ships RRTMG band tables for six canonical species (BR has no band
    // file in the snapshot, so it is intentionally absent and reports NotInSource).
    exaero::GocartPackage package;
    package.initialize(kMieSixYaml);

    auto& store = exaero::MieTableStore::instance();
    // (1) Exactly the pinned baked set is bound (data, not a literal in the query path).
    assert(store.num_bound_species() == 6);

    // (2) Every bound species answers microphysical + spectral at grid points.
    double v = 0.0;
    exaero::ProvenanceInfo p{};
    for (int i = 0; i < 6; ++i) {
        auto st = package.queryAttribute(i, exaero::AttributeCategory::Microphysical,
                                         EFFECTIVE_RADIUS, 0.5, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        assert(v > 0.0 && std::isfinite(v));
        st = package.queryAttribute(i, exaero::AttributeCategory::SpectralOptical,
                                    EXTINCTION_EFFICIENCY, 0.5, 3.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        assert(p.num_lambda == 30); // the full RRTMG band set is active by default
    }

    // (3) Polarized moments are OFF in the default mask (FR-010/FR-017): NotActivated.
    auto st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                     exaero::polarized_moment_indices::PHASE_FUNCTION_MOMENT,
                                     0.5, 3.0, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotActivated);

    // (4) Even with the category force-activated, the baked curves carry no pmom field:
    //     explicit NotInSource -- moments require the runtime file (FR-008/C9).
    package.setAttributeActivation({},
        exaero::attribute_category_bit(exaero::AttributeCategory::Microphysical) |
        exaero::attribute_category_bit(exaero::AttributeCategory::SpectralOptical) |
        exaero::attribute_category_bit(exaero::AttributeCategory::PolarizedMoment), "");
    st = package.queryAttribute(0, exaero::AttributeCategory::PolarizedMoment,
                                exaero::polarized_moment_indices::PHASE_FUNCTION_MOMENT,
                                0.5, 3.0, &v, nullptr);
    assert(st == exaero::AttributeStatus::NotInSource);
    int n_pol = -1, n_mom = -1;
    package.momentCounts(0, &n_pol, &n_mom);
    assert(n_pol == 0 && n_mom == 0); // no moment data => no slots (FR-003 data-driven)

    std::cout << "MIE Default Baked Set (FR-017/C9/T037): PASS" << std::endl;
}

void test_mie_config_curves() {
    using namespace exaero::microphysical_indices;
    using namespace exaero::spectral_optical_indices;
    auto& store = exaero::MieTableStore::instance();
    const std::string dir(EXAERO_TEST_DATA_DIR);

    // ---- C11 config-invariance (research R9/R10, ADR-003) ---------------------------
    // A config species bound to DU that declares its OWN radius nodes -- the five source
    // nodes PLUS one inserted node (1.8e-6, between source bins 1 and 2) -- must reproduce
    // the baked curve bit-for-bit AT every source node (resample weight is exactly 0 there),
    // while still carrying the extra node. Delivery stays baked (resample is not an
    // override), so a source-node query reports Available, never Config/file.
    {
        std::string yaml = R"YAML(species:
      - name: "dust_cfg"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
        mie_table:
          source: DU
          radius_nodes: [6.358845325848961e-07, 1.324423010373721e-06, 1.8e-06,
                         2.301213726241258e-06, 4.1672033148643095e-06, 7.670712875551544e-06]
)YAML";
        exaero::GocartPackage package;
        package.initialize(yaml);

        // The inserted node sits at new index 2, so the five source nodes map to new
        // indices {0,1,3,4,5} (a data-derived remap, not a hardcoded count).
        const int src_to_new[5] = {0, 1, 3, 4, 5};
        for (int b = 0; b < 5; ++b) {
            const int nb = src_to_new[b];
            // Microphysical (rank-2 radius,rh) at the dry-ish grid RH 0.5.
            double cfg_v = 0.0, du_v = -1.0;
            auto s1 = store.query("dust_cfg", exaero::AttributeCategory::Microphysical,
                                  EFFECTIVE_RADIUS, 0.5, 0.0, nb, &cfg_v, nullptr);
            auto s0 = store.query("DU", exaero::AttributeCategory::Microphysical,
                                  EFFECTIVE_RADIUS, 0.5, 0.0, b, &du_v, nullptr);
            assert(s1 == exaero::AttributeStatus::Available); // still baked delivery
            assert(s0 == exaero::AttributeStatus::Available);
            check_close(cfg_v, du_v, "C11 reff source-node invariance");
            // Spectral (rank-3 radius,rh,lambda) at the band-3 grid point.
            double cfg_s = 0.0, du_s = -1.0;
            s1 = store.query("dust_cfg", exaero::AttributeCategory::SpectralOptical,
                             EXTINCTION_EFFICIENCY, 0.5, 3.0, nb, &cfg_s, nullptr);
            s0 = store.query("DU", exaero::AttributeCategory::SpectralOptical,
                             EXTINCTION_EFFICIENCY, 0.5, 3.0, b, &du_s, nullptr);
            assert(s1 == exaero::AttributeStatus::Available);
            check_close(cfg_s, du_s, "C11 qext source-node invariance");
        }
        // The public surface at the default solver node (0) is likewise unchanged.
        double v = 0.0;
        exaero::ProvenanceInfo p{};
        auto st = package.queryAttribute(0, exaero::AttributeCategory::Microphysical,
                                         EFFECTIVE_RADIUS, 0.0, 0.0, &v, &p);
        assert(st == exaero::AttributeStatus::Available);
        check_close(v, 6.358845325848961e-07, "C11 queryAttribute bin0 baked");
        assert(p.num_radius == 6); // the curve really carries the extra node (data, R10)
    }

    // ---- C12 no-hardcoding gate (invariant C12, R10) --------------------------------
    // A file-only species "DUST7" (absent from the baked set) declares a PRIME radius
    // count (7) and a NON-SOURCE RH length (7) over 3 bands -- none of these match the
    // baked DU 5x36x30 layout. Provenance extents must report the file's own counts and
    // every query path must resolve them, proving no bin/RH/band count is baked in.
    {
        std::string yaml = std::string(R"YAML(species:
      - name: "d7"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 1.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
        mie_table: {source: DUST7, solver_radius_node: 3}
activation:
  data_file: )YAML") + dir + "/mie_dust7.txt\n";
        exaero::GocartPackage package;
        package.initialize(yaml); // the 7x7x3 pool wiring + fail-fast length checks pass

        double v = 0.0;
        exaero::ProvenanceInfo p{};
        // rh=0.495 (index 3) and band 2.0 (index 1) are exact grid points; radius node 3.
        // qext[b,h,l] = (b+1) + 0.1*h + 0.01*l  =>  qext[3,3,1] = 4 + 0.3 + 0.01 = 4.31.
        auto st = package.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                         EXTINCTION_EFFICIENCY, 0.495, 2.0, &v, &p);
        assert(st == exaero::AttributeStatus::AvailableFile);
        check_close(v, 4.31, "C12 DUST7 qext grid point");
        // The reported extents are the FILE's prime counts, not any baked literal.
        assert(p.num_radius == 7);
        assert(p.num_rh == 7);
        assert(p.num_lambda == 3);
        assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::RuntimeFile));

        // The bulk path resolves the same curve across all three bands (device-block
        // consumers read the identical pool the store produced).
        const int num_cells = 1, num_levels = 1;
        double temp[1] = {298.0}, pres[1] = {101325.0}, dens[1] = {1.2};
        double rh_raw[1] = {0.495}, thick[1] = {100.0}, state_raw[1] = {1.0e-6};
        exaero::View2D<const double> temperature(temp, num_cells, num_levels);
        exaero::View2D<const double> pressure(pres, num_cells, num_levels);
        exaero::View2D<const double> air_density(dens, num_cells, num_levels);
        exaero::View2D<const double> relative_humidity(rh_raw, num_cells, num_levels);
        exaero::View2D<const double> layer_thickness(thick, num_cells, num_levels);
        exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                           relative_humidity, layer_thickness};
        exaero::View3D<const double> state(state_raw, num_cells, num_levels, 1);
        double wl_raw[3] = {1.0, 2.0, 3.0};
        exaero::View1D<const double> wavelengths(wl_raw, 3);
        constexpr int kSpecAttr = exaero::spectral_optical_indices::NUM_ATTRIBUTES;
        const int slots = 3 * kSpecAttr;
        double attrs[3 * 11];
        for (double& x : attrs) x = -777.0;
        exaero::View3D<double> attributes_out(attrs, num_cells, num_levels, slots);
        package.computeAttributes(env, state, 0, exaero::AttributeCategory::SpectralOptical,
                                  wavelengths, attributes_out, nullptr);
        // Band 2 (index 1) at radius node 3, rh index 3: qext = 4.31 (grid point).
        check_close(attrs[1 * kSpecAttr + EXTINCTION_EFFICIENCY], 4.31,
                    "C12 computeAttributes band2 qext");
    }

    // ---- C13 override precedence config > file > baked (research R9, FR-012) ---------
    // The SAME attribute (DU bext) is set at all three layers; the winning layer must be
    // reported by provenance.delivery_source at each stage.
    {
        // Stage 1 -- baked only: DU mass-extinction at (bin0, rh0.5, band3).
        double baked_v = 0.0;
        {
            exaero::GocartPackage pkg;
            pkg.initialize(kMieSixYaml);
            exaero::ProvenanceInfo p{};
            auto st = pkg.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                         MASS_EXTINCTION, 0.5, 3.0, &baked_v, &p);
            assert(st == exaero::AttributeStatus::Available);
            assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::BakedIn));
        }
        // Stage 2 -- file overrides baked: runtime file redefines bext = 1234.5.
        {
            exaero::GocartPackage pkg;
            pkg.initialize(mie_file_yaml(dir + "/mie_dust_override.txt"));
            double v = 0.0;
            exaero::ProvenanceInfo p{};
            auto st = pkg.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                         MASS_EXTINCTION, 0.5, 3.0, &v, &p);
            assert(st == exaero::AttributeStatus::AvailableFile);
            assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::RuntimeFile));
            check_close(v, 1234.5, "C13 file beats baked");
            assert(v != baked_v); // file value really differs from the baked golden
        }
        // Stage 3 -- config overrides file: the SAME file is loaded, then a config
        // override (bext = 999 at all 5 radius nodes) wins over the file's 1234.5.
        {
            std::string yaml = R"YAML(species:
      - name: "du_ov"
        dry_density: 2600.0
        molecular_weight: 100.0
        dry_particle_diameter: 2.0e-6
        hygroscopicity: 0.1
        lognormal_sigma: 1.5
        lognormal_dg: 1.0e-6
        refractive_index_real: 1.55
        refractive_index_imag: 0.002
        mie_table:
          source: DU
          overrides:
            bext: [999.0, 999.0, 999.0, 999.0, 999.0]
activation:
  data_file: __OVERRIDE_FILE__
)YAML";
            yaml = std::regex_replace(yaml, std::regex("__OVERRIDE_FILE__"),
                                      dir + "/mie_dust_override.txt");
            exaero::GocartPackage pkg;
            pkg.initialize(yaml);
            double v = 0.0;
            exaero::ProvenanceInfo p{};
            auto st = pkg.queryAttribute(0, exaero::AttributeCategory::SpectralOptical,
                                         MASS_EXTINCTION, 0.5, 3.0, &v, &p);
            assert(st == exaero::AttributeStatus::AvailableConfig);
            assert(p.delivery_source == static_cast<int>(exaero::DeliverySource::Config));
            check_close(v, 999.0, "C13 config beats file beats baked");
        }
    }

    std::cout << "MIE Config Curves: invariance + no-hardcoding + precedence "
                 "(R9/R10/C11/C12/C13/T044): PASS" << std::endl;
}

// --- T038 (SC-008 / invariant C10): coarse per-timestep overhead gate -----------------
// The curve-bound hot path must add <=1% to the per-timestep diagnostics+optics step
// relative to the identical analytical (ADT/Kohler) path, and the timestep loop must
// perform ZERO host->device transfers: the flat curve pool is uploaded exactly once in
// create_solver_state() (initialize()), and every step only wraps caller-owned host
// pointers in unmanaged Kokkos views (GocartSolver.cpp). A per-step H2D of the pool
// would show up both as a >1% A/B delta AND as a first-steps-vs-median spike, which the
// two assertions below detect. Timing is min-of-rounds (the standard noise-robust
// estimator) with A/B interleaved in the same process so machine drift hits both arms.
static double time_step_loop(exaero::GocartPackage& pkg,
                             const exaero::EnvironmentalStateView& env,
                             const exaero::View3D<const double>& state,
                             const exaero::View1D<const double>& wavelengths,
                             int num_cells, int num_levels, int steps,
                             std::vector<double>* step_ms) {
    using clock = std::chrono::steady_clock;
    std::vector<double> diags_raw(static_cast<std::size_t>(num_cells) * num_levels *
                                      exaero::diagnostic_indices::NUM_DIAGNOSTICS, 0.0);
    exaero::View3D<double> diagnostics_out(diags_raw.data(), num_cells, num_levels,
                                           exaero::diagnostic_indices::NUM_DIAGNOSTICS);
    std::vector<double> optics_raw(static_cast<std::size_t>(num_cells) * num_levels *
                                       wavelengths.extent(0) * exaero::optical_indices::NUM_OPTICS, 0.0);
    exaero::View4D<double> optics_out(optics_raw.data(), num_cells, num_levels,
                                      static_cast<int>(wavelengths.extent(0)),
                                      exaero::optical_indices::NUM_OPTICS);
    const auto t0 = clock::now();
    for (int s = 0; s < steps; ++s) {
        const auto ts0 = clock::now();
        pkg.computeDerivedDiagnostics(env, state, diagnostics_out);
        pkg.computeOptics(env, state, wavelengths, optics_out);
        exaero::fence_environment();
        if (step_ms) {
            const auto ts1 = clock::now();
            step_ms->push_back(std::chrono::duration<double, std::milli>(ts1 - ts0).count());
        }
    }
    const auto t1 = clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void test_mie_hot_path_overhead(bool verbose) {
    constexpr int kCells = 64, kLevels = 8, kSteps = 200, kRounds = 7;
    // Identical single-species grids; the ONLY difference is the curve binding, so the
    // A/B delta isolates the table-read overhead (same launch counts, same memory).
    auto make_yaml = [](const char* name, bool bound) {
        std::string y = std::string("species:\n  - name: \"") + name + R"("
    dry_density: 2600.0
    molecular_weight: 100.0
    dry_particle_diameter: 2.0e-6
    hygroscopicity: 0.1
    lognormal_sigma: 1.5
    lognormal_dg: 1.0e-6
    refractive_index_real: 1.55
    refractive_index_imag: 0.002
)";
        if (bound) y += "    mie_table: {source: DU}\n";
        return y;
    };
    exaero::GocartPackage adt_pkg, curve_pkg;
    adt_pkg.initialize(make_yaml("dust_adt", false));   // analytical path (no curve)
    curve_pkg.initialize(make_yaml("dust_curve", true)); // device table path

    // Shared host inputs (unmanaged views: neither package owns or copies them).
    std::vector<double> temp(static_cast<std::size_t>(kCells) * kLevels, 298.0);
    std::vector<double> pres(static_cast<std::size_t>(kCells) * kLevels, 101325.0);
    std::vector<double> dens(static_cast<std::size_t>(kCells) * kLevels, 1.2);
    std::vector<double> rh(static_cast<std::size_t>(kCells) * kLevels);
    std::vector<double> thick(static_cast<std::size_t>(kCells) * kLevels, 100.0);
    for (int i = 0; i < kCells * kLevels; ++i) rh[i] = 0.3 + 0.4 * (i % 37) / 36.0;
    std::vector<double> state_raw(static_cast<std::size_t>(kCells) * kLevels, 1.0e-6);
    exaero::View2D<const double> temperature(temp.data(), kCells, kLevels);
    exaero::View2D<const double> pressure(pres.data(), kCells, kLevels);
    exaero::View2D<const double> air_density(dens.data(), kCells, kLevels);
    exaero::View2D<const double> relative_humidity(rh.data(), kCells, kLevels);
    exaero::View2D<const double> layer_thickness(thick.data(), kCells, kLevels);
    exaero::EnvironmentalStateView env{temperature, pressure, air_density,
                                       relative_humidity, layer_thickness};
    exaero::View3D<const double> state(state_raw.data(), kCells, kLevels, 1);
    // Band coordinates INSIDE the DU curve domain (indices 1..30): the table read is
    // ACTIVE (replaces the ADT series), which is the default-activation hot path FR-015.
    double wl_raw[3] = {3.0, 4.0, 5.0};
    exaero::View1D<const double> wavelengths(wl_raw, 3);

    // Warmup both arms (JIT, page faults, OpenMP thread spin-up).
    time_step_loop(adt_pkg, env, state, wavelengths, kCells, kLevels, 25, nullptr);
    time_step_loop(curve_pkg, env, state, wavelengths, kCells, kLevels, 25, nullptr);

    double best_adt = 1e300, best_curve = 1e300;
    std::vector<double> curve_steps;
    for (int r = 0; r < kRounds; ++r) {
        const double a = time_step_loop(adt_pkg, env, state, wavelengths, kCells, kLevels,
                                        kSteps, nullptr);
        const bool collect = (r == kRounds - 1);
        const double c = time_step_loop(curve_pkg, env, state, wavelengths, kCells, kLevels,
                                        kSteps, collect ? &curve_steps : nullptr);
        best_adt = std::min(best_adt, a);
        best_curve = std::min(best_curve, c);
    }

    // (1) SC-008/C10: curve path within 1% of the analytical baseline (min-of-rounds).
    // The strict gate is meaningful in OPTIMIZED builds (the spec's performance target is
    // the production step; measurements here: ratio ~0.78 Release). In a Debug build the
    // un-inlined ADT Mie series and the table path's binary search have distorted relative
    // costs (measured ~1.06), so Debug enforces a coarse structural bound instead; the
    // zero-transfer proxy (2) is build-independent. assert() would compile out under
    // NDEBUG, so this gate throws explicitly in both build types.
    if (verbose) {
        std::cout << "  overhead bench: adt=" << best_adt << " ms / curve=" << best_curve
                  << " ms over " << kSteps << " steps (ratio "
                  << best_curve / best_adt << ")" << std::endl;
    }
#ifdef NDEBUG
    constexpr double kRatioLimit = 1.01; // SC-008 strict (Release)
#else
    constexpr double kRatioLimit = 1.25; // Debug: order-of-regression guard only
#endif
    if (!(best_curve <= best_adt * kRatioLimit)) {
        throw std::runtime_error("FATAL ERROR: SC-008 overhead gate: curve path " +
            std::to_string(best_curve) + " ms vs analytical " + std::to_string(best_adt) +
            " ms (ratio " + std::to_string(best_curve / best_adt) + " > " +
            std::to_string(kRatioLimit) + ")");
    }

    // (2) Zero-transfer proxy: a hidden per-step (or first-step lazy) H2D would inflate
    //     the opening steps far above the median. Median-of-steps vs mean-of-first-20.
    std::vector<double> sorted = curve_steps;
    std::sort(sorted.begin(), sorted.end());
    const double median = sorted[sorted.size() / 2];
    double head = 0.0;
    for (int i = 0; i < 20; ++i) head += curve_steps[i];
    head /= 20.0;
    if (verbose) {
        std::cout << "  step flatness: median=" << median << " ms, first20 mean=" << head
                  << " ms" << std::endl;
    }
    if (!(head <= std::max(median * 3.0, median + 0.05))) {
        throw std::runtime_error("FATAL ERROR: SC-008 zero-transfer proxy: first-20-step mean " +
            std::to_string(head) + " ms vs median " + std::to_string(median) +
            " ms (per-step transfer suspected)");
    }

    std::cout << "MIE Hot-Path Overhead (SC-008/C10/T038): PASS" << std::endl;
}

// A package that does not override the structured overload must throw
// std::logic_error (additive-surface pattern, contract §2).
struct UnsupportedPkg : public exaero::IAerosolPackage {
    using IAerosolPackage::initialize; // un-hide the structured overload for concrete-type calls
    void initialize(const std::string&) override {}
    void executeMicrophysics(const exaero::EnvironmentalStateView&,
                             exaero::View3D<double>&, double) override {}
    void computeDerivedDiagnostics(const exaero::EnvironmentalStateView&,
                                   const exaero::View3D<const double>&,
                                   exaero::View3D<double>&) override {}
    void computeEmissions(const exaero::EnvironmentalStateView&,
                          const exaero::EmissionsInputView&,
                          exaero::View3D<double>&) override {}
    void computeOptics(const exaero::EnvironmentalStateView&,
                       const exaero::View3D<const double>&,
                       const exaero::View1D<const double>&,
                       exaero::View4D<double>&) override {}
    void computeCCN(const exaero::EnvironmentalStateView&,
                    const exaero::View3D<const double>&,
                    const exaero::View1D<const double>&,
                    exaero::View4D<double>&) override {}
    int getSpeciesIndex(const std::string&) const override { return -1; }
    std::string getSpeciesName(int) const override { return {}; }
};

// Explicit throw-based gate: bare assert() compiles out under NDEBUG and would
// make Release runs a no-op (T038 precedent). ALL structured-config tests use it.
static void gate(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(std::string("FATAL ERROR: structured-config test: ") + what);
}

void test_structured_config_unsupported_package() {
    UnsupportedPkg pkg;
    exaero::GocartConfig cfg;
    cfg.species.push_back(exaero::GocartSpeciesConfig{});
    bool threw = false;
    try { pkg.initialize(cfg); }
    catch (const std::logic_error& e) {
        threw = std::string(e.what()).find("does not support structured configuration") != std::string::npos;
    }
    gate(threw, "default structured initialize must throw logic_error with the contract message");
    std::cout << "Structured config unsupported-package guard: PASS" << std::endl;
}

// Build the GocartConfig equivalent of kMieSixYaml (six species, no optional blocks).
static exaero::GocartConfig kMieSixConfig() {
    auto sp = [](const char* n, double dens, double mw, double dpg, double kap,
                 double sig, double dg, double nr, double ni) {
        exaero::GocartSpeciesConfig c;
        c.name = n; c.dry_density = dens; c.molecular_weight = mw;
        c.dry_particle_diameter = dpg; c.hygroscopicity = kap;
        c.lognormal_sigma = sig; c.lognormal_dg = dg;
        c.refractive_index_real = nr; c.refractive_index_imag = ni;
        return c;
    };
    exaero::GocartConfig cfg;
    cfg.species = {
        sp("DU", 2600.0, 100.0, 2.0e-6, 0.1, 1.5, 1.0e-6, 1.55, 0.002),
        sp("SS", 1800.0, 98.0, 0.2e-6, 0.5, 2.0, 0.15e-6, 1.43, 1.0e-8),
        sp("SU", 1800.0, 98.0, 0.2e-6, 0.5, 2.0, 0.15e-6, 1.43, 1.0e-8),
        sp("BC", 1800.0, 12.0, 0.1e-6, 0.0, 1.8, 0.05e-6, 1.85, 0.75),
        sp("OC", 1300.0, 150.0, 0.1e-6, 0.1, 1.8, 0.05e-6, 1.55, 0.0),
        sp("NI", 2150.0, 85.0, 0.3e-6, 0.6, 1.6, 0.1e-6, 1.52, 0.01),
    };
    return cfg;
}

// Task 3 (design 2026-09-10): the structured entry point must produce a package
// state bit-identical to the YAML path for the same configuration — same scalars,
// same lookup/curve/micro/spectral/pmom pool wiring. Throw-based gates (not assert)
// so Release builds enforce parity too.
void test_structured_config_species_parity() {
    exaero::GocartPackage pkg_yaml, pkg_cfg;
    pkg_yaml.initialize(kMieSixYaml);
    pkg_cfg.initialize(kMieSixConfig());

    gate(pkg_cfg.get_num_species() == 6, "six species configured");
    for (int i = 0; i < 6; ++i) {
        auto a = pkg_yaml.get_species_params(i);
        auto b = pkg_cfg.get_species_params(i);
        gate(pkg_cfg.getSpeciesName(i) == pkg_yaml.getSpeciesName(i), "species name parity");
        // Bit-identical scalars and pool wiring (offsets/extents included).
        gate(b.dry_density == a.dry_density && b.molecular_weight == a.molecular_weight,
             "density/molecular-weight parity");
        gate(b.dry_particle_diameter == a.dry_particle_diameter, "dpg parity");
        gate(b.hygroscopicity == a.hygroscopicity && b.lognormal_sigma == a.lognormal_sigma,
             "kappa/sigma parity");
        gate(b.lognormal_dg == a.lognormal_dg, "dg parity");
        gate(b.refractive_index_real == a.refractive_index_real, "n parity");
        gate(b.refractive_index_imag == a.refractive_index_imag, "k parity");
        gate(b.has_optics_lookup == a.has_optics_lookup, "lookup flag parity");
        gate(b.curve_offset == a.curve_offset && b.n_rh == a.n_rh, "legacy curve wiring parity");
        gate(b.micro_offset == a.micro_offset && b.spec_offset == a.spec_offset,
             "micro/spectral wiring parity");
        gate(b.pmom_offset == a.pmom_offset, "pmom wiring parity");
    }
    std::cout << "Structured config species parity: PASS" << std::endl;
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

void test_gocart_emissions_mapping() {
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
      - name: "Sulfate_1"
        dry_density: 1800.0
        molecular_weight: 98.0
        dry_particle_diameter: 0.20e-6
        hygroscopicity: 0.50
        lognormal_sigma: 1.6
        lognormal_dg: 0.20e-6
        refractive_index_real: 1.43
        refractive_index_imag: 1.0e-7
    emissions_mapping:
      - raw_name: "CECE_Dust"
        mappings:
          - target_species: "Dust_1"
            mass_split_fraction: 0.40
      - raw_name: "CECE_Sulfate"
        mappings:
          - target_species: "Sulfate_1"
            mass_split_fraction: 0.60
            is_modal_mode: true
            emitted_particle_diameter: 0.25e-6
            lognormal_sigma: 1.8
    )";

    exaero::GocartPackage package;
    package.initialize(yaml_string);

    // Inputs (1 cell, 1 level)
    double temp_raw[1] = { 298.0 };
    double pres_raw[1] = { 101325.0 };
    double dens_raw[1] = { 1.2 };
    double rh_raw[1] = { 0.50 };
    double thick_raw[1] = { 50.0 }; // Δz = 50 m

    exaero::View2D<const double> temperature(temp_raw, 1, 1);
    exaero::View2D<const double> pressure(pres_raw, 1, 1);
    exaero::View2D<const double> air_density(dens_raw, 1, 1);
    exaero::View2D<const double> relative_humidity(rh_raw, 1, 1);
    exaero::View2D<const double> layer_thickness(thick_raw, 1, 1);

    exaero::EnvironmentalStateView env{temperature, pressure, air_density, relative_humidity, layer_thickness};

    // Raw fluxes: 2 CECE species. 
    // CECE_Dust = 1.0e-4 kg/m²/s (AREA_FLUX)
    // CECE_Sulfate = 1.0e-4 kg/m³/s (MASS_CONCENTRATION_RATE)
    double flux_raw[2] = { 1.0e-4, 1.0e-4 };
    exaero::View3D<const double> flux(flux_raw, 1, 1, 2);

    // Outputs: size (1, 1, 2)
    double out_raw[2] = { 0.0, 0.0 };
    exaero::View3D<double> emissions_out(out_raw, 1, 1, 2);

    // 1. Run Area-Flux scaling test (CECE_Dust)
    exaero::EmissionsInputView em_in_area{flux, exaero::FluxType::AREA_FLUX};
    package.computeEmissions(env, em_in_area, emissions_out);

    // Expect: (1.0e-4 / 50.0) * 0.40 = 8.0e-7 kg/m³/s
    double mass_dust_rate = emissions_out(0, 0, 0);
    assert(std::abs(mass_dust_rate - 8.0e-7) < 1e-12);

    // 2. Run Mass-Concentration-Rate & Modal-Number conversion test (CECE_Sulfate)
    // Clear outputs and zero out CECE_Dust flux to isolate CECE_Sulfate
    flux_raw[0] = 0.0;
    out_raw[0] = 0.0; out_raw[1] = 0.0;
    exaero::EmissionsInputView em_in_rate{flux, exaero::FluxType::MASS_CONCENTRATION_RATE};
    package.computeEmissions(env, em_in_rate, emissions_out);

    // Expect: target mass rate = 1.0e-4 * 0.60 = 6.0e-5 kg/m³/s.
    // Vol factor = (M_PI/6) * 1800 * (0.25e-6)^3 * exp(4.5 * ln^2(1.8)) ≈ 6.97103e-17
    // Emitted Number rate = 6.0e-5 / 6.97103e-17 ≈ 8.60705e11
    double mass_sulfate_rate = emissions_out(0, 0, 0); // Dust is 0 because CECE_Sulfate has no mapping to Dust_1
    double num_sulfate_rate = emissions_out(0, 0, 1);
    
    assert(mass_sulfate_rate == 0.0);
    assert(std::abs(num_sulfate_rate - 8.60705e11) / 8.60705e11 < 1e-6);

    // 3. Boundary Clamping & Crash fuzzer test
    double nan_flux_raw[2] = { NAN, -5.0 }; // Pass NaN and negatives
    exaero::View3D<const double> nan_flux(nan_flux_raw, 1, 1, 2);
    exaero::EmissionsInputView em_nan{nan_flux, exaero::FluxType::AREA_FLUX};
    out_raw[0] = 1.234; out_raw[1] = 5.678; // non-zero default
    
    // Execute fuzzer-condition
    package.computeEmissions(env, em_nan, emissions_out);
    
    // NaNs and negatives must be clamped defensively to 0.0, yielding 0.0 outputs
    assert(emissions_out(0, 0, 0) == 0.0);
    assert(emissions_out(0, 0, 1) == 0.0);

    std::cout << "GOCART Parallel GPU Emissions Mapping Solver Tests: PASS" << std::endl;
}

void test_lut_generator_stub() {
    std::string yaml_string = R"(
    generation_grid:
      rh_bins: [0.0, 0.50, 0.99]
      wavelengths: [550.0e-9]
      legendre_moments: 16
    )";
    exaero::LutGenerator generator(yaml_string);
    assert(generator.num_rh() == 3);
    assert(generator.num_bands() == 1);
    assert(generator.num_moments() == 16);
    std::cout << "LutGenerator Stub test: PASS" << std::endl;
}

void test_spheroid_database_interpolation() {
    const auto& db = exaero::SpheroidDatabase::instance();

    // 1. In-bounds query: Exact match at grid point (1.53, 0.003, 1.0)
    auto p_exact = db.interpolate(1.53, 0.003, 1.0);
    double expected_ext = 2.0 * (1.0 - std::exp(-0.5)) * (1.53 / 1.53); // 2.0 * (1.0 - exp(-0.5)) ≈ 0.786938
    assert(std::abs(p_exact.ext_efficiency - expected_ext) < 1e-6);
    assert(p_exact.sca_efficiency == p_exact.ext_efficiency * (1.0 - 0.003 * 10.0));

    // 2. Linear Interpolation check: midpoint between 1.53 and 1.56 at (1.545, 0.003, 1.0)
    auto p_mid = db.interpolate(1.545, 0.003, 1.0);
    double expected_mid_ext = 2.0 * (1.0 - std::exp(-0.5)) * (1.545 / 1.53); // midpoint ext
    assert(std::abs(p_mid.ext_efficiency - expected_mid_ext) < 1e-6);

    // 3. Boundary Clamping check: pass values exceeding grid limits (e.g. n_real = 1.60, size = 15.0)
    // It must clamp n_real to 1.56, n_imag to 0.008, and size to 10.0 defensively
    auto p_clamp = db.interpolate(1.60, 0.010, 15.0);
    auto p_limit = db.interpolate(1.56, 0.008, 10.0);
    assert(p_clamp.ext_efficiency == p_limit.ext_efficiency);
    assert(p_clamp.sca_efficiency == p_limit.sca_efficiency);
    assert(p_clamp.moments[5] == p_limit.moments[5]);

    std::cout << "Dubovik Spheroid Database Trilinear Interpolation & Bounds Clamping: PASS" << std::endl;
}

int main(int argc, char** argv) {
    // The overhead bench is opt-in (quickstart §7: --bench-optics-defaults); the coarse
    // assertion still runs by default so CI keeps the SC-008 gate, with verbose output
    // when explicitly requested.
    const bool bench_verbose = (argc > 1 && std::string(argv[1]) == "--bench-optics-defaults");
    exaero::initialize_environment();
    
    test_gocart_yaml_parsing();
    test_gocart_yaml_emissions_parsing();
    test_gocart_emissions_mapping();
    test_lut_generator_stub();
    test_spheroid_database_interpolation();
    test_zero_copy_mapping();
    test_gocart_optics();
    test_gocart_arbitrary_lookup_length();
    test_gocart_lookup_length_mismatch_fails();
    test_mie_baked_in_microphysical();
    test_mie_provenance();
    test_mie_not_available();
    test_mie_rrtmg_spectral();
    test_mie_interpolation();
    test_mie_compute_attributes_multiband();
    test_mie_monochromatic_file();
    test_mie_moments_not_available();
    test_mie_file_override_and_failfast();
    test_mie_default_set();
    test_mie_config_curves();
    test_mie_hot_path_overhead(bench_verbose);
    test_structured_config_unsupported_package();
    test_structured_config_species_parity();
    test_optical_precision();
    test_gocart_ccn();
    
    exaero::finalize_environment();
    return 0;
}
