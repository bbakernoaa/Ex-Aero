#pragma once
#include <string>
#include <vector>
#include <exaero/AerosolIndices.hpp>
#include <exaero/AttributeQuery.hpp>

// Force backport to compile under exaero_mdspan namespace to prevent redefinition conflicts with system Kokkos
#define MDSPAN_IMPL_STANDARD_NAMESPACE exaero_mdspan
#include <experimental/mdspan>

namespace exaero {

    template <typename T>
    using View1D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent>>;

    // Transition 2D, 3D, and 4D layouts to column-major (layout_left) for zero-copy Fortran/CCPP memory alignments
    template <typename T>
    using View2D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent>, exaero_mdspan::layout_left>;

    template <typename T>
    using View3D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, exaero_mdspan::layout_left>;

    template <typename T>
    using View4D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>, exaero_mdspan::layout_left>;

    struct EnvironmentalStateView {
        View2D<const double> temperature;       // [K] (cell, level)
        View2D<const double> pressure;          // [Pa] (cell, level)
        View2D<const double> air_density;       // [kg/m³] (cell, level)
        View2D<const double> relative_humidity;  // [fraction 0-1] (cell, level)
        View2D<const double> layer_thickness;   // [m] grid layer vertical thickness Delta_z (cell, level)
    };

    enum class FluxType {
        MASS_CONCENTRATION_RATE,  // [kg/m³/s]
        AREA_FLUX                 // [kg/m²/s]
    };

    struct EmissionsInputView {
        View3D<const double> flux;
        FluxType flux_type;
    };

    class IAerosolPackage {
    public:
        virtual ~IAerosolPackage() = default;

        // Dynamic Initialization (Loads YAML config)
        virtual void initialize(const std::string& config_yaml) = 0;

        // Passive microphysics step
        virtual void executeMicrophysics(
            const EnvironmentalStateView& env,
            View3D<double>& state,          // (cell, level, species)
            double delta_time_sec) = 0;

        // Active Diagnostics Step (PM2.5, PM10, Derived Number, Column Mass, ALW, etc.)
        virtual void computeDerivedDiagnostics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& diagnostics_out) = 0; // (cell, level, diagnostic_index)

        // Active Optical Properties Step (Extinction, Scattering, Lidar Backscatter, AOT, etc. across multiple bands)
        virtual void computeOptics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            const View1D<const double>& wavelengths, // [m] queried wavelengths (band_index)
            View4D<double>& optics_out) = 0;      // (cell, level, band_index, optics_index)

        // Active Cloud Microphysics CCN Activation Spectra (liquid activation over supersaturations)
        virtual void computeCCN(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,
            const View1D<const double>& supersaturations, // [fraction 0-1] queried supersaturations (S_index)
            View4D<double>& ccn_out) = 0;                 // (cell, level, S_index, activated_ccn_number_concentration)

        // Maps raw incoming emissions to package-specific mass and number emission rates
        virtual void computeEmissions(
            const EnvironmentalStateView& env,
            const EmissionsInputView& emissions_in,
            View3D<double>& emissions_out) = 0;

        // Dynamic Species-to-Index mapping query APIs to avoid hardcoding on the host side
        virtual int getSpeciesIndex(const std::string& name) const = 0;
        virtual std::string getSpeciesName(int index) const = 0;

        // --- GEOSmie MIE table attribute surface (additive, backward compatible) ---
        // Defaults return NotActivated / do nothing so existing packages (MAM4xx wrapper,
        // tests) compile and behave unchanged (Principle II, contract §2).

        // Hot-path, device-resident bulk query mirroring computeOptics (FR-015).
        // attributes_out: (cell, level, attribute_index) for one species + category.
        // wavelengths is ignored for Microphysical. status_out optionally receives the
        // AttributeStatus code per (cell, level, attribute_index).
        virtual void computeAttributes(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,
            int species_index,
            AttributeCategory category,
            const View1D<const double>& wavelengths,
            View3D<double>& attributes_out,
            View3D<int>* status_out = nullptr) {
            (void)env; (void)state; (void)species_index; (void)category;
            (void)wavelengths; (void)attributes_out; (void)status_out;
        }

        // Scalar query for init-time / retrieval consumers (FR-001..FR-003).
        // wavelength_m is ignored (pass NaN) for Microphysical. Returns the availability
        // status; value_out/provenance_out are written only when the status is available
        // or interpolated (FR-008: never a silent 0).
        virtual AttributeStatus queryAttribute(
            int species_index,
            AttributeCategory category,
            int attribute_index,
            double rh,
            double wavelength_m,
            double* value_out,
            ProvenanceInfo* provenance_out = nullptr) const {
            (void)species_index; (void)category; (void)attribute_index;
            (void)rh; (void)wavelength_m; (void)value_out; (void)provenance_out;
            return AttributeStatus::NotActivated;
        }

        // Activation control (FR-010) with optional runtime extension/override file
        // (FR-012). A file that fails validation aborts initialization with a
        // "FATAL ERROR:" diagnostic -- no silent fallback (FR-009).
        virtual void setAttributeActivation(
            const std::vector<std::string>& species,
            int categories_mask,
            const std::string& runtime_file_path = "") {
            (void)species; (void)categories_mask; (void)runtime_file_path;
        }

        // Curve mapping for runtime-configurable species sets (CATChem; research R9,
        // ADR-003). Resolved during initialize(); bad config aborts with "FATAL ERROR:".
        virtual void setSpeciesCurveConfig(const std::vector<SpeciesCurveConfig>& curves) {
            (void)curves;
        }

        // Polarized-moment slot counts for one species (FR-003/FR-008): the number of
        // (element, moment) slots its curve carries, sized from data. Both outputs are 0
        // when the species ships no moments or the category is not activated (FR-010).
        // Defaults to zeros so packages without moments compile unchanged.
        virtual void momentCounts(int species_index, int* num_pol_out, int* num_moment_out) const {
            (void)species_index;
            if (num_pol_out) *num_pol_out = 0;
            if (num_moment_out) *num_moment_out = 0;
        }
    };

} // namespace exaero
