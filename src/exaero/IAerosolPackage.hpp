#pragma once
#include <string>
#include <exaero/AerosolIndices.hpp>

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

        // Dynamic Species-to-Index mapping query APIs to avoid hardcoding on the host side
        virtual int getSpeciesIndex(const std::string& name) const = 0;
        virtual std::string getSpeciesName(int index) const = 0;
    };

} // namespace exaero
