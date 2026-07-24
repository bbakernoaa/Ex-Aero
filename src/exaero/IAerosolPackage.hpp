#pragma once
#include <string>
#include <exaero/AerosolIndices.hpp>

// Force backport to compile under exaero_mdspan namespace to prevent redefinition conflicts with system Kokkos
#define MDSPAN_IMPL_STANDARD_NAMESPACE exaero_mdspan
#include <experimental/mdspan>

namespace exaero {

    template <typename T>
    using View3D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>>;

    template <typename T>
    using View2D = exaero_mdspan::mdspan<T, exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent>>;

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

        // Active Diagnostics Step (PM2.5, PM10, Derived Number, Column Mass, etc.)
        virtual void computeDerivedDiagnostics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& diagnostics_out) = 0; // (cell, level, diagnostic_index)

        // Active Optical Properties Step (Extinction, Scattering, AOT, etc.)
        virtual void computeOptics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& optics_out) = 0;      // (cell, level, optics_index)
    };

} // namespace exaero
