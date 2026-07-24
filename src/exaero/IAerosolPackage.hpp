#pragma once
#include <string>

// Force backport to compile under exaero_mdspan namespace to prevent redefinition conflicts with system Kokkos
#define MDSPAN_IMPL_STANDARD_NAMESPACE exaero_mdspan
#include <exaero_experimental/mdspan>

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

    namespace diagnostic_indices {
        // 3D Mass and Size Diagnostics (cell, level, index)
        constexpr int MASS_CONCENTRATION = 0;       // 3D mass concentration per species [kg/m³]
        constexpr int PM2_5_CONCENTRATION = 1;      // 3D PM2.5 mass concentration [kg/m³]
        constexpr int PM10_CONCENTRATION = 2;       // 3D PM10 mass concentration [kg/m³]
        constexpr int NUMBER_CONCENTRATION = 3;     // Derived number concentration [particles/m³]
        constexpr int SURFACE_AREA_DENSITY = 4;     // SAD [m²/m³]
        
        // 2D Column-Integrated and Surface Mass Diagnostics (cell, level=0, index)
        constexpr int SURFACE_MASS = 5;             // Surface layer mass concentration [kg/m³]
        constexpr int COLUMN_MASS = 6;              // Column-integrated mass density [kg/m²]
        constexpr int SURFACE_PM2_5_MASS = 7;       // Surface layer PM2.5 mass concentration [kg/m³]
        constexpr int COLUMN_PM2_5_MASS = 8;        // Column-integrated PM2.5 mass density [kg/m²]
        
        constexpr int NUM_DIAGNOSTICS = 9;
    }

    namespace optical_indices {
        // 3D Optical Properties (cell, level, index)
        constexpr int EXTINCTION_COEFF = 0;         // 3D extinction coefficient [m⁻¹]
        constexpr int SCATTERING_COEFF = 1;         // 3D scattering coefficient [m⁻¹]
        constexpr int BACKSCATTER_COEFF = 2;        // 3D backscatter coefficient [m⁻¹ sr⁻¹]
        constexpr int ASYMMETRY_FACTOR = 3;         // g asymmetry parameter [fraction -1 to 1]
        
        // 2D Column-Integrated Optical Diagnostics (Column AOT / AOD at 550nm) (cell, level=0, index)
        constexpr int EXTINCTION_AOT = 4;           // Total Extinction Aerosol Optical Thickness (AOT)
        constexpr int SCATTERING_AOT = 5;           // Total Scattering Aerosol Optical Thickness (AOT)
        constexpr int FINE_MODE_EXTINCTION_AOT = 6;  // Fine mode (sub-micron) extinction AOT
        constexpr int FINE_MODE_SCATTERING_AOT = 7;  // Fine mode (sub-micron) scattering AOT
        constexpr int PM2_5_EXTINCTION_AOT = 8;     // PM2.5 extinction AOT
        constexpr int PM2_5_SCATTERING_AOT = 9;     // PM2.5 scattering AOT
        
        constexpr int NUM_OPTICS = 10;
    }

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
