#pragma once
#include <string>
#include <experimental/mdspan>

namespace exaero {

    template <typename T>
    using View3D = std::experimental::mdspan<T, std::experimental::extents<size_t, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>>;

    template <typename T>
    using View2D = std::experimental::mdspan<T, std::experimental::extents<size_t, std::dynamic_extent, std::dynamic_extent>>;

    struct EnvironmentalStateView {
        View2D<const double> temperature;       // [K] (cell, level)
        View2D<const double> pressure;          // [Pa] (cell, level)
        View2D<const double> air_density;       // [kg/m³] (cell, level)
        View2D<const double> relative_humidity;  // [fraction 0-1] (cell, level)
    };

    namespace diagnostic_indices {
        constexpr int PM2_5 = 0;
        constexpr int PM10 = 1;
        constexpr int NUMBER_CONCENTRATION = 2; // Derived number concentration [particles/m³]
        constexpr int SURFACE_AREA_DENSITY = 3; // [m²/m³]
        constexpr int CCN = 4;                  // Cloud Condensation Nuclei
        constexpr int INP = 5;                  // Ice Nucleating Particles
        constexpr int NUM_DIAGNOSTICS = 6;
    }

    namespace optical_indices {
        constexpr int EXTINCTION = 0;               // Extinction coefficient [m⁻¹]
        constexpr int SINGLE_SCATTERING_ALBEDO = 1; // SSA [fraction 0-1]
        constexpr int ASYMMETRY_FACTOR = 2;        // g [fraction -1 to 1]
        constexpr int NUM_OPTICS = 3;
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

        // Active Diagnostics Step
        virtual void computeDerivedDiagnostics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& diagnostics_out) = 0; // (cell, level, diagnostic_index)

        // Active Optical Properties Step
        virtual void computeOptics(
            const EnvironmentalStateView& env,
            const View3D<const double>& state,    // (cell, level, species)
            View3D<double>& optics_out) = 0;      // (cell, level, optics_index)
    };

} // namespace exaero
