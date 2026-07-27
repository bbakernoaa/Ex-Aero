#pragma once

namespace exaero {

namespace diagnostic_indices {
// 3D Mass and Size Diagnostics (cell, level, index)
constexpr int MASS_CONCENTRATION =
    0; // 3D mass concentration per species [kg/m³]
constexpr int PM2_5_CONCENTRATION = 1; // 3D PM2.5 mass concentration [kg/m³]
constexpr int PM10_CONCENTRATION = 2;  // 3D PM10 mass concentration [kg/m³]
constexpr int NUMBER_CONCENTRATION =
    3; // Derived number concentration [particles/m³]
constexpr int SURFACE_AREA_DENSITY = 4;            // SAD [m²/m³]
constexpr int AEROSOL_LIQUID_WATER = 5;            // ALW [kg/m³]
constexpr int GRAVITATIONAL_SETTLING_VELOCITY = 6; // vg settling velocity [m/s]

// 2D Column-Integrated and Surface Mass Diagnostics (cell, level=0, index)
constexpr int SURFACE_MASS = 7; // Surface layer mass concentration [kg/m³]
constexpr int COLUMN_MASS = 8;  // Column-integrated mass density [kg/m²]
constexpr int SURFACE_PM2_5_MASS =
    9; // Surface layer PM2.5 mass concentration [kg/m³]
constexpr int COLUMN_PM2_5_MASS =
    10; // Column-integrated PM2.5 mass density [kg/m²]

constexpr int NUM_DIAGNOSTICS = 11;
} // namespace diagnostic_indices

namespace optical_indices {
// 3D Optical Properties (cell, level, index)
constexpr int EXTINCTION_COEFF = 0; // 3D extinction coefficient [m⁻¹]
constexpr int SCATTERING_COEFF = 1; // 3D scattering coefficient [m⁻¹]
constexpr int BACKSCATTER_COEFF = 2; // 3D backscatter coefficient [m⁻¹ sr⁻¹]
constexpr int ASYMMETRY_FACTOR = 3; // g asymmetry parameter [fraction -1 to 1]
constexpr int LIDAR_BACKSCATTER = 4; // 3D lidar 180-deg backscatter [m⁻¹ sr⁻¹]

// 2D Column-Integrated Optical Diagnostics (Column AOT at 550nm) (cell,
// level=0, index)
constexpr int EXTINCTION_AOT =
    5; // Total Extinction Aerosol Optical Thickness (AOT)
constexpr int SCATTERING_AOT =
    6; // Total Scattering Aerosol Optical Thickness (AOT)
constexpr int FINE_MODE_EXTINCTION_AOT =
    7; // Fine mode (sub-micron) extinction AOT
constexpr int FINE_MODE_SCATTERING_AOT =
    8;                                  // Fine mode (sub-micron) scattering AOT
constexpr int PM2_5_EXTINCTION_AOT = 9; // PM2.5 extinction AOT
constexpr int PM2_5_SCATTERING_AOT = 10; // PM2.5 scattering AOT

constexpr int NUM_OPTICS = 11;
} // namespace optical_indices

} // namespace exaero
