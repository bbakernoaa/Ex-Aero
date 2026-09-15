#pragma once
/// @file AerosolIndices.hpp
/// @brief Normative index codes for the diagnostics and optics slabs.
///
/// These constants fix the slot ordering of the (cell, level, index)
/// output arrays written by computeDerivedDiagnostics() and
/// computeOptics(). The ordering is part of the public ABI: host models
/// address output by these codes, so values must never be reordered or
/// reused — only appended (with NUM_* bumped).
namespace exaero {

/// @brief Slot indices for the derived-diagnostics output slab.
///
/// 3D mass/size fields fill (cell, level, index); the column-integrated
/// and surface fields are written at level 0 only.
namespace diagnostic_indices {
constexpr int MASS_CONCENTRATION =
    0; ///< 3D mass concentration per species [kg m^-3]
constexpr int PM2_5_CONCENTRATION =
    1; ///< 3D PM2.5 mass concentration [kg m^-3] (aerodynamic-diameter cut)
constexpr int PM10_CONCENTRATION =
    2; ///< 3D PM10 mass concentration [kg m^-3] (aerodynamic-diameter cut)
constexpr int NUMBER_CONCENTRATION =
    3; ///< Derived number concentration [particles m^-3]
constexpr int SURFACE_AREA_DENSITY =
    4; ///< Surface area density (SAD) [m^2 m^-3]
constexpr int AEROSOL_LIQUID_WATER =
    5; ///< Aerosol liquid water (ALW) [kg m^-3]
constexpr int GRAVITATIONAL_SETTLING_VELOCITY =
    6; ///< Mass-weighted Stokes settling velocity v_g [m s^-1]

constexpr int SURFACE_MASS =
    7; ///< Surface layer mass concentration [kg m^-3] (level 0)
constexpr int COLUMN_MASS =
    8; ///< Column-integrated mass density [kg m^-2] (level 0)
constexpr int SURFACE_PM2_5_MASS =
    9; ///< Surface layer PM2.5 mass concentration [kg m^-3] (level 0)
constexpr int COLUMN_PM2_5_MASS =
    10; ///< Column-integrated PM2.5 mass density [kg m^-2] (level 0)

constexpr int NUM_DIAGNOSTICS = 11; ///< Total diagnostic slots (slab extent).
} // namespace diagnostic_indices

/// @brief Slot indices for the optical-properties output slab.
///
/// 3D volume coefficients fill (cell, level, index); the column AOT
/// values are integrated over @f$ \Delta z @f$ and written at level 0.
namespace optical_indices {
constexpr int EXTINCTION_COEFF = 0; ///< 3D extinction coefficient [m^-1]
constexpr int SCATTERING_COEFF = 1; ///< 3D scattering coefficient [m^-1]
constexpr int BACKSCATTER_COEFF =
    2; ///< 3D hemispheric backscatter coefficient [m^-1 sr^-1]
constexpr int ASYMMETRY_FACTOR =
    3; ///< g asymmetry parameter [fraction, -1 to 1]
constexpr int LIDAR_BACKSCATTER =
    4; ///< 3D lidar 180-degree backscatter coefficient [m^-1 sr^-1]

constexpr int EXTINCTION_AOT =
    5; ///< Total extinction aerosol optical thickness (AOT), column-integrated
constexpr int SCATTERING_AOT =
    6; ///< Total scattering aerosol optical thickness (AOT), column-integrated
constexpr int FINE_MODE_EXTINCTION_AOT =
    7; ///< Fine mode (sub-micron) extinction AOT
constexpr int FINE_MODE_SCATTERING_AOT =
    8; ///< Fine mode (sub-micron) scattering AOT
constexpr int PM2_5_EXTINCTION_AOT = 9;  ///< PM2.5-cut extinction AOT
constexpr int PM2_5_SCATTERING_AOT = 10; ///< PM2.5-cut scattering AOT

constexpr int NUM_OPTICS = 11; ///< Total optics slots (slab extent).
} // namespace optical_indices

} // namespace exaero
