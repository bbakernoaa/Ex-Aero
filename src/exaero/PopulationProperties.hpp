#pragma once
/// @file PopulationProperties.hpp
/// @brief Physical population properties, independent of optical tables.
#include <exaero/HostCoupling.hpp>
namespace exaero {
enum class PopulationProperty : int {
  Mass=0, Number, DryNumberMedianRadius, WetNumberMedianRadius,
  DryScavengingDiameter, WetScavengingDiameter, DryVolume, DrySurfaceArea,
  Kappa, DryDensity, WetDensity, GeometricStdDev, LowerDryRadius,
  UpperDryRadius, EffectiveDryRadius, Count
};
enum class PropertyAvailability : int { Unavailable=0, Configured=1, Derived=2 };
enum class PopulationStatus : int { Empty=0, Populated=1 };
/// @brief Populate (cell,level,P,Count) values/availability and (cell,level,P) status.
/// @details Mass kg/m3, number particles/m3, radii/diameters m, volume m3/m3,
/// dry area m2/m3, densities kg/m3 of particle volume; kappa/width dimensionless.
/// Empty populations retain configured dry reference geometry; empty modal
/// composition and wet properties are unavailable. Numerical unavailable slots
/// are zero and must be interpreted with availability. No state is mutated;
/// disjoint outputs publish atomically, including on unsupported-package errors.
void populationProperties(IAerosolPackage&,const EnvironmentalStateView&,
    const View3D<const double>& state,View4D<double>& values,
    View4D<int>& availability,View3D<int>& population_status);
}
