#pragma once
/// @file HostCoupling.hpp
/// @brief Initialized population metadata and optional host-helper capabilities.
/// @details Value snapshots; re-query after reinitialization. No state update or I/O.
#include <exaero/IAerosolPackage.hpp>
namespace exaero {
enum class NumberKind : int { Diagnosed = 0, Prognostic = 1 };
enum class HostOperation : int {
  PopulationProperties = 0, CloudIce = 1, Activation = 2, RemovalRates = 3,
  RemovalVelocity = 4, Resuspension = 5, DepositionRebin = 6,
  InspectNumberBounds = 7, EnforceNumberBounds = 8
};
struct PopulationInfo {
  int index;
  std::string name;
  NumberKind number_kind;
  int number_field; ///< -1 for diagnosed bulk number.
};
struct ConstituentInfo {
  int population_index, species_index, mass_field;
  std::string name;
  bool permitted; ///< False for allocated modal material slots prohibited by configuration.
};
struct HostCouplingInfo {
  AerosolRepresentation representation;
  int field_count;
  std::vector<PopulationInfo> populations;
  std::vector<ConstituentInfo> constituents;
};
/// @brief Query zero-based fields without parsing names; throws if unsupported/uninitialized.
HostCouplingInfo hostCouplingInfo(IAerosolPackage &package);
/// @brief Known unsupported operations return false; an invalid operation code throws.
bool supportsHostOperation(IAerosolPackage &package, HostOperation operation);
} // namespace exaero
