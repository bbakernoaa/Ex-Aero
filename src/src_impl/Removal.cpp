/// @file Removal.cpp
/// @brief Representation dispatch for host-owned particle and constituent removal.
#include "HostCouplingInternal.hpp"
#include <exaero/BulkRemoval.hpp>
#include <exaero/ModalRemoval.hpp>
#include <exaero/Removal.hpp>
namespace exaero {
void removalRateTendencies(IAerosolPackage &p,const View3D<const double> &state,
    const View3D<const double> &mass_rate,const View3D<const double> *number_rate,
    View3D<double> &mass,View3D<double> &number) {
  const auto info=hostCouplingInfo(p);
  if(info.representation==AerosolRepresentation::Bulk) {
    host_detail::require(!number_rate,"bulk number loss follows mass loss; omit number-rate input");
    bulkRemovalRateTendencies(p,state,mass_rate,mass,number);
  } else {
    host_detail::require(number_rate,"modal removal requires a mode number-rate input");
    modalSpeciesRemovalRateTendencies(p,state,mass_rate,*number_rate,mass,number);
  }
}
void removalVelocityTendencies(IAerosolPackage &p,const View2D<const double> &dz,
    const View3D<const double> &state,const View3D<const double> &velocity,
    View3D<double> &mass,View3D<double> &number) {
  const auto info=hostCouplingInfo(p);
  if(info.representation==AerosolRepresentation::Bulk)
    bulkRemovalVelocityTendencies(p,dz,state,velocity,mass,number);
  else modalRemovalVelocityTendencies(p,dz,state,velocity,mass,number);
}
}
