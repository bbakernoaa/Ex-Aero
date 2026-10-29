/// @file CloudIceCoupling.cpp
/// @brief Common cloud/ice dispatch; host-owned reservoirs, unchanged kernels.
#include "HostCouplingInternal.hpp"
#include <exaero/CloudIceCoupling.hpp>
namespace exaero {
void cloudIce(IAerosolPackage &p,const EnvironmentalStateView &env,
    const View3D<const double> &ambient,const View3D<const double> *cloudborne,
    IceCategory category,bool ice,View4D<double> &a,View4D<double> &b,
    View3D<int> &as,View3D<int> &bs) {
  host_detail::require(supportsHostOperation(p,HostOperation::CloudIce),"package lacks host cloud/ice support");
  if(p.representation()==AerosolRepresentation::Bulk)bulkCloudIce(p,env,ambient,cloudborne,category,ice,a,b,as,bs);
  else modalCloudIce(p,env,ambient,cloudborne,category,ice,a,b,as,bs);
}
void activationFractions(IAerosolPackage &p,const EnvironmentalStateView &env,
    const View3D<const double> &state,const View1D<const double> &ss,View4D<double> &out) {
  host_detail::require(supportsHostOperation(p,HostOperation::Activation),"package lacks host activation support");
  if(p.representation()==AerosolRepresentation::Bulk)bulkActivationFractions(p,env,state,ss,out);
  else modalActivationFractions(p,env,state,ss,out);
}
}
