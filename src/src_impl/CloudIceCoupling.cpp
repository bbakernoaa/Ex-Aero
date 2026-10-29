/// @file CloudIceCoupling.cpp
/// @brief Common cloud/ice dispatch; host-owned reservoirs, unchanged kernels.
#include "HostCouplingInternal.hpp"
#include <exaero/CloudIceCoupling.hpp>
namespace exaero {
void cloudIce(IAerosolPackage &p, const EnvironmentalStateView &env,
              const View3D<const double> &ambient,
              const View3D<const double> *cloudborne, IceCategory category,
              bool ice, View4D<double> &a, View4D<double> &b, View3D<int> &as,
              View3D<int> &bs) {
  host_detail::require(supportsHostOperation(p, HostOperation::CloudIce),
                       "package lacks host cloud/ice support");
  const auto info = hostCouplingInfo(p);
  const auto nc = ambient.extent(0), nl = ambient.extent(1),
             np = info.populations.size();
  host_detail::shape(a, {nc, nl, np, cloud_ice_properties::COUNT});
  host_detail::shape(b, {nc, nl, np, cloud_ice_properties::COUNT});
  host_detail::shape(as, {nc, nl, np});
  host_detail::shape(bs, {nc, nl, np});
  host_detail::environmentInputs(env, nc, nl, a, b, as, bs);
  if (p.representation() == AerosolRepresentation::Bulk)
    bulkCloudIce(p, env, ambient, cloudborne, category, ice, a, b, as, bs);
  else
    modalCloudIce(p, env, ambient, cloudborne, category, ice, a, b, as, bs);
}
void activationFractions(IAerosolPackage &p, const EnvironmentalStateView &env,
                         const View3D<const double> &state,
                         const View1D<const double> &ss, View4D<double> &out) {
  host_detail::require(supportsHostOperation(p, HostOperation::Activation),
                       "package lacks host activation support");
  const auto info = hostCouplingInfo(p);
  const auto nc = state.extent(0), nl = state.extent(1),
             np = info.populations.size();
  host_detail::require(ss.extent(0) > 0 &&
                           ss.extent(0) <=
                               std::numeric_limits<std::size_t>::max() / 2,
                       "invalid activation supersaturation count");
  host_detail::shape(out, {nc, nl, np, 2 * ss.extent(0)});
  host_detail::environmentInputs(env, nc, nl, out);
  if (p.representation() == AerosolRepresentation::Bulk)
    bulkActivationFractions(p, env, state, ss, out);
  else
    modalActivationFractions(p, env, state, ss, out);
}
} // namespace exaero
