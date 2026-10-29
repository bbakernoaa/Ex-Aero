/// @file HostCoupling.cpp
/// @brief Representation-aware population metadata, without changing package
/// vtables.
#include "HostCouplingInternal.hpp"
#include <exaero/BulkRemoval.hpp>
#include <exaero/CloudIceCoupling.hpp>
#include <exaero/ModalCoupling.hpp>
namespace exaero {
HostCouplingInfo hostCouplingInfo(IAerosolPackage &p) {
  using host_detail::require;
  const auto fields = p.stateFields();
  require(!fields.empty() &&
              fields.size() <= std::size_t(std::numeric_limits<int>::max()),
          "initialize a supported package before querying host metadata");
  HostCouplingInfo info{p.representation(), int(fields.size()), {}, {}};
  if (info.representation == AerosolRepresentation::Bulk) {
    require(dynamic_cast<IBulkRemovalMetadata *>(&p),
            "package lacks bulk physical metadata");
    for (int f = 0; f < info.field_count; ++f) {
      require(fields[f].quantity == StateQuantity::Mass &&
                  fields[f].mode == -1 && !fields[f].cloudborne,
              "invalid bulk field metadata");
      info.populations.push_back(
          {f, fields[f].name, NumberKind::Diagnosed, -1});
      info.constituents.push_back({f, f, f, fields[f].name, true});
    }
  } else {
    require(info.representation == AerosolRepresentation::Modal,
            "unknown aerosol representation");
    const auto &cfg = modalCoupling(p).modalConfig();
    const auto ns = cfg.species.size(), nm = cfg.modes.size();
    require(ns > 0 && nm > 0 && ns < fields.size() &&
                nm == fields.size() / (ns + 1) &&
                nm * (ns + 1) == fields.size(),
            "unsupported modal state layout");
    for (std::size_t m = 0; m < nm; ++m) {
      const int number = int(m * (ns + 1) + ns);
      require(fields[number].quantity == StateQuantity::Number &&
                  fields[number].mode == int(m) && !fields[number].cloudborne,
              "invalid modal number metadata");
      info.populations.push_back(
          {int(m), cfg.modes[m].name, NumberKind::Prognostic, number});
      for (std::size_t s = 0; s < ns; ++s) {
        const int f = int(m * (ns + 1) + s);
        require(fields[f].quantity == StateQuantity::Mass &&
                    fields[f].mode == int(m) && !fields[f].cloudborne,
                "invalid modal constituent metadata");
        const auto &allowed = cfg.modes[m].species;
        const bool permitted =
            allowed.empty() || std::find(allowed.begin(), allowed.end(),
                                         cfg.species[s].name) != allowed.end();
        info.constituents.push_back(
            {int(m), int(s), f, cfg.species[s].name, permitted});
      }
    }
  }
  return info;
}
bool supportsHostOperation(IAerosolPackage &p, HostOperation op) {
  host_detail::require(int(op) >= 0 && int(op) <= 8,
                       "invalid host operation code");
  try {
    const auto info = hostCouplingInfo(p);
    if (op == HostOperation::CloudIce)
      return info.representation == AerosolRepresentation::Modal ||
             dynamic_cast<IBulkCloudIceMetadata *>(&p) != nullptr;
    if (op == HostOperation::Resuspension ||
        op == HostOperation::DepositionRebin)
      return true;
    if (op == HostOperation::PopulationProperties)
      return true;
    if (op == HostOperation::Activation)
      return true;
    if (op == HostOperation::RemovalRates ||
        op == HostOperation::RemovalVelocity)
      return true;
    return (op == HostOperation::InspectNumberBounds ||
            op == HostOperation::EnforceNumberBounds) &&
           info.representation == AerosolRepresentation::Modal;
  } catch (const std::logic_error &) {
    return false;
  }
}
} // namespace exaero
