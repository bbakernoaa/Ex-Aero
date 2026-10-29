/// @file DepositionCoupling.cpp
/// @brief Snapshot-based material routing and dispatch to existing deposition
/// kernels.
#include "HostCouplingInternal.hpp"
#include <exaero/BulkRemoval.hpp>
#include <exaero/DepositionCoupling.hpp>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemoval.hpp>
namespace exaero {
std::vector<ResuspensionFieldRoute>
configuredResuspensionRoutes(IAerosolPackage &p) {
  const auto info = hostCouplingInfo(p);
  std::vector<ResuspensionFieldRoute> result;
  if (info.representation == AerosolRepresentation::Bulk)
    return result;
  const auto &cfg = modalCoupling(p).modalConfig();
  for (const auto &route : cfg.resuspension_routes) {
    auto field = [&](const std::string &mode) {
      const auto pop =
          std::find_if(info.populations.begin(), info.populations.end(),
                       [&](const auto &v) { return v.name == mode; });
      host_detail::require(pop != info.populations.end(),
                           "invalid configured resuspension population");
      const auto material =
          std::find_if(info.constituents.begin(), info.constituents.end(),
                       [&](const auto &v) {
                         return v.population_index == pop->index &&
                                v.name == route.species && v.permitted;
                       });
      host_detail::require(material != info.constituents.end(),
                           "invalid configured resuspension material");
      return std::size_t(material->mass_field);
    };
    result.push_back({field(route.source), field(route.destination)});
  }
  return result;
}
void redistributeResuspension(IAerosolPackage &p,
                              const View3D<const double> &rates,
                              std::span<const ResuspensionFieldRoute> routes,
                              View3D<double> &out) {
  using namespace host_detail;
  const auto info = hostCouplingInfo(p);
  const auto nc = rates.extent(0), nl = rates.extent(1),
             nf = std::size_t(info.field_count);
  require(nc > 0 && nl > 0 &&
              nc <= std::numeric_limits<std::size_t>::max() / nl / nf /
                        sizeof(double),
          "resuspension extent overflow");
  shape(rates, {nc, nl, nf});
  shape(out, {nc, nl, nf});
  disjoint(rates, out);
  if (!routes.empty())
    separate(routes.data(), routes.size_bytes(), out.data_handle(),
             out.size() * sizeof(double));
  std::vector<const ConstituentInfo *> fields(nf, nullptr);
  for (const auto &v : info.constituents)
    fields[v.mass_field] = &v;
  for (std::size_t f = 0; f < nf; ++f)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto value = rates(c, k, f);
        require(std::isfinite(value) && value >= 0 &&
                    (value == 0 || (fields[f] && fields[f]->permitted)),
                "invalid resuspension mass rate or number slot");
      }
  std::vector<int> destination(nf, -1);
  for (const auto &r : routes) {
    require(r.source_field < nf && r.destination_field < nf &&
                r.source_field != r.destination_field,
            "invalid resuspension field route");
    const auto *a = fields[r.source_field], *b = fields[r.destination_field];
    require(a && b && a->permitted && b->permitted &&
                destination[r.source_field] < 0,
            "duplicate or prohibited resuspension route");
    require(info.representation == AerosolRepresentation::Bulk ||
                a->species_index == b->species_index,
            "resuspension cannot change modal material identity");
    destination[r.source_field] = int(r.destination_field);
  }
  for (std::size_t f = 0; f < nf; ++f) {
    std::vector<bool> visited(nf, false);
    for (int next = int(f); next >= 0; next = destination[next]) {
      require(!visited[next], "resuspension routes contain a cycle");
      visited[next] = true;
    }
  }
  std::vector<double> scratch(rates.data_handle(),
                              rates.data_handle() + rates.size());
  View3D<double> result(scratch.data(), nc, nl, nf);
  // Remove all original routed sources before accumulating incoming material.
  // This avoids artificial intermediate overflow for a route chain.
  for (const auto &r : routes)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c)
        result(c, k, r.source_field) = 0;
  for (const auto &r : routes)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const long double value =
            static_cast<long double>(result(c, k, r.destination_field)) +
            rates(c, k, r.source_field);
        require(std::isfinite(value) &&
                    value <= std::numeric_limits<double>::max(),
                "resuspension destination overflow");
        result(c, k, r.destination_field) = static_cast<double>(value);
      }
  std::copy(scratch.begin(), scratch.end(), out.data_handle());
}
void rebinDeposition(IAerosolPackage &p, const View3D<const double> &budget,
                     std::size_t species, const View1D<const double> &edges,
                     View3D<double> &bins, View3D<double> &tails) {
  const auto info = hostCouplingInfo(p);
  if (info.representation == AerosolRepresentation::Bulk)
    rebinBulkDeposition(p, budget, species, edges, bins, tails);
  else
    rebinModalDeposition(p, budget, species, edges, bins, tails);
}
} // namespace exaero
