/// @file BulkRemoval.cpp
/// @brief Transactional host utilities for mass-only bulk removal coupling.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exaero/BulkRemoval.hpp>
#include <exaero/ModalRemovalMath.hpp>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace exaero {
namespace {
void require(bool condition, const char *message) {
  if (!condition)
    throw std::invalid_argument(std::string("FATAL ERROR: ") + message);
}
template <class V>
void shape(const V &v, std::initializer_list<std::size_t> sizes) {
  require(v.data_handle() != nullptr, "null bulk removal array");
  std::size_t dim = 0;
  for (const auto size : sizes)
    require(v.extent(dim++) == size, "bulk removal array shape mismatch");
}
template <class A, class B> void disjoint(const A &a, const B &b) {
  const auto an = a.size() * sizeof(*a.data_handle());
  const auto bn = b.size() * sizeof(*b.data_handle());
  const auto ap = reinterpret_cast<std::uintptr_t>(a.data_handle());
  const auto bp = reinterpret_cast<std::uintptr_t>(b.data_handle());
  require(ap <= std::numeric_limits<std::uintptr_t>::max() - an &&
              bp <= std::numeric_limits<std::uintptr_t>::max() - bn &&
              (ap + an <= bp || bp + bn <= ap),
          "bulk removal arrays overlap or address overflow");
}
IBulkRemovalMetadata &metadata(IAerosolPackage &package) {
  require(package.representation() == AerosolRepresentation::Bulk,
          "bulk package required");
  auto *bulk = dynamic_cast<IBulkRemovalMetadata *>(&package);
  require(bulk != nullptr, "package has no bulk removal metadata");
  return *bulk;
}
std::size_t speciesCount(IAerosolPackage &package) {
  const auto fields = package.stateFields();
  require(!fields.empty(), "initialize bulk package before removal utilities");
  for (const auto &field : fields)
    require(field.quantity == StateQuantity::Mass && field.mode < 0,
            "bulk removal requires mass-only species fields");
  return fields.size();
}
void stateShape(IAerosolPackage &package, const View3D<const double> &state) {
  require(state.data_handle() && state.extent(0) && state.extent(1) &&
              state.extent(2) == speciesCount(package),
          "invalid bulk state shape");
}
std::size_t index3(std::size_t c, std::size_t k, std::size_t s, std::size_t nc,
                   std::size_t nl) {
  return c + nc * (k + nl * s);
}
void validateMass(const View3D<const double> &input) {
  for (std::size_t i = 0; i < input.size(); ++i)
    require(std::isfinite(input.data_handle()[i]) &&
                input.data_handle()[i] >= 0,
            "bulk material must be finite and nonnegative");
}
/// @brief Fixed lognormal mean particle mass; independent of optical tables.
double particleMass(const BulkSpeciesSize &size) {
  require(std::isfinite(size.dry_density) && size.dry_density > 0 &&
              std::isfinite(size.lognormal_dg) && size.lognormal_dg > 0 &&
              std::isfinite(size.lognormal_sigma) && size.lognormal_sigma >= 1,
          "invalid fixed bulk number size parameters");
  const double l = std::log(size.lognormal_sigma);
  const double mass = (std::acos(-1.0) / 6.0) * size.dry_density *
                      size.lognormal_dg * size.lognormal_dg *
                      size.lognormal_dg * std::exp(4.5 * l * l);
  require(std::isfinite(mass) && mass > 0,
          "unrepresentable fixed bulk particle mass");
  return mass;
}
/// @brief Validate a staged loss before conversion to a signed double sink.
double negativeSink(long double loss) {
  require(std::isfinite(loss) && loss >= 0 &&
              loss <= std::numeric_limits<double>::max(),
          "bulk removal tendency overflow");
  return loss == 0 ? 0 : -static_cast<double>(loss);
}
/// @brief Common synchronous rate/velocity conversion with atomic publication.
void removalTendencies(IAerosolPackage &package,
                       const View3D<const double> &state,
                       const View3D<const double> &rate,
                       const View2D<const double> *thickness,
                       View3D<double> &mass, View3D<double> &number) {
  auto &meta = metadata(package);
  stateShape(package, state);
  const auto nc = state.extent(0), nl = state.extent(1), ns = state.extent(2);
  shape(rate, {nc, nl, ns});
  shape(mass, {nc, nl, ns});
  shape(number, {nc, nl, ns});
  disjoint(state, mass);
  disjoint(state, number);
  disjoint(rate, mass);
  disjoint(rate, number);
  disjoint(mass, number);
  if (thickness) {
    shape(*thickness, {nc, nl});
    disjoint(*thickness, mass);
    disjoint(*thickness, number);
    for (std::size_t i = 0; i < thickness->size(); ++i)
      require(std::isfinite(thickness->data_handle()[i]) &&
                  thickness->data_handle()[i] > 0,
              "bulk removal thickness must be finite and positive");
  }
  validateMass(state);
  validateMass(rate);
  std::vector<double> masses(mass.size()), numbers(number.size());
  View3D<double> mass_result(masses.data(), nc, nl, ns),
      number_result(numbers.data(), nc, nl, ns);
  for (std::size_t s = 0; s < ns; ++s) {
    const double particle_mass = particleMass(meta.bulkSpeciesSize(s));
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const long double loss = static_cast<long double>(rate(c, k, s)) *
                                 state(c, k, s) /
                                 (thickness ? (*thickness)(c, k) : 1.);
        mass_result(c, k, s) = negativeSink(loss);
        number_result(c, k, s) = negativeSink(loss / particle_mass);
      }
  }
  std::copy(masses.begin(), masses.end(), mass.data_handle());
  std::copy(numbers.begin(), numbers.end(), number.data_handle());
}
void edgesValid(const View1D<const double> &edges) {
  require(edges.data_handle() && edges.extent(0) >= 2,
          "invalid bulk diameter edges");
  const auto n = edges.extent(0);
  for (std::size_t i = 0; i < n; ++i) {
    const double v = edges(i);
    require((i == 0 && v == 0) ||
                (i == n - 1 && v == std::numeric_limits<double>::infinity()) ||
                (std::isfinite(v) && v > 0),
            "invalid bulk diameter edge");
    if (i)
      require(edges(i) > edges(i - 1), "bulk diameter edges must increase");
  }
}
} // namespace

void bulkFixedNumber(IAerosolPackage &package,
                     const View3D<const double> &state,
                     View3D<double> &number) {
  auto &meta = metadata(package);
  stateShape(package, state);
  const auto nc = state.extent(0), nl = state.extent(1), ns = state.extent(2);
  shape(number, {nc, nl, ns});
  disjoint(state, number);
  validateMass(state);
  std::vector<double> out(number.size());
  for (std::size_t s = 0; s < ns; ++s) {
    const double particle_mass = particleMass(meta.bulkSpeciesSize(s));
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c)
        out[index3(c, k, s, nc, nl)] = state(c, k, s) / particle_mass;
  }
  for (double value : out)
    require(std::isfinite(value), "bulk number concentration overflow");
  std::copy(out.begin(), out.end(), number.data_handle());
}

void bulkRemovalRateTendencies(IAerosolPackage &package,
                               const View3D<const double> &state,
                               const View3D<const double> &mass_scavenging_rate,
                               View3D<double> &mass_tendency,
                               View3D<double> &number_tendency) {
  removalTendencies(package, state, mass_scavenging_rate, nullptr,
                    mass_tendency, number_tendency);
}

void bulkRemovalVelocityTendencies(IAerosolPackage &package,
                                   const View2D<const double> &layer_thickness,
                                   const View3D<const double> &state,
                                   const View3D<const double> &velocity,
                                   View3D<double> &mass_tendency,
                                   View3D<double> &number_tendency) {
  removalTendencies(package, state, velocity, &layer_thickness, mass_tendency,
                    number_tendency);
}

void bulkReferenceRadii(IAerosolPackage &package, View2D<double> &radii,
                        View2D<int> &status) {
  auto &meta = metadata(package);
  const auto ns = speciesCount(package);
  shape(radii, {ns, 3});
  shape(status, {ns, 3});
  disjoint(radii, status);
  std::vector<double> values(radii.size());
  std::vector<int> codes(status.size());
  for (std::size_t s = 0; s < ns; ++s) {
    const auto size = meta.bulkSpeciesSize(s);
    require(std::isfinite(size.lognormal_dg) && size.lognormal_dg > 0 &&
                std::isfinite(size.lognormal_sigma) &&
                size.lognormal_sigma >= 1,
            "invalid configured bulk size distribution");
    require((size.dry_radius_lower == 0 && size.dry_radius_upper == 0) ||
                (std::isfinite(size.dry_radius_lower) &&
                 std::isfinite(size.dry_radius_upper) &&
                 size.dry_radius_lower > 0 &&
                 size.dry_radius_lower < size.dry_radius_upper),
            "invalid configured bulk radius bounds");
    require(size.dry_effective_radius == 0 ||
                (std::isfinite(size.dry_effective_radius) &&
                 size.dry_effective_radius > 0),
            "invalid configured bulk effective radius");
    if (size.dry_radius_lower > 0) {
      values[s] = size.dry_radius_lower;
      values[s + ns] = size.dry_radius_upper;
      codes[s] = codes[s + ns] = static_cast<int>(BulkRadiusStatus::Configured);
    }
    const double l = std::log(size.lognormal_sigma);
    const double derived = size.lognormal_dg / 2 * std::exp(2.5 * l * l);
    const double effective =
        size.dry_effective_radius > 0 ? size.dry_effective_radius : derived;
    require(std::isfinite(effective) && effective > 0,
            "unrepresentable bulk effective radius");
    values[s + 2 * ns] = effective;
    codes[s + 2 * ns] = static_cast<int>(size.dry_effective_radius > 0
                                             ? BulkRadiusStatus::Configured
                                             : BulkRadiusStatus::Derived);
  }
  std::copy(values.begin(), values.end(), radii.data_handle());
  std::copy(codes.begin(), codes.end(), status.data_handle());
}

void bulkScavengingSizes(IAerosolPackage &package,
                         const View2D<const double> &rh,
                         const View3D<const double> &state,
                         View4D<double> &diameters, View3D<int> &size_status) {
  auto &meta = metadata(package);
  stateShape(package, state);
  const auto nc = state.extent(0), nl = state.extent(1), ns = state.extent(2);
  shape(rh, {nc, nl});
  shape(diameters, {nc, nl, ns, 2});
  shape(size_status, {nc, nl, ns});
  disjoint(state, diameters);
  disjoint(state, size_status);
  disjoint(rh, diameters);
  disjoint(rh, size_status);
  disjoint(diameters, size_status);
  validateMass(state);
  for (std::size_t i = 0; i < rh.size(); ++i)
    require(std::isfinite(rh.data_handle()[i]) && rh.data_handle()[i] >= 0 &&
                rh.data_handle()[i] <= 1,
            "invalid bulk relative humidity");
  std::vector<double> out(diameters.size());
  std::vector<int> statuses(size_status.size(), BULK_SIZE_EMPTY);
  const auto points = nc * nl;
  for (std::size_t s = 0; s < ns; ++s) {
    const auto size = meta.bulkSpeciesSize(s);
    require(std::isfinite(size.dry_particle_diameter) &&
                size.dry_particle_diameter > 0 &&
                std::isfinite(size.dry_density) && size.dry_density > 0 &&
                std::isfinite(size.hygroscopicity) && size.hygroscopicity >= 0,
            "invalid configured bulk dry size, density or hygroscopicity");
    bool populated = false;
    for (std::size_t i = 0; i < points; ++i)
      populated |= state.data_handle()[i + points * s] > 0;
    if (!populated)
      continue;
    for (std::size_t i = 0; i < points; ++i) {
      if (state.data_handle()[i + points * s] == 0)
        continue;
      const auto capped_rh = std::min(rh.data_handle()[i], .99);
      const auto growth =
          std::cbrt(1 + size.hygroscopicity * capped_rh / (1 - capped_rh));
      const auto wet = size.dry_particle_diameter * growth;
      require(std::isfinite(wet) && wet >= size.dry_particle_diameter,
              "invalid configured bulk wet size");
      statuses[i + points * s] = BULK_SIZE_CONFIGURED;
      out[i + points * s] = size.dry_particle_diameter;
      out[i + points * (s + ns)] = wet;
    }
  }
  for (const double value : out)
    require(std::isfinite(value), "bulk diameter overflow");
  std::copy(out.begin(), out.end(), diameters.data_handle());
  std::copy(statuses.begin(), statuses.end(), size_status.data_handle());
}

void redistributeBulkResuspension(IAerosolPackage &package,
                                  const View3D<const double> &rates,
                                  std::span<const BulkResuspensionRoute> routes,
                                  View3D<double> &routed_rates) {
  (void)metadata(package);
  stateShape(package, rates);
  shape(routed_rates, {rates.extent(0), rates.extent(1), rates.extent(2)});
  disjoint(rates, routed_rates);
  validateMass(rates);
  const auto nc = rates.extent(0), nl = rates.extent(1), ns = rates.extent(2);
  std::vector<int> destination(ns, -1);
  for (const auto &route : routes) {
    require(route.source < ns && route.destination < ns &&
                route.source != route.destination &&
                destination[route.source] < 0,
            "invalid or duplicate bulk resuspension route");
    destination[route.source] = static_cast<int>(route.destination);
  }
  for (std::size_t s = 0; s < ns; ++s) {
    std::vector<bool> visited(ns);
    for (int next = static_cast<int>(s); next >= 0; next = destination[next]) {
      require(!visited[next], "bulk resuspension routes contain a cycle");
      visited[next] = true;
    }
  }
  std::vector<double> out(rates.data_handle(),
                          rates.data_handle() + rates.size());
  for (std::size_t s = 0; s < ns; ++s) {
    if (destination[s] < 0)
      continue;
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto src = index3(c, k, s, nc, nl);
        const auto dst = index3(c, k, destination[s], nc, nl);
        out[src] -= rates.data_handle()[src];
        out[dst] += rates.data_handle()[src];
        require(std::isfinite(out[dst]),
                "bulk resuspension destination overflow");
      }
  }
  std::copy(out.begin(), out.end(), routed_rates.data_handle());
}

void rebinBulkDeposition(IAerosolPackage &package,
                         const View3D<const double> &budget,
                         std::size_t species_index,
                         const View1D<const double> &edges,
                         View3D<double> &bins, View3D<double> &tails) {
  auto &meta = metadata(package);
  stateShape(package, budget);
  require(species_index < budget.extent(2), "invalid bulk deposition species");
  edgesValid(edges);
  const auto nc = budget.extent(0), nl = budget.extent(1),
             nb = edges.extent(0) - 1;
  shape(bins, {nc, nl, nb});
  shape(tails, {nc, nl, 2});
  disjoint(budget, bins);
  disjoint(budget, tails);
  disjoint(edges, bins);
  disjoint(edges, tails);
  disjoint(bins, tails);
  validateMass(budget);
  const auto size = meta.bulkSpeciesSize(species_index);
  require(std::isfinite(size.lognormal_dg) && size.lognormal_dg > 0 &&
              std::isfinite(size.lognormal_sigma) && size.lognormal_sigma >= 1,
          "invalid bulk deposited size distribution");
  std::vector<double> out(bins.size()), remainder(tails.size());
  for (std::size_t k = 0; k < nl; ++k)
    for (std::size_t c = 0; c < nc; ++c) {
      const auto mass = budget(c, k, species_index);
      if (mass == 0)
        continue;
      if (size.lognormal_sigma == 1) {
        if (size.lognormal_dg < edges(0))
          remainder[index3(c, k, 0, nc, nl)] = mass;
        else if (size.lognormal_dg >= edges(nb))
          remainder[index3(c, k, 1, nc, nl)] = mass;
        else {
          const auto j = std::upper_bound(edges.data_handle(),
                                          edges.data_handle() + nb + 1,
                                          size.lognormal_dg) -
                         edges.data_handle() - 1;
          out[index3(c, k, j, nc, nl)] = mass;
        }
      } else {
        const auto radius = size.lognormal_dg / 2;
        const auto sigma = size.lognormal_sigma;
        remainder[index3(c, k, 0, nc, nl)] =
            mass * modalMassFractionBelow(edges(0), radius, sigma);
        remainder[index3(c, k, 1, nc, nl)] =
            mass * modalMassFractionAbove(edges(nb), radius, sigma);
        for (std::size_t j = 0; j < nb; ++j) {
          const double low = modalMassFractionBelow(edges(j), radius, sigma);
          const double high =
              modalMassFractionBelow(edges(j + 1), radius, sigma);
          const double fraction =
              low > .5 ? modalMassFractionAbove(edges(j), radius, sigma) -
                             modalMassFractionAbove(edges(j + 1), radius, sigma)
                       : high - low;
          out[index3(c, k, j, nc, nl)] = mass * std::max(0.0, fraction);
        }
      }
    }
  for (double value : out)
    require(std::isfinite(value), "bulk deposition bin overflow");
  for (double value : remainder)
    require(std::isfinite(value), "bulk deposition tail overflow");
  std::copy(out.begin(), out.end(), bins.data_handle());
  std::copy(remainder.begin(), remainder.end(), tails.data_handle());
}
} // namespace exaero
