/// @file ModalRemoval.cpp
/// @brief Transactional host removal utilities for configured modal aerosol.
/// @details No process integration, file I/O, or implicit state mutation.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exaero/ModalCoupling.hpp>
#include <exaero/ModalRemoval.hpp>
#include <exaero/ModalRemovalMath.hpp>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace exaero {
namespace {
void check(bool ok, const char *why) {
  if (!ok)
    throw std::invalid_argument(std::string("FATAL ERROR: ") + why);
}
template <class V>
void checkShape(const V &v, std::initializer_list<std::size_t> extents) {
  check(v.data_handle() != nullptr, "null modal removal array");
  std::size_t d = 0;
  for (auto n : extents)
    check(v.extent(d++) == n, "modal removal array shape mismatch");
}
template <class A, class B> void disjoint(const A &a, const B &b) {
  const auto ab = reinterpret_cast<std::uintptr_t>(a.data_handle());
  const auto bb = reinterpret_cast<std::uintptr_t>(b.data_handle());
  const auto as = a.size() * sizeof(*a.data_handle());
  const auto bs = b.size() * sizeof(*b.data_handle());
  check(ab <= std::numeric_limits<std::uintptr_t>::max() - as &&
            bb <= std::numeric_limits<std::uintptr_t>::max() - bs,
        "modal removal address overflow");
  check(ab + as <= bb || bb + bs <= ab, "modal removal arrays overlap");
}
const ModalConfig &configuration(IAerosolPackage &package) {
  check(package.representation() == AerosolRepresentation::Modal,
        "generic modal package required");
  const auto &config = modalCoupling(package).modalConfig();
  check(!config.modes.empty() && !config.species.empty(),
        "initialize modal package before removal utilities");
  return config;
}
void stateShape(const ModalConfig &cfg, const View3D<const double> &state) {
  check(state.data_handle() && state.extent(0) > 0 && state.extent(1) > 0 &&
            state.extent(2) == cfg.modes.size() * (cfg.species.size() + 1),
        "invalid modal state shape");
}
struct Moment {
  double volume, number, radius;
  bool empty;
};
Moment moment(const ModalConfig &cfg, const View3D<const double> &state,
              std::size_t c, std::size_t k, std::size_t m) {
  const auto ns = cfg.species.size();
  double volume = 0, mass = 0;
  for (std::size_t s = 0; s < ns; ++s) {
    const double value = state(c, k, m * (ns + 1) + s);
    check(std::isfinite(value) && value >= 0, "invalid modal material mass");
    if (!cfg.modes[m].species.empty())
      check(value == 0 ||
                std::find(cfg.modes[m].species.begin(),
                          cfg.modes[m].species.end(),
                          cfg.species[s].name) != cfg.modes[m].species.end(),
            "disallowed material populated");
    mass += value;
    volume += value / cfg.species[s].density;
  }
  const double n = state(c, k, m * (ns + 1) + ns);
  check(std::isfinite(n) && n >= 0 && std::isfinite(mass) &&
            std::isfinite(volume) && (mass > 0) == (n > 0),
        "invalid modal mass and number pair");
  if (n == 0)
    return {0, 0, 0, true};
  return {volume, n, modalRadius(volume, n, cfg.modes[m].geometric_std_dev),
          false};
}
template <class V>
void copyScratch(const std::vector<typename V::element_type> &scratch, V &out) {
  std::copy(scratch.begin(), scratch.end(), out.data_handle());
}
std::size_t index3(std::size_t c, std::size_t k, std::size_t m, std::size_t nc,
                   std::size_t nl) {
  return c + nc * (k + nl * m);
}
/// @brief Check the final signed sink without overflowing double intermediates.
double removalSink(double coefficient, double concentration, double depth) {
  const long double loss =
      static_cast<long double>(coefficient) * concentration / depth;
  check(std::isfinite(loss) && loss >= 0 &&
            loss <= std::numeric_limits<double>::max(),
        "modal removal tendency overflow");
  return loss == 0 ? 0 : -static_cast<double>(loss);
}
/// @brief Translate per-mode moment coefficients; publish both buffers at end.
void removalTendencies(IAerosolPackage &package,
                       const View3D<const double> &state,
                       const View3D<const double> &mass_rate,
                       const View3D<const double> &number_rate,
                       const View2D<const double> *thickness,
                       View3D<double> &mass, View3D<double> &number,
                       bool constituent_rates = false) {
  const auto &cfg = configuration(package);
  stateShape(cfg, state);
  const auto nc = state.extent(0), nl = state.extent(1), nf = state.extent(2),
             nm = cfg.modes.size(), ns = cfg.species.size();
  checkShape(mass_rate, {nc, nl, constituent_rates ? nf : nm});
  checkShape(number_rate, {nc, nl, nm});
  checkShape(mass, {nc, nl, nf});
  checkShape(number, {nc, nl, nm});
  disjoint(state, mass);
  disjoint(state, number);
  disjoint(mass_rate, mass);
  disjoint(mass_rate, number);
  disjoint(number_rate, mass);
  disjoint(number_rate, number);
  disjoint(mass, number);
  if (thickness) {
    checkShape(*thickness, {nc, nl});
    disjoint(*thickness, mass);
    disjoint(*thickness, number);
    for (std::size_t i = 0; i < thickness->size(); ++i)
      check(std::isfinite(thickness->data_handle()[i]) &&
                thickness->data_handle()[i] > 0,
            "modal removal thickness must be finite and positive");
  }
  for (const auto &rates : {mass_rate, number_rate})
    for (std::size_t i = 0; i < rates.size(); ++i)
      check(std::isfinite(rates.data_handle()[i]) &&
                rates.data_handle()[i] >= 0,
            "modal removal coefficient must be finite and nonnegative");
  if (constituent_rates)
    for (std::size_t m=0;m<nm;++m)
      for (std::size_t k=0;k<nl;++k)
        for (std::size_t c=0;c<nc;++c)
          check(mass_rate(c,k,m*(ns+1)+ns)==0,
                "modal mass-rate number slots must be zero");
  std::vector<double> masses(mass.size(), 0), numbers(number.size());
  View3D<double> mass_result(masses.data(), nc, nl, nf),
      number_result(numbers.data(), nc, nl, nm);
  for (std::size_t m = 0; m < nm; ++m) {
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto a = moment(cfg, state, c, k, m);
        number_result(c, k, m) =
            removalSink(number_rate(c, k, m), a.number,
                        thickness ? (*thickness)(c, k) : 1.);
      }
    for (std::size_t s = 0; s < ns; ++s)
      for (std::size_t k = 0; k < nl; ++k)
        for (std::size_t c = 0; c < nc; ++c)
          mass_result(c, k, m * (ns + 1) + s) =
              removalSink(mass_rate(c, k, constituent_rates ? m*(ns+1)+s : m), state(c, k, m * (ns + 1) + s),
                          thickness ? (*thickness)(c, k) : 1.);
  }
  copyScratch(masses, mass);
  copyScratch(numbers, number);
}
} // namespace

void modalRemovalRateTendencies(
    IAerosolPackage &package, const View3D<const double> &state,
    const View3D<const double> &mass_scavenging_rate,
    const View3D<const double> &number_scavenging_rate,
    View3D<double> &mass_tendency, View3D<double> &number_tendency) {
  removalTendencies(package, state, mass_scavenging_rate,
                    number_scavenging_rate, nullptr, mass_tendency,
                    number_tendency);
}

void modalSpeciesRemovalRateTendencies(
    IAerosolPackage &package, const View3D<const double> &state,
    const View3D<const double> &mass_rate,
    const View3D<const double> &number_rate,
    View3D<double> &mass, View3D<double> &number) {
  removalTendencies(package,state,mass_rate,number_rate,nullptr,mass,number,true);
}

void modalRemovalVelocityTendencies(IAerosolPackage &package,
                                    const View2D<const double> &layer_thickness,
                                    const View3D<const double> &state,
                                    const View3D<const double> &velocity,
                                    View3D<double> &mass_tendency,
                                    View3D<double> &number_tendency) {
  removalTendencies(package, state, velocity, velocity, &layer_thickness,
                    mass_tendency, number_tendency);
}

void inspectModalBounds(IAerosolPackage &p, const View3D<const double> &state,
                        View4D<double> &report, View3D<int> &status) {
  const auto &cfg = configuration(p);
  stateShape(cfg, state);
  const auto nc = state.extent(0), nl = state.extent(1), nm = cfg.modes.size();
  checkShape(report, {nc, nl, nm, 3});
  checkShape(status, {nc, nl, nm});
  disjoint(state, report);
  disjoint(state, status);
  disjoint(report, status);
  std::vector<double> r(report.size());
  std::vector<int> st(status.size());
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto a = moment(cfg, state, c, k, m);
        const auto i = index3(c, k, m, nc, nl);
        r[i] = a.radius;
        if (a.empty) {
          st[i] = static_cast<int>(ModalBoundStatus::Empty);
          continue;
        }
        const auto &mode = cfg.modes[m];
        if (!mode.dry_geometric_radius_min) {
          st[i] = static_cast<int>(ModalBoundStatus::Disabled);
          continue;
        }
        const double lo = modalNumberAtRadius(
            a.volume, *mode.dry_geometric_radius_max, mode.geometric_std_dev);
        const double hi = modalNumberAtRadius(
            a.volume, *mode.dry_geometric_radius_min, mode.geometric_std_dev);
        r[i + nc * nl * nm] = lo;
        r[i + 2 * nc * nl * nm] = hi;
        st[i] = static_cast<int>(a.radius < *mode.dry_geometric_radius_min
                                     ? ModalBoundStatus::BelowRadius
                                 : a.radius > *mode.dry_geometric_radius_max
                                     ? ModalBoundStatus::AboveRadius
                                     : ModalBoundStatus::Within);
      }
  copyScratch(r, report);
  copyScratch(st, status);
}

void enforceModalBounds(IAerosolPackage &p, View3D<double> &state,
                        View3D<double> &delta_number, View3D<int> &status) {
  const auto &cfg = configuration(p);
  const View3D<const double> input(state.data_handle(), state.extent(0),
                                   state.extent(1), state.extent(2));
  stateShape(cfg, input);
  const auto nc = state.extent(0), nl = state.extent(1), nm = cfg.modes.size(),
             ns = cfg.species.size();
  checkShape(delta_number, {nc, nl, nm});
  checkShape(status, {nc, nl, nm});
  disjoint(input, delta_number);
  disjoint(input, status);
  disjoint(delta_number, status);
  std::vector<double> report(nc * nl * nm * 3), delta(delta_number.size());
  std::vector<int> st(status.size());
  View4D<double> rv(report.data(), nc, nl, nm, 3);
  View3D<int> sv(st.data(), nc, nl, nm);
  inspectModalBounds(p, input, rv, sv);
  std::vector<double> updated(state.size());
  std::copy(state.data_handle(), state.data_handle() + state.size(),
            updated.begin());
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto i = index3(c, k, m, nc, nl),
                   f = index3(c, k, m * (ns + 1) + ns, nc, nl);
        if (cfg.modes[m].dry_geometric_radius_min &&
            st[i] != static_cast<int>(ModalBoundStatus::Empty)) {
          const double n = std::clamp(updated[f], report[i + nc * nl * nm],
                                      report[i + 2 * nc * nl * nm]);
          delta[i] = n - updated[f];
          updated[f] = n;
          check(std::isfinite(delta[i]) && n > 0,
                "unrepresentable number correction");
        }
      }
  std::copy(updated.begin(), updated.end(), state.data_handle());
  copyScratch(delta, delta_number);
  copyScratch(st, status);
}

void modalScavengingSizes(IAerosolPackage &p, const EnvironmentalStateView &env,
                          const View3D<const double> &state,
                          View4D<double> &diameters, View3D<int> &populated) {
  const auto &cfg = configuration(p);
  stateShape(cfg, state);
  const auto nc = state.extent(0), nl = state.extent(1), nm = cfg.modes.size();
  checkShape(diameters, {nc, nl, nm, 2});
  checkShape(populated, {nc, nl, nm});
  disjoint(state, diameters);
  disjoint(state, populated);
  disjoint(diameters, populated);
  for (const auto &field : {env.temperature, env.relative_humidity}) {
    if (field.data_handle()) {
      disjoint(field, diameters);
      disjoint(field, populated);
    }
  }
  std::vector<double> properties(nc * nl * nm * modal_properties::COUNT),
      out(diameters.size());
  std::vector<int> pop(populated.size());
  View4D<double> pv(properties.data(), nc, nl, nm, modal_properties::COUNT);
  modalCoupling(p).modalProperties(env, state, pv, nullptr);
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto i = index3(c, k, m, nc, nl);
        if (!moment(cfg, state, c, k, m).empty) {
          pop[i] = 1;
          const double dry = 2 * pv(c, k, m, modal_properties::DRY_RADIUS);
          const double wet = 2 * pv(c, k, m, modal_properties::WET_RADIUS);
          check(std::isfinite(dry) && dry > 0 && std::isfinite(wet) &&
                    wet >= dry,
                "invalid modal scavenging diameter");
          out[i] = dry;
          out[i + nc * nl * nm] = wet;
        }
      }
  copyScratch(out, diameters);
  copyScratch(pop, populated);
}

void redistributeModalResuspension(IAerosolPackage &p,
                                   const View3D<const double> &rates,
                                   View3D<double> &out) {
  const auto &cfg = configuration(p);
  stateShape(cfg, rates);
  checkShape(out, {rates.extent(0), rates.extent(1), rates.extent(2)});
  disjoint(rates, out);
  const auto nc = rates.extent(0), nl = rates.extent(1),
             ns = cfg.species.size();
  std::map<std::string, std::size_t> species, modes;
  for (std::size_t s = 0; s < ns; ++s)
    species[cfg.species[s].name] = s;
  for (std::size_t m = 0; m < cfg.modes.size(); ++m)
    modes[cfg.modes[m].name] = m;
  std::vector<double> next(rates.size());
  std::copy(rates.data_handle(), rates.data_handle() + rates.size(),
            next.begin());
  for (std::size_t f = 0; f < rates.extent(2); ++f)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto v = rates(c, k, f);
        check(std::isfinite(v) && v >= 0 && ((f % (ns + 1) != ns) || v == 0),
              "invalid resuspension mass rate or nonzero number slot");
        if (f % (ns + 1) != ns && v > 0) {
          const auto m = f / (ns + 1), s = f % (ns + 1);
          const auto &allowed = cfg.modes[m].species;
          check(allowed.empty() ||
                    std::find(allowed.begin(), allowed.end(),
                              cfg.species[s].name) != allowed.end(),
                "resuspension source contains disallowed material");
        }
      }
  for (const auto &route : cfg.resuspension_routes) {
    const auto s = species.at(route.species), a = modes.at(route.source),
               b = modes.at(route.destination);
    const auto af = a * (ns + 1) + s, bf = b * (ns + 1) + s;
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto v = rates(c, k, af);
        const auto ai = index3(c, k, af, nc, nl), bi = index3(c, k, bf, nc, nl);
        next[ai] -= v;
        next[bi] += v;
        check(std::isfinite(next[bi]), "resuspension destination overflow");
      }
  }
  copyScratch(next, out);
}

void rebinModalDeposition(IAerosolPackage &p,
                          const View3D<const double> &budget,
                          std::size_t species_index,
                          const View1D<const double> &edges,
                          View3D<double> &bins, View3D<double> &tails) {
  const auto &cfg = configuration(p);
  stateShape(cfg, budget);
  check(species_index < cfg.species.size() && edges.data_handle() &&
            edges.extent(0) >= 2,
        "invalid deposited species or diameter edges");
  const auto nc = budget.extent(0), nl = budget.extent(1),
             nb = edges.extent(0) - 1, ns = cfg.species.size();
  checkShape(bins, {nc, nl, nb});
  checkShape(tails, {nc, nl, 2});
  disjoint(budget, bins);
  disjoint(budget, tails);
  disjoint(bins, tails);
  disjoint(edges, bins);
  disjoint(edges, tails);
  for (std::size_t i = 0; i <= nb; ++i) {
    const double x = edges(i);
    check((i == 0 && x == 0) ||
              (i == nb && x == std::numeric_limits<double>::infinity()) ||
              (std::isfinite(x) && x > 0),
          "invalid diameter edge");
    if (i)
      check(edges(i) > edges(i - 1), "diameter edges must increase");
  }
  std::vector<double> b(bins.size()), t(tails.size());
  for (std::size_t m = 0; m < cfg.modes.size(); ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto a = moment(cfg, budget, c, k, m);
        if (a.empty)
          continue;
        const double mass = budget(c, k, m * (ns + 1) + species_index);
        const double first = modalMassFractionBelow(
            edges(0), a.radius, cfg.modes[m].geometric_std_dev);
        t[index3(c, k, 0, nc, nl)] += mass * first;
        t[index3(c, k, 1, nc, nl)] +=
            mass * modalMassFractionAbove(edges(nb), a.radius,
                                          cfg.modes[m].geometric_std_dev);
        for (std::size_t j = 0; j < nb; ++j) {
          const double l = modalMassFractionBelow(
              edges(j), a.radius, cfg.modes[m].geometric_std_dev);
          const double h = modalMassFractionBelow(
              edges(j + 1), a.radius, cfg.modes[m].geometric_std_dev);
          const double fraction =
              l > .5
                  ? modalMassFractionAbove(edges(j), a.radius,
                                           cfg.modes[m].geometric_std_dev) -
                        modalMassFractionAbove(edges(j + 1), a.radius,
                                               cfg.modes[m].geometric_std_dev)
                  : h - l;
          b[index3(c, k, j, nc, nl)] += mass * std::max(0.0, fraction);
        }
      }
  for (double v : b)
    check(std::isfinite(v), "deposition bin overflow");
  for (double v : t)
    check(std::isfinite(v), "deposition tail overflow");
  copyScratch(b, bins);
  copyScratch(t, tails);
}
} // namespace exaero
