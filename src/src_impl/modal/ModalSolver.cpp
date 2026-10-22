/// @file ModalSolver.cpp
/// @brief Generic fixed-width modal Kokkos kernels, isolated from mdspan.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: ModalPackage pointer bridge; no file I/O.
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <exaero/AerosolIndices.hpp>
#include <modal/ModalSolver.hpp>
#include <stdexcept>
#ifdef EXAERO_ENABLE_MAM4XX
#include <mam4xx/water_uptake.hpp>
#endif

namespace exaero::modal {
namespace {
constexpr double pi = 3.14159265358979323846, kb = 1.380649e-23;
using D1 = Kokkos::View<double *>;
using D2 = Kokkos::View<double **, Kokkos::LayoutLeft>;
using D3 = Kokkos::View<double ***, Kokkos::LayoutLeft>;
using D4 = Kokkos::View<double ****, Kokkos::LayoutLeft>;
using I1 = Kokkos::View<int *>;
using H3 = Kokkos::View<const double ***, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryUnmanaged>;
using H2 = Kokkos::View<const double **, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryUnmanaged>;
/// @brief Upload a small immutable array.
template <class T> Kokkos::View<T *> upload(const std::vector<T> &values) {
  Kokkos::View<T *> result("modal_config", values.size());
  auto host = Kokkos::create_mirror_view(result);
  std::copy(values.begin(), values.end(), host.data());
  Kokkos::deep_copy(result, host);
  return result;
}
/// @brief Guard compute-space lifecycle before allocating or launching kernels.
void ready() {
  if (!Kokkos::is_initialized())
    throw std::logic_error(
        "FATAL ERROR: initialize EX-aero before modal compute");
}
/// @brief Validate all results, then commit to caller memory in one pass.
template <class View>
void commit(const View &view, double *output, bool signed_values = false) {
  const auto host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
  for (std::size_t i = 0; i < host.size(); ++i)
    if (!std::isfinite(host.data()[i]) ||
        (!signed_values && host.data()[i] < 0))
      throw std::runtime_error(
          "FATAL ERROR: modal result is nonfinite or negative");
  std::copy_n(host.data(), host.size(), output);
}
/// @brief Standard-normal eight-node Gauss-Hermite abscissae and weights.
KOKKOS_INLINE_FUNCTION double node(int i) {
  constexpr double x[] = {-4.144547186125894, -2.802485861287542,
                          -1.636519042435108, -0.539079811351375,
                          0.539079811351375,  1.636519042435108,
                          2.802485861287542,  4.144547186125894};
  return x[i];
}
KOKKOS_INLINE_FUNCTION double weight(int i) {
  constexpr double w[] = {0.000112614538375, 0.009635220120789,
                          0.117239907661759, 0.373012257679077,
                          0.373012257679077, 0.117239907661759,
                          0.009635220120789, 0.000112614538375};
  return w[i];
}
/// @brief Modal moments and their local differential.
struct Moments {
  double r, wet_r, n, mass, volume, hyg_volume, kappa, wet_volume, density,
      water, area;
};
/// @brief Immutable table device descriptor; flat payload offsets are private.
struct TableDescriptor {
  int kind, nr, nh, nl, r, h, l, ext, sca, asym, back;
};
} // namespace

/// @brief Persistent immutable parameters in the default Kokkos memory space.
class Solver {
public:
  int ns, nm, nf;
  bool coagulation, mam4xx_water;
  double gravity, water_molecular_weight, gas_constant, surface_tension;
  D1 density, kappa, radius, sigma, crystallization_rh, deliquescence_rh,
      transfer_radius, pool;
  I1 allowed, destinations, transfer_source, transfer_destination;
  Kokkos::View<TableDescriptor *> tables;
  /// @brief Construct and upload the adapter's validated snapshots.
  Solver(const Parameters &p, const std::vector<Table> &t)
      : ns(p.density.size()), nm(p.radius.size()), nf(nm * (ns + 1)),
        coagulation(p.coagulation), mam4xx_water(p.mam4xx_water),
        gravity(p.gravity), water_molecular_weight(p.water_molecular_weight),
        gas_constant(p.gas_constant), surface_tension(p.surface_tension),
        density(upload(p.density)), kappa(upload(p.kappa)),
        radius(upload(p.radius)), sigma(upload(p.sigma)),
        crystallization_rh(upload(p.crystallization_rh)),
        deliquescence_rh(upload(p.deliquescence_rh)),
        transfer_radius(upload(p.transfer_radius)), allowed(upload(p.allowed)),
        destinations(upload(p.destinations)),
        transfer_source(upload(p.transfer_source)),
        transfer_destination(upload(p.transfer_destination)),
        tables("modal_tables", t.size()) {
    std::vector<double> values;
    auto ht = Kokkos::create_mirror_view(tables);
    for (std::size_t i = 0; i < t.size(); ++i) {
      const auto append = [&](const auto &v) {
        const int offset = values.size();
        values.insert(values.end(), v.begin(), v.end());
        return offset;
      };
      const auto &a = t[i];
      ht(i) = {a.kind,
               static_cast<int>(a.radius.size()),
               static_cast<int>(a.rh.size()),
               static_cast<int>(a.wavelength.size()),
               append(a.radius),
               append(a.rh),
               append(a.wavelength),
               append(a.extinction),
               append(a.scattering),
               append(a.asymmetry),
               append(a.backscatter)};
    }
    pool = upload(values);
    Kokkos::deep_copy(tables, ht);
  }
};

namespace {
/// @brief Upload host state and environment; temperature and RH always
/// required.
struct Snapshot {
  D3 state;
  D2 temperature, pressure, rh, thickness;
  Snapshot(std::size_t nc, std::size_t nl, int nf, const Environment &e,
           const double *input)
      : state("modal_state", nc, nl, nf), temperature("modal_t", nc, nl),
        pressure("modal_p", nc, nl), rh("modal_rh", nc, nl),
        thickness("modal_dz", nc, nl) {
    Kokkos::deep_copy(state, H3(input, nc, nl, nf));
    Kokkos::deep_copy(temperature, H2(e.temperature, nc, nl));
    Kokkos::deep_copy(rh, H2(e.rh, nc, nl));
    if (e.pressure)
      Kokkos::deep_copy(pressure, H2(e.pressure, nc, nl));
    if (e.thickness)
      Kokkos::deep_copy(thickness, H2(e.thickness, nc, nl));
  }
};
/// @brief Mass/number-constrained spherical lognormal moments, RH capped at
/// .98.
KOKKOS_INLINE_FUNCTION Moments moments(const Solver &p, const Snapshot &d,
                                       int c, int k, int m) {
  Moments a{};
  const int base = m * (p.ns + 1);
  a.n = d.state(c, k, base + p.ns);
  for (int s = 0; s < p.ns; ++s) {
    const double mass = d.state(c, k, base + s), v = mass / p.density(s);
    a.mass += mass;
    a.volume += v;
    a.hyg_volume += v * p.kappa(s);
  }
  const double f =
      Kokkos::min(d.rh(c, k), .98) / (1 - Kokkos::min(d.rh(c, k), .98));
  const double ln = Kokkos::log(p.sigma(m));
  const double factor = 4 * pi / 3 * Kokkos::exp(4.5 * ln * ln);
  a.kappa = a.volume > 0 ? a.hyg_volume / a.volume : 0;
  a.wet_volume = a.volume + a.hyg_volume * f;
  a.r = a.n > 0 ? Kokkos::cbrt(a.volume / a.n / factor) : p.radius(m);
  a.wet_r = a.n > 0 ? Kokkos::cbrt(a.wet_volume / a.n / factor) : a.r;
  a.water = 1000 * a.hyg_volume * f;
#ifdef EXAERO_ENABLE_MAM4XX
  if (p.mam4xx_water && a.n > 0) {
    const double dry_mean_volume = a.volume / a.n;
    const double dry_mean_radius = Kokkos::cbrt(3 * dry_mean_volume / (4 * pi));
    const double rh = Kokkos::min(d.rh(c, k), .98);
    double wet_mean_radius = dry_mean_radius;
    mam4::water_uptake::modal_aero_kohler(dry_mean_radius, a.kappa, rh,
                                          wet_mean_radius);
    wet_mean_radius = Kokkos::max(wet_mean_radius, dry_mean_radius);
    double water_volume = Kokkos::max(
        0., 4 * pi / 3 * wet_mean_radius * wet_mean_radius * wet_mean_radius -
                dry_mean_volume);
    const double ramp =
        rh < p.crystallization_rh(m) ? 0.
        : rh >= p.deliquescence_rh(m)
            ? 1.
            : (rh - p.crystallization_rh(m)) /
                  Kokkos::max(1e-5,
                              p.deliquescence_rh(m) - p.crystallization_rh(m));
    water_volume *= ramp;
    a.wet_volume = a.volume + a.n * water_volume;
    a.water = 1000 * a.n * water_volume;
    a.wet_r = a.r * Kokkos::cbrt(a.wet_volume / a.volume);
  }
#endif
  a.density = a.wet_volume > 0 ? (a.mass + a.water) / a.wet_volume : 1000;
  a.area = 4 * pi * a.n * a.wet_r * a.wet_r * Kokkos::exp(2 * ln * ln);
  return a;
}
/// @brief Pack the nine public properties without including public view
/// headers.
KOKKOS_INLINE_FUNCTION double value(const Moments &a, int q) {
  switch (q) {
  case 0:
    return a.r;
  case 1:
    return a.wet_r;
  case 2:
    return a.n;
  case 3:
    return a.mass;
  case 4:
    return a.volume;
  case 5:
    return a.kappa;
  case 6:
    return a.area;
  case 7:
    return a.density;
  default:
    return a.water;
  }
}
/// @brief Analytic differential with respect to one local mass/number field.
KOKKOS_INLINE_FUNCTION double derivative(const Solver &p, const Snapshot &d,
                                         const Moments &a, int c, int k, int m,
                                         int q, int field) {
  const int s = field - m * (p.ns + 1);
  if (s < 0 || s > p.ns)
    return 0;
  const double dn = s == p.ns ? 1 : 0, dm = s < p.ns ? 1 : 0;
  const double dv = s < p.ns ? 1 / p.density(s) : 0;
  const double dh = s < p.ns ? p.kappa(s) / p.density(s) : 0;
  const double rh = Kokkos::min(d.rh(c, k), .98), f = rh / (1 - rh);
  const double dw = 1000 * f * dh, dwv = dv + f * dh;
  if (q == 2)
    return dn;
  if (q == 3)
    return dm;
  if (q == 4)
    return dv;
  if (q == 8)
    return dw;
  if (a.n == 0)
    return 0;
  switch (q) {
  case 0:
    return a.r / 3 * (dv / a.volume - dn / a.n);
  case 1:
    return a.wet_r / 3 * (dwv / a.wet_volume - dn / a.n);
  case 5:
    return (dh - a.kappa * dv) / a.volume;
  case 6:
    return a.area * (dn / (3 * a.n) + 2 * dwv / (3 * a.wet_volume));
  case 7:
    return (dm + dw - a.density * dwv) / a.wet_volume;
  default:
    return 0;
  }
}
/// @brief Lognormal third-moment fraction below a dry or aerodynamic radius.
KOKKOS_INLINE_FUNCTION double fraction(double r, double sigma, double cut) {
  const double ln = Kokkos::log(sigma);
  return .5 * Kokkos::erfc((Kokkos::log(r / cut) + 3 * ln * ln) /
                           (Kokkos::sqrt(2.) * ln));
}
/// @brief Sutherland viscosity [Pa s].
KOKKOS_INLINE_FUNCTION double viscosity(double t) {
  return 1.716e-5 * Kokkos::pow(t / 273.15, 1.5) * (273.15 + 111) / (t + 111);
}
/// @brief Fuchs flux-matching Brownian kernel, with Cunningham slip [m3/s].
/// See Seinfeld & Pandis (2016), Brownian coagulation transition regime.
KOKKOS_INLINE_FUNCTION double brownian(double r1, double r2, double rho1,
                                       double rho2, double t, double pressure) {
  const double mu = viscosity(t);
  const double mean_path = 6.65e-8 * (t / 293.15) * (101325 / pressure);
  const auto particle = [&](double r, double rho, double &diff, double &speed,
                            double &g) {
    const double kn = mean_path / r;
    const double slip = 1 + kn * (1.257 + .4 * Kokkos::exp(-1.1 / kn));
    diff = kb * t * slip / (6 * pi * mu * r);
    speed = Kokkos::sqrt(8 * kb * t / (pi * (4 * pi / 3 * r * r * r * rho)));
    const double l = 8 * diff / (pi * speed);
    g = l / r < 1e-3 ? l / 2
                     : (Kokkos::pow(2 * r + l, 3) -
                        Kokkos::pow(4 * r * r + l * l, 1.5)) /
                               (6 * r * l) -
                           2 * r;
  };
  double d1, v1, g1, d2, v2, g2;
  particle(r1, rho1, d1, v1, g1);
  particle(r2, rho2, d2, v2, g2);
  const double r = r1 + r2, diff = d1 + d2;
  return 4 * pi * r * diff /
         (r / (r + Kokkos::sqrt(g1 * g1 + g2 * g2)) +
          4 * diff / (r * Kokkos::sqrt(v1 * v1 + v2 * v2)));
}
/// @brief Integrate collision rate and each participant's volume-weighted rate.
KOKKOS_INLINE_FUNCTION void kernels(const Solver &p, const Moments &a,
                                    const Moments &b, int i, int j, double t,
                                    double pressure, double &number,
                                    double &mass_i, double &mass_j) {
  number = mass_i = mass_j = 0;
  const double l1 = Kokkos::log(p.sigma(i)), l2 = Kokkos::log(p.sigma(j));
  const double mean1 = Kokkos::exp(4.5 * l1 * l1),
               mean2 = Kokkos::exp(4.5 * l2 * l2);
  for (int u = 0; u < 8; ++u)
    for (int v = 0; v < 8; ++v) {
      const double f1 = Kokkos::exp(l1 * node(u)),
                   f2 = Kokkos::exp(l2 * node(v));
      const double rate = weight(u) * weight(v) *
                          brownian(a.wet_r * f1, b.wet_r * f2, a.density,
                                   b.density, t, pressure);
      number += rate;
      mass_i += rate * f1 * f1 * f1 / mean1;
      mass_j += rate * f2 * f2 * f2 / mean2;
    }
}
} // namespace

std::shared_ptr<Solver> create(const Parameters &p,
                               const std::vector<Table> &t) {
  ready();
  return std::make_shared<Solver>(p, t);
}

void properties(const Solver &p, std::size_t nc, std::size_t nl,
                const Environment &env, const double *state, double *output,
                double *jac) {
  ready();
  const Snapshot d(nc, nl, p.nf, env, state);
  D4 result("modal_properties", nc, nl, p.nm, 9);
  D4 derivatives("modal_property_derivatives", nc, nl, p.nm,
                 jac ? 9 * p.nf : 0);
  Kokkos::parallel_for(
      "modal_properties", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        for (int m = 0; m < p.nm; ++m) {
          const auto a = moments(p, d, c, k, m);
          for (int q = 0; q < 9; ++q) {
            result(c, k, m, q) = value(a, q);
            if (jac)
              for (int f = 0; f < p.nf; ++f)
                derivatives(c, k, m, q * p.nf + f) =
                    derivative(p, d, a, c, k, m, q, f);
          }
        }
      });
  // Validate both arrays before publishing either one.
  std::vector<double> out(result.size()), dj(derivatives.size());
  commit(result, out.data());
  if (jac)
    commit(derivatives, dj.data(), true);
  std::copy(out.begin(), out.end(), output);
  if (jac)
    std::copy(dj.begin(), dj.end(), jac);
}

void advance(const Solver &p, std::size_t nc, std::size_t nl,
             const Environment &env, double *state, double dt) {
  ready();
  const Snapshot d(nc, nl, p.nf, env, state);
  D3 rates("modal_rates", nc, nl, p.nf);
  I1 failed("modal_integration_failed", nc * nl);
  Kokkos::parallel_for(
      "modal_microphysics", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        double remaining = dt;
        int steps = 0;
        while (p.coagulation && remaining > 0 && steps++ < 100000) {
          for (int f = 0; f < p.nf; ++f)
            rates(c, k, f) = 0;
          for (int i = 0; i < p.nm; ++i) {
            const auto a = moments(p, d, c, k, i);
            if (a.n == 0)
              continue;
            const int bi = i * (p.ns + 1);
            for (int j = i; j < p.nm; ++j) {
              const auto b = moments(p, d, c, k, j);
              if (b.n == 0)
                continue;
              const int bj = j * (p.ns + 1),
                        dest = p.destinations(i * p.nm + j) * (p.ns + 1);
              double kn, ki, kj;
              kernels(p, a, b, i, j, d.temperature(c, k), d.pressure(c, k), kn,
                      ki, kj);
              if (i == j) {
                rates(c, k, bi + p.ns) -= .5 * kn * a.n * a.n;
                continue;
              }
              const double collisions = kn * a.n * b.n;
              rates(c, k, bi + p.ns) -= collisions;
              rates(c, k, bj + p.ns) -= collisions;
              rates(c, k, dest + p.ns) += collisions;
              for (int s = 0; s < p.ns; ++s) {
                const double mi = d.state(c, k, bi + s) * b.n * ki;
                const double mj = d.state(c, k, bj + s) * a.n * kj;
                rates(c, k, bi + s) -= mi;
                rates(c, k, bj + s) -= mj;
                rates(c, k, dest + s) += mi + mj;
              }
            }
          }
          double step = remaining;
          for (int f = 0; f < p.nf; ++f)
            if (rates(c, k, f) < 0)
              step =
                  Kokkos::min(step, .1 * d.state(c, k, f) / (-rates(c, k, f)));
          if (!(step > 0) || !Kokkos::isfinite(step)) {
            failed(point) = 1;
            break;
          }
          for (int f = 0; f < p.nf; ++f)
            d.state(c, k, f) += step * rates(c, k, f);
          const double next = remaining - step;
          if (next == remaining) {
            failed(point) = 1;
            break;
          }
          remaining = next;
        }
        if (p.coagulation && remaining > 0)
          failed(point) = 1;
        for (int route = 0;
             route < static_cast<int>(p.transfer_radius.extent(0)); ++route) {
          const int i = p.transfer_source(route),
                    j = p.transfer_destination(route);
          const auto a = moments(p, d, c, k, i);
          if (a.n == 0 || a.r <= p.transfer_radius(route))
            continue;
          const double ln = Kokkos::log(p.sigma(i)),
                       x = Kokkos::log(p.transfer_radius(route) / a.r);
          double fn = .5 * Kokkos::erfc(x / (Kokkos::sqrt(2.) * ln));
          double fm =
              .5 * Kokkos::erfc((x - 3 * ln * ln) / (Kokkos::sqrt(2.) * ln));
          // At numerical saturation move the whole mode, avoiding ghost
          // particles.
          if (fn == 1 || fm == 1)
            fn = fm = 1;
          const int bi = i * (p.ns + 1), bj = j * (p.ns + 1);
          for (int s = 0; s <= p.ns; ++s) {
            const double moved = d.state(c, k, bi + s) * (s == p.ns ? fn : fm);
            d.state(c, k, bi + s) -= moved;
            d.state(c, k, bj + s) += moved;
          }
        }
      });
  const auto hf =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, failed);
  for (std::size_t i = 0; i < hf.size(); ++i)
    if (hf(i))
      throw std::runtime_error(
          "FATAL ERROR: modal coagulation substep limit or underflow");
  commit(d.state, state);
}

void diagnostics(const Solver &p, std::size_t nc, std::size_t nl,
                 const Environment &env, const double *state, double *output) {
  ready();
  const Snapshot d(nc, nl, p.nf, env, state);
  D3 result("modal_diagnostics", nc, nl, diagnostic_indices::NUM_DIAGNOSTICS);
  Kokkos::parallel_for(
      "modal_diagnostics", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        double settle = 0;
        using namespace diagnostic_indices;
        for (int m = 0; m < p.nm; ++m) {
          const auto a = moments(p, d, c, k, m);
          const double ln = Kokkos::log(p.sigma(m));
          const double aero = a.wet_r * Kokkos::sqrt(a.density / 1000);
          result(c, k, MASS_CONCENTRATION) += a.mass;
          result(c, k, NUMBER_CONCENTRATION) += a.n;
          result(c, k, SURFACE_AREA_DENSITY) += a.area;
          result(c, k, AEROSOL_LIQUID_WATER) += a.water;
          result(c, k, PM1_CONCENTRATION) +=
              a.mass * fraction(a.r, p.sigma(m), .5e-6);
          result(c, k, PM2_5_CONCENTRATION) +=
              a.mass * fraction(aero, p.sigma(m), 1.25e-6);
          result(c, k, PM10_CONCENTRATION) +=
              a.mass * fraction(aero, p.sigma(m), 5e-6);
          settle += a.mass * 2 * a.density * p.gravity * a.wet_r * a.wet_r *
                    Kokkos::exp(8 * ln * ln) /
                    (9 * viscosity(d.temperature(c, k)));
        }
        result(c, k, GRAVITATIONAL_SETTLING_VELOCITY) =
            result(c, k, MASS_CONCENTRATION) > 0
                ? settle / result(c, k, MASS_CONCENTRATION)
                : 0;
      });
  Kokkos::parallel_for(
      "modal_columns", nc, KOKKOS_LAMBDA(int c) {
        using namespace diagnostic_indices;
        for (std::size_t k = 0; k < nl; ++k) {
          result(c, 0, COLUMN_MASS) +=
              result(c, k, MASS_CONCENTRATION) * d.thickness(c, k);
          result(c, 0, COLUMN_PM2_5_MASS) +=
              result(c, k, PM2_5_CONCENTRATION) * d.thickness(c, k);
        }
        result(c, 0, SURFACE_MASS) = result(c, 0, MASS_CONCENTRATION);
        result(c, 0, SURFACE_PM2_5_MASS) = result(c, 0, PM2_5_CONCENTRATION);
      });
  commit(result, output);
}

void ccn(const Solver &p, std::size_t nc, std::size_t nl,
         const Environment &env, const double *state, const double *axis,
         std::size_t count, double *output) {
  ready();
  const Snapshot d(nc, nl, p.nf, env, state);
  const auto ss = upload(std::vector<double>(axis, axis + count));
  D3 result("modal_ccn", nc, nl, count);
  Kokkos::parallel_for(
      "modal_ccn", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        const double kelvin = 4 * p.surface_tension * p.water_molecular_weight /
                              (p.gas_constant * d.temperature(c, k) * 1000);
        for (int m = 0; m < p.nm; ++m) {
          const auto a = moments(p, d, c, k, m);
          if (a.n == 0 || a.kappa == 0)
            continue;
          const double critical =
              Kokkos::sqrt(4 * kelvin * kelvin * kelvin /
                           (27 * a.kappa * Kokkos::pow(2 * a.r, 3)));
          for (std::size_t s = 0; s < count; ++s)
            if (ss(s) > 0)
              result(c, k, s) +=
                  a.n * .5 *
                  Kokkos::erfc(
                      2 * Kokkos::log(critical / ss(s)) /
                      (3 * Kokkos::sqrt(2.) * Kokkos::log(p.sigma(m))));
        }
      });
  commit(result, output);
}

void uptake(const Solver &p, std::size_t nc, std::size_t nl,
            const Environment &env, const double *state,
            const std::vector<Vapor> &vapors, double *output, double *jac) {
  ready();
  const Snapshot d(nc, nl, p.nf, env, state);
  const auto v = upload(vapors);
  const int nv = vapors.size();
  D4 result("modal_uptake", nc, nl, p.nm, nv);
  D4 derivatives("modal_uptake_jac", nc, nl, p.nm, jac ? nv * p.nf : 0);
  Kokkos::parallel_for(
      "modal_uptake", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        for (int m = 0; m < p.nm; ++m) {
          const auto a = moments(p, d, c, k, m);
          if (a.n == 0)
            continue;
          for (int gas = 0; gas < nv; ++gas) {
            const double speed =
                Kokkos::sqrt(8 * p.gas_constant * d.temperature(c, k) /
                             (pi * v(gas).molecular_weight));
            const double path = 3 * v(gas).diffusivity / speed;
            const double alpha = v(gas).accommodation;
            double integral = 0, dr = 0;
            for (int q = 0; q < 8; ++q) {
              const double scale =
                  Kokkos::exp(Kokkos::log(p.sigma(m)) * node(q));
              const double r = a.r * scale, kn = path / r;
              const double aa = 4 / (3 * alpha), bb = .377 + aa;
              const double den = 1 + bb * kn + aa * kn * kn;
              const double f = (1 + kn) / den;
              const double df =
                  (den - (1 + kn) * (bb + 2 * aa * kn)) / (den * den);
              integral += weight(q) * 4 * pi * v(gas).diffusivity * r * f;
              dr += weight(q) * 4 * pi * v(gas).diffusivity * scale *
                    (f - kn * df);
            }
            result(c, k, m, gas) = a.n * integral;
            if (jac)
              for (int field = 0; field < p.nf; ++field)
                derivatives(c, k, m, gas * p.nf + field) =
                    (field == m * (p.ns + 1) + p.ns ? integral : 0) +
                    a.n * dr * derivative(p, d, a, c, k, m, 0, field);
          }
        }
      });
  std::vector<double> out(result.size()), dj(derivatives.size());
  commit(result, out.data());
  if (jac)
    commit(derivatives, dj.data(), true);
  std::copy(out.begin(), out.end(), output);
  if (jac)
    std::copy(dj.begin(), dj.end(), jac);
}

namespace {
/// @brief Locate within a validated increasing axis; no extrapolation.
KOKKOS_INLINE_FUNCTION int locate(const D1 &v, int offset, int count, double x,
                                  double &w, bool logarithmic = false) {
  const double tolerance =
      64 * 2.220446049250313e-16 *
      Kokkos::max(Kokkos::abs(v(offset)), Kokkos::abs(v(offset + count - 1)));
  if (x < v(offset) - tolerance || x > v(offset + count - 1) + tolerance)
    return -1;
  x = Kokkos::max(v(offset), Kokkos::min(v(offset + count - 1), x));
  if (count == 1) {
    w = 0;
    return 0;
  }
  int lo = 0;
  while (lo < count - 2 && x > v(offset + lo + 1))
    ++lo;
  const double a = v(offset + lo), b = v(offset + lo + 1);
  w = logarithmic ? Kokkos::log(x / a) / Kokkos::log(b / a) : (x - a) / (b - a);
  return lo;
}
/// @brief Trilinear table evaluation; g is interpolated as scattering*g.
KOKKOS_INLINE_FUNCTION double lookup(const Solver &p, const TableDescriptor &t,
                                     int field, int r, int h, int l, double wr,
                                     double wh, double wl) {
  const int offset = field == 0   ? t.ext
                     : field == 1 ? t.sca
                     : field == 2 ? t.asym
                                  : t.back;
  double sum = 0;
  for (int ir = 0; ir < (t.nr > 1 ? 2 : 1); ++ir)
    for (int ih = 0; ih < (t.nh > 1 ? 2 : 1); ++ih)
      for (int il = 0; il < (t.nl > 1 ? 2 : 1); ++il) {
        const int index = ((r + ir) * t.nh + h + ih) * t.nl + l + il;
        const double w = (t.nr == 1 ? 1
                          : ir      ? wr
                                    : 1 - wr) *
                         (t.nh == 1 ? 1
                          : ih      ? wh
                                    : 1 - wh) *
                         (t.nl == 1 ? 1
                          : il      ? wl
                                    : 1 - wl);
        double val = p.pool(offset + index);
        if (field == 2)
          val *= p.pool(t.sca + index);
        sum += w * val;
      }
  return sum;
}
/// @brief Integrate a linear-in-radius coefficient against lognormal dry mass.
/// Return endpoint basis weights over [lo, min(hi,cut)]. This integrates the
/// interpolated curve exactly, avoiding sparse quadrature across Mie
/// resonances.
KOKKOS_INLINE_FUNCTION void opticalWeights(double rg, double ln, double lo,
                                           double hi, double cut, double &w0,
                                           double &w1) {
  const double end = Kokkos::min(hi, cut);
  w0 = w1 = 0;
  if (end <= lo)
    return;
  const double root = Kokkos::sqrt(2.) * ln;
  const double l = Kokkos::log(lo / rg), u = Kokkos::log(end / rg);
  const double p3 = .5 * (Kokkos::erfc((l - 3 * ln * ln) / root) -
                          Kokkos::erfc((u - 3 * ln * ln) / root));
  const double p4 = rg * Kokkos::exp(3.5 * ln * ln) * .5 *
                    (Kokkos::erfc((l - 4 * ln * ln) / root) -
                     Kokkos::erfc((u - 4 * ln * ln) / root));
  w0 = Kokkos::max(0., (hi * p3 - p4) / (hi - lo));
  w1 = Kokkos::max(0., (p4 - lo * p3) / (hi - lo));
}
} // namespace

void optics(const Solver &p, std::size_t nc, std::size_t nl,
            const Environment &env, const double *state, const double *axis,
            std::size_t count, double *output) {
  ready();
  const Snapshot d(nc, nl, p.nf, env, state);
  const auto bands = upload(std::vector<double>(axis, axis + count));
  Kokkos::View<double *****, Kokkos::LayoutLeft> result(
      "modal_optical_integrals", nc, nl, p.nm, count, 9);
  I1 failed("modal_optical_domain_error", nc * nl);
  Kokkos::parallel_for(
      "modal_table_optics", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        for (int m = 0; m < p.nm; ++m) {
          const auto a = moments(p, d, c, k, m);
          if (a.n == 0)
            continue;
          const double ln = Kokkos::log(p.sigma(m));
          for (int s = 0; s < p.ns; ++s) {
            const double mass = d.state(c, k, m * (p.ns + 1) + s);
            if (mass == 0)
              continue;
            const auto t = p.tables(m * p.ns + s);
            if (t.kind == 0)
              continue; // provider-owned mode
            for (std::size_t b = 0; b < count; ++b) {
              double wh = 0, wl = 0;
              const int h = locate(p.pool, t.h, t.nh, d.rh(c, k), wh);
              const int l = locate(p.pool, t.l, t.nl, bands(b), wl, true);
              if (h < 0 || l < 0) {
                failed(point) = 1;
                continue;
              }
              const int nq = t.kind == 1 ? t.nr - 1 : 1;
              if (t.kind == 1) {
                const double coverage =
                    fraction(a.r, p.sigma(m), p.pool(t.r + t.nr - 1)) -
                    fraction(a.r, p.sigma(m), p.pool(t.r));
                if (coverage < 1 - 1e-8) {
                  failed(point) = 1;
                  continue;
                }
              }
              for (int q = 0; q < nq; ++q) {
                double wr = 0;
                const int r =
                    t.kind == 1 ? q : locate(p.pool, t.r, t.nr, a.r, wr);
                if (r < 0) {
                  failed(point) = 1;
                  continue;
                }
                double w0 = 1, w1 = 0, fw0 = 0, fw1 = 0, pw0 = 0, pw1 = 0;
                const double aero =
                    a.wet_r / a.r * Kokkos::sqrt(a.density / 1000);
                if (t.kind == 1) {
                  const double lo = p.pool(t.r + r), hi = p.pool(t.r + r + 1);
                  opticalWeights(a.r, ln, lo, hi, hi, w0, w1);
                  opticalWeights(a.r, ln, lo, hi, .5e-6, fw0, fw1);
                  opticalWeights(a.r, ln, lo, hi, 1.25e-6 / aero, pw0, pw1);
                } else {
                  fw0 = fraction(a.r, p.sigma(m), .5e-6);
                  pw0 = fraction(a.r * aero, p.sigma(m), 1.25e-6);
                }
                const auto integral = [&](int field, double left,
                                          double right) {
                  if (t.kind != 1)
                    return mass * left *
                           lookup(p, t, field, r, h, l, wr, wh, wl);
                  return mass *
                         (left * lookup(p, t, field, r, h, l, 0, wh, wl) +
                          right * lookup(p, t, field, r, h, l, 1, wh, wl));
                };
                const double ext = integral(0, w0, w1),
                             sca = integral(1, w0, w1);
                result(c, k, m, b, 0) += ext;
                result(c, k, m, b, 1) += sca;
                result(c, k, m, b, 2) += integral(2, w0, w1);
                // Source bbck is the lidar (180 degree) coefficient. Use the
                // HG hemispheric integral for the separate hemispheric
                // diagnostic, normalized by 2*pi to retain the public
                // per-steradian unit.
                const double weighted_g = integral(2, w0, w1);
                const double g = sca > 0 ? weighted_g / sca : 0;
                const double hemi =
                    g >= 1    ? 0
                    : g <= -1 ? 1
                    : Kokkos::abs(g) < 1e-6
                        ? .5
                        : (1 - g * g) / (2 * g) *
                              (1 / Kokkos::sqrt(1 + g * g) - 1 / (1 + g));
                result(c, k, m, b, 3) += sca * hemi / (2 * pi);
                result(c, k, m, b, 4) += integral(3, w0, w1);
                result(c, k, m, b, 5) += integral(0, fw0, fw1);
                result(c, k, m, b, 6) += integral(1, fw0, fw1);
                result(c, k, m, b, 7) += integral(0, pw0, pw1);
                result(c, k, m, b, 8) += integral(1, pw0, pw1);
              }
            }
          }
          for (std::size_t b = 0; b < count; ++b)
            result(c, k, m, b, 2) =
                result(c, k, m, b, 1) > 0
                    ? result(c, k, m, b, 2) / result(c, k, m, b, 1)
                    : 0;
        }
      });
  const auto hf =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, failed);
  for (std::size_t i = 0; i < hf.size(); ++i)
    if (hf(i))
      throw std::runtime_error("FATAL ERROR: modal optical table does not "
                               "cover radius/RH/wavelength query");
  commit(result, output, true);
}
} // namespace exaero::modal
