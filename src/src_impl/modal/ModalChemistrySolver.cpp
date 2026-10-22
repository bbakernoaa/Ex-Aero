/// @file ModalChemistrySolver.cpp
/// @brief Device-parallel balanced exchange and coating-shell aging RHS.
/// @author EX-aero contributors
/// @details History: 2026-10. Host column-major buffers; no file
/// I/O/integration.
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <modal/ModalChemistrySolver.hpp>
#include <stdexcept>

namespace exaero::modal {
namespace {
using D3 = Kokkos::View<double ***, Kokkos::LayoutLeft>;
using D4 = Kokkos::View<double ****, Kokkos::LayoutLeft>;
using H3 = Kokkos::View<const double ***, Kokkos::LayoutLeft, Kokkos::HostSpace,
                        Kokkos::MemoryUnmanaged>;
using H4 = Kokkos::View<const double ****, Kokkos::LayoutLeft,
                        Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
/// @brief Guard lifecycle before staging host data.
void ready() {
  if (!Kokkos::is_initialized())
    throw std::logic_error(
        "FATAL ERROR: initialize EX-aero before modal chemistry");
}
/// @brief Small immutable metadata uploaded to the default memory space.
template <class T> Kokkos::View<T *> upload(const std::vector<T> &data) {
  Kokkos::View<T *> out("chemistry_config", data.size());
  auto host = Kokkos::create_mirror_view(out);
  std::copy(data.begin(), data.end(), host.data());
  Kokkos::deep_copy(out, host);
  return out;
}
/// @brief Stage a validated host field slab.
D3 stage(const double *input, std::size_t nc, std::size_t nl, int nf) {
  D3 out("chemistry_input", nc, nl, nf);
  Kokkos::deep_copy(out, H3(input, nc, nl, nf));
  return out;
}
/// @brief Stage optional signed derivatives; absent arrays remain zero.
D4 derivatives(const double *input, std::size_t nc, std::size_t nl, int nr,
               int ni) {
  D4 out("chemistry_derivatives", nc, nl, input ? nr : 0, input ? ni : 0);
  if (input)
    Kokkos::deep_copy(out, H4(input, nc, nl, nr, ni));
  return out;
}
/// @brief Validate every signed output before publishing any buffer.
template <class View> std::vector<double> checked(const View &out) {
  const auto host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, out);
  std::vector<double> result(host.size());
  for (std::size_t i = 0; i < host.size(); ++i) {
    if (!std::isfinite(host.data()[i]))
      throw std::runtime_error(
          "FATAL ERROR: modal chemistry result is nonfinite");
    result[i] = host.data()[i];
  }
  return result;
}
} // namespace

void exchange(std::size_t nc, std::size_t nl, int ns, int nm,
              const std::vector<int> &binding,
              const std::vector<int> &permitted, int ni, const double *state,
              const double *gas, const double *equilibrium,
              const double *uptake, const double *uptake_jac,
              const double *equilibrium_jac, double *aerosol_rates,
              double *gas_rates, double *jac) {
  ready();
  const int nv = binding.size(), nf = nm * (ns + 1), nq = nf + nv;
  const auto species = upload(binding), allowed = upload(permitted);
  const auto a = stage(state, nc, nl, nf), g = stage(gas, nc, nl, nv);
  D4 ce("exchange_equilibrium", nc, nl, nm, nv),
      ku("exchange_uptake", nc, nl, nm, nv);
  Kokkos::deep_copy(ce, H4(equilibrium, nc, nl, nm, nv));
  Kokkos::deep_copy(ku, H4(uptake, nc, nl, nm, nv));
  const auto dk = derivatives(uptake_jac, nc, nl, nm, nv * nf);
  const auto de = derivatives(equilibrium_jac, nc, nl, nm * nv, ni);
  D3 ar("exchange_aerosol_rhs", nc, nl, nf), gr("exchange_gas_rhs", nc, nl, nv);
  D4 jr("exchange_jacobian", nc, nl, jac ? nq : 0, jac ? ni : 0);
  const bool has_equilibrium_jac = equilibrium_jac != nullptr;
  Kokkos::parallel_for(
      "modal_phase_exchange", nc * nl, KOKKOS_LAMBDA(int point) {
        const int c = point % nc, k = point / nc;
        for (int m = 0; m < nm; ++m)
          for (int v = 0; v < nv; ++v) {
            const int s = species(v), field = m * (ns + 1) + s;
            if (!allowed(m * ns + s))
              continue;
            const double delta = g(c, k, v) - ce(c, k, m, v),
                         flux = ku(c, k, m, v) * delta;
            if (delta < 0 && a(c, k, field) == 0)
              continue;
            ar(c, k, field) += flux;
            gr(c, k, v) -= flux;
            if (jac)
              for (int column = 0; column < ni; ++column) {
                double value =
                    column < nf ? dk(c, k, m, v * nf + column) * delta : 0;
                if (column == nf + v)
                  value += ku(c, k, m, v);
                if (has_equilibrium_jac)
                  value -= ku(c, k, m, v) * de(c, k, m * nv + v, column);
                jr(c, k, field, column) += value;
                jr(c, k, nf + v, column) -= value;
              }
          }
      });
  const auto ah = checked(ar), gh = checked(gr), jh = checked(jr);
  std::copy(ah.begin(), ah.end(), aerosol_rates);
  std::copy(gh.begin(), gh.end(), gas_rates);
  if (jac)
    std::copy(jh.begin(), jh.end(), jac);
}

void aging(std::size_t nc, std::size_t nl, const AgingParameters &p, int ni,
           const double *state, const double *condensation,
           const double *coagulation, const double *condensation_jac,
           int condensation_rows, const double *coagulation_jac,
           int coagulation_rows, double *rates, double *jac) {
  ready();
  const int ns = p.density.size(), nm = p.sigma.size(), nf = nm * (ns + 1),
            nr = p.source.size();
  const auto density = upload(p.density), sigma = upload(p.sigma),
             thickness = upload(p.thickness);
  const auto source = upload(p.source), destination = upload(p.destination),
             coating = upload(p.coating);
  const auto a = stage(state, nc, nl, nf),
             cond = stage(condensation, nc, nl, nf),
             coag = stage(coagulation, nc, nl, nf);
  const auto dc = derivatives(condensation_jac, nc, nl, condensation_rows, ni);
  const auto dg = derivatives(coagulation_jac, nc, nl, coagulation_rows, ni);
  D3 out("aging_rhs", nc, nl, nf);
  D4 dj("aging_jacobian", nc, nl, jac ? nf : 0, jac ? ni : 0);
  Kokkos::View<int *> failed("aging_domain_failure", nc * nl);
  const bool has_dc = condensation_jac != nullptr,
             has_dg = coagulation_jac != nullptr;
  Kokkos::parallel_for(
      "modal_coating_aging", nc * nl, KOKKOS_LAMBDA(int point) {
        constexpr double pi = 3.14159265358979323846;
        const int c = point % nc, k = point / nc;
        for (int route = 0; route < nr; ++route) {
          const int m = source(route), base = m * (ns + 1),
                    dest = destination(route) * (ns + 1);
          const double number = a(c, k, base + ns);
          if (number == 0)
            continue;
          double volume = 0, budget = 0;
          for (int s = 0; s < ns; ++s) {
            volume += a(c, k, base + s) / density(s);
            if (coating(route * ns + s))
              budget +=
                  (cond(c, k, base + s) + coag(c, k, base + s)) / density(s);
          }
          if (budget <= 0)
            continue;
          const double l = Kokkos::log(sigma(m)), h = thickness(route);
          const double radius = Kokkos::cbrt(volume / number / (4 * pi / 3) /
                                             Kokkos::exp(4.5 * l * l));
          const double e2 = Kokkos::exp(2 * l * l),
                       e1 = Kokkos::exp(.5 * l * l);
          const double shell_moment =
              radius * radius * e2 * h + radius * e1 * h * h + h * h * h / 3;
          const double shell = 4 * pi * number * shell_moment;
          if (!(radius > 0) || !(shell > 0) || !Kokkos::isfinite(shell) ||
              !Kokkos::isfinite(budget)) {
            failed(point) = 1;
            continue;
          }
          const double frequency = budget / shell;
          if (!Kokkos::isfinite(frequency)) {
            failed(point) = 1;
            continue;
          }
          for (int s = 0; s <= ns; ++s) {
            const double value = frequency * a(c, k, base + s);
            out(c, k, base + s) -= value;
            out(c, k, dest + s) += value;
          }
          if (jac)
            for (int column = 0; column < ni; ++column) {
              const int local = column - base;
              const double dn = local == ns ? 1. : 0.,
                           dv = local >= 0 && local < ns ? 1 / density(local)
                                                         : 0.;
              const double dr = radius / 3 * (dv / volume - dn / number);
              const double dshell =
                  4 * pi *
                  (dn * shell_moment +
                   number * (2 * radius * e2 * h + e1 * h * h) * dr);
              double dbudget = 0;
              for (int s = 0; s < ns; ++s)
                if (coating(route * ns + s)) {
                  if (has_dc)
                    dbudget += dc(c, k, base + s, column) / density(s);
                  if (has_dg)
                    dbudget += dg(c, k, base + s, column) / density(s);
                }
              const double dfrequency =
                  dbudget / shell - frequency * dshell / shell;
              for (int s = 0; s <= ns; ++s) {
                const double value = dfrequency * a(c, k, base + s) +
                                     (column == base + s ? frequency : 0.);
                dj(c, k, base + s, column) -= value;
                dj(c, k, dest + s, column) += value;
              }
            }
        }
      });
  const auto failure =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, failed);
  for (std::size_t i = 0; i < failure.size(); ++i)
    if (failure(i))
      throw std::runtime_error(
          "FATAL ERROR: unrepresentable modal aging shell/budget");
  const auto rh = checked(out), jh = checked(dj);
  std::copy(rh.begin(), rh.end(), rates);
  if (jac)
    std::copy(jh.begin(), jh.end(), jac);
}
} // namespace exaero::modal
