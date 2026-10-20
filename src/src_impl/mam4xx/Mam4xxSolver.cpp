/// @file Mam4xxSolver.cpp
/// @brief Kokkos-only implementation using pinned upstream MAM4xx kernels.
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <exaero/AerosolIndices.hpp>
#include <exaero/Mam4xxConfig.hpp>
#include <mam4xx/Mam4xxSolver.hpp>
#include <mam4xx/coagulation.hpp>
#include <mam4xx/conversions.hpp>
#include <mam4xx/water_uptake.hpp>
#include <stdexcept>
#include <type_traits>

namespace exaero {
namespace {
static_assert(std::is_same_v<mam4::Real, double>,
              "EX-aero requires double-precision MAM4xx");
static_assert(mam4::AeroConfig::num_modes() == mam4xx_layout::modes);
static_assert(mam4::AeroConfig::num_aerosol_ids() == mam4xx_layout::species);
using Device3 = Kokkos::View<double ***, Kokkos::LayoutLeft>;
using Device2 = Kokkos::View<double **, Kokkos::LayoutLeft>;
using Host3 =
    Kokkos::View<const double ***, Kokkos::LayoutLeft, Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using Host2 =
    Kokkos::View<const double **, Kokkos::LayoutLeft, Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

/// @brief Own the device snapshot; public arrays are always host memory.
struct Snapshot {
  Device3 state;
  Device2 temperature, pressure, rh, thickness;
  mam4::AeroSpeciesView species;
  /// @brief Upload the immutable host boundary for this call.
  Snapshot(std::size_t cells, std::size_t levels, const Mam4xxEnvironment &env,
           const double *input)
      : state("mam4xx_state", cells, levels, mam4xx_layout::fields),
        temperature("mam4xx_temperature", cells, levels),
        pressure("mam4xx_pressure", cells, levels),
        rh("mam4xx_rh", cells, levels),
        thickness("mam4xx_thickness", cells, levels),
        species(mam4::aero_species_on_device(mam4::default_aero_species())) {
    Kokkos::deep_copy(state,
                      Host3(input, cells, levels, mam4xx_layout::fields));
    Kokkos::deep_copy(temperature, Host2(env.temperature, cells, levels));
    // The legacy C CCN surface supplies only temperature and RH; pressure
    // and thickness are not read by the CCN kernel.
    if (env.pressure)
      Kokkos::deep_copy(pressure, Host2(env.pressure, cells, levels));
    Kokkos::deep_copy(rh, Host2(env.rh, cells, levels));
    if (env.thickness)
      Kokkos::deep_copy(thickness, Host2(env.thickness, cells, levels));
  }
};

/// @brief Modal integral properties; all extensive quantities are per m3 of
/// air.
struct ModalProperties {
  double number[4], mass[4], volume[4], kappa[4], diameter[4], wet_diameter[4];
  double wet_density[4], water[4];
};

/// @brief Derive internally mixed sizes from mass AND prognostic number.
/// Uses upstream material properties and water-uptake hysteresis for each
/// phase.
KOKKOS_INLINE_FUNCTION ModalProperties properties(const Snapshot &data,
                                                  std::size_t c, std::size_t k,
                                                  bool cloudborne) {
  ModalProperties p{};
  double dry_radius[4], dry_volume[4], dry_mass[4], crystal[4],
      deliquescence[4];
  double wet_radius[4], wet_volume[4], water_volume[4];
  for (int m = 0; m < 4; ++m) {
    p.number[m] = data.state(c, k, mam4xx_layout::number(m, cloudborne));
    for (int s = 0; s < 7; ++s) {
      const double mass =
          data.state(c, k, mam4xx_layout::mass(m, s, cloudborne));
      const double volume = mass / data.species(s).density;
      p.mass[m] += mass;
      p.volume[m] += volume;
      p.kappa[m] += volume * data.species(s).hygroscopicity;
    }
    p.kappa[m] = p.volume[m] > 0 ? p.kappa[m] / p.volume[m] : 0;
    const auto mode = mam4::modes(m);
    dry_volume[m] = p.number[m] > 0
                        ? p.volume[m] / p.number[m]
                        : mam4::conversions::mean_particle_volume_from_diameter(
                              mode.nom_diameter, mode.mean_std_dev);
    p.diameter[m] = mam4::conversions::mean_particle_diameter_from_volume(
        dry_volume[m], mode.mean_std_dev);
    dry_radius[m] = Kokkos::cbrt(3 * dry_volume[m] / (4 * mam4::Constants::pi));
    dry_mass[m] = p.number[m] > 0 ? p.mass[m] / p.number[m] : 0;
    crystal[m] = mode.crystallization_pt;
    deliquescence[m] = mode.deliquescence_pt;
  }
  mam4::water_uptake::modal_aero_water_uptake_wetaer(
      crystal, deliquescence, p.diameter, dry_radius, p.kappa,
      Kokkos::min(data.rh(c, k), 0.98), p.number, dry_volume, wet_radius,
      wet_volume, water_volume, p.wet_diameter, p.water);
  for (int m = 0; m < 4; ++m)
    p.wet_density[m] =
        p.number[m] > 0
            ? (dry_mass[m] + mam4::Constants::density_h2o * water_volume[m]) /
                  wet_volume[m]
            : 1000;
  return p;
}

/// @brief Validate device results before changing any caller-owned output.
void copyBack(const Device3 &view, double *output) {
  const auto host =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, view);
  for (std::size_t i = 0; i < host.size(); ++i)
    if (!std::isfinite(host.data()[i]) || host.data()[i] < 0)
      throw std::runtime_error(
          "FATAL ERROR: MAM4xx produced a nonfinite or negative result");
  std::copy_n(host.data(), host.size(), output);
}

/// @brief Lognormal mass fraction below a diameter, with the third-moment
/// shift.
KOKKOS_INLINE_FUNCTION double massFraction(double diameter, double sigma,
                                           double cut) {
  const double lns = Kokkos::log(sigma);
  return 0.5 * Kokkos::erfc((Kokkos::log(diameter / cut) + 3 * lns * lns) /
                            (Kokkos::sqrt(2.0) * lns));
}
} // namespace

void mam4xxMicrophysics(std::size_t cells, std::size_t levels,
                        const Mam4xxEnvironment &env, double *state,
                        double dt) {
  if (!Kokkos::is_initialized())
    throw std::logic_error("FATAL ERROR: initialize the EX-aero environment "
                           "before MAM4xx compute");
  const Snapshot data(cells, levels, env, state);
  // Read-only input and separate result avoid races in modal property queries.
  Device3 result("mam4xx_advanced_state", cells, levels, mam4xx_layout::fields);
  Kokkos::deep_copy(result, data.state);
  Kokkos::parallel_for(
      "mam4xx_coagulation", Kokkos::RangePolicy<>(0, cells * levels),
      KOKKOS_LAMBDA(int point) {
        const int c = point % cells, k = point / cells;
        const auto p = properties(data, c, k, false);
        const double air_moles = data.pressure(c, k) / (mam4::Constants::r_gas *
                                                        data.temperature(c, k));
        double number[4], mass[7][4],
            transfer[7][mam4::AeroConfig::max_agepair()]{};
        for (int m = 0; m < 4; ++m) {
          // Upstream coefficients are multiplied by mol-air/m3: number must be
          // #/mol-air.
          number[m] = p.number[m] / air_moles;
          for (int s = 0; s < 7; ++s)
            mass[s][m] = data.state(c, k, mam4xx_layout::mass(m, s));
        }
        // Mass transfer is linear in mass; concentration units are retained
        // exactly.
        mam4::coagulation::mam_coag_1subarea(
            dt, data.temperature(c, k), data.pressure(c, k), air_moles,
            p.wet_diameter, p.wet_density, number, mass, transfer);
        for (int m = 0; m < 4; ++m) {
          result(c, k, mam4xx_layout::number(m)) = number[m] * air_moles;
          for (int s = 0; s < 7; ++s)
            result(c, k, mam4xx_layout::mass(m, s)) = mass[s][m];
        }
      });
  copyBack(result, state);
}

void mam4xxDiagnostics(std::size_t cells, std::size_t levels,
                       const Mam4xxEnvironment &env, const double *state,
                       double *output) {
  if (!Kokkos::is_initialized())
    throw std::logic_error("FATAL ERROR: initialize the EX-aero environment "
                           "before MAM4xx compute");
  const Snapshot data(cells, levels, env, state);
  Device3 result("mam4xx_diagnostics", cells, levels,
                 diagnostic_indices::NUM_DIAGNOSTICS);
  Kokkos::parallel_for(
      "mam4xx_diagnostics", Kokkos::RangePolicy<>(0, cells * levels),
      KOKKOS_LAMBDA(int point) {
        const int c = point % cells, k = point / cells;
        double mass = 0, number = 0, area = 0, water = 0, pm1 = 0, pm25 = 0,
               pm10 = 0, settling = 0;
        for (int phase = 0; phase < 2; ++phase) {
          const auto p = properties(data, c, k, phase);
          for (int m = 0; m < 4; ++m) {
            const double sigma = mam4::modes(m).mean_std_dev;
            const double ln2 = Kokkos::pow(Kokkos::log(sigma), 2);
            const double aerodynamic =
                p.wet_diameter[m] * Kokkos::sqrt(p.wet_density[m] / 1000.0);
            mass += p.mass[m];
            number += p.number[m];
            water += p.water[m];
            area += mam4::Constants::pi * p.number[m] * p.wet_diameter[m] *
                    p.wet_diameter[m] * Kokkos::exp(2 * ln2);
            pm1 += p.mass[m] * massFraction(p.diameter[m], sigma, 1e-6);
            pm25 += p.mass[m] * massFraction(aerodynamic, sigma, 2.5e-6);
            pm10 += p.mass[m] * massFraction(aerodynamic, sigma, 10e-6);
            // Stokes approximation with Sutherland viscosity and mass-weighted
            // second moment.
            const double viscosity =
                1.716e-5 * Kokkos::pow(data.temperature(c, k) / 273.15, 1.5) *
                (273.15 + 111) / (data.temperature(c, k) + 111);
            settling += p.mass[m] * p.wet_density[m] *
                        mam4::Constants::gravity * p.wet_diameter[m] *
                        p.wet_diameter[m] * Kokkos::exp(8 * ln2) /
                        (18 * viscosity);
          }
        }
        using namespace diagnostic_indices;
        result(c, k, MASS_CONCENTRATION) = mass;
        result(c, k, NUMBER_CONCENTRATION) = number;
        result(c, k, SURFACE_AREA_DENSITY) = area;
        result(c, k, AEROSOL_LIQUID_WATER) = water;
        result(c, k, PM1_CONCENTRATION) = pm1;
        result(c, k, PM2_5_CONCENTRATION) = pm25;
        result(c, k, PM10_CONCENTRATION) = pm10;
        result(c, k, GRAVITATIONAL_SETTLING_VELOCITY) =
            mass > 0 ? settling / mass : 0;
      });
  Kokkos::parallel_for(
      "mam4xx_column_integrals", Kokkos::RangePolicy<>(0, cells),
      KOKKOS_LAMBDA(int c) {
        using namespace diagnostic_indices;
        double column = 0, pm25 = 0;
        for (std::size_t k = 0; k < levels; ++k) {
          column += result(c, k, MASS_CONCENTRATION) * data.thickness(c, k);
          pm25 += result(c, k, PM2_5_CONCENTRATION) * data.thickness(c, k);
        }
        result(c, 0, COLUMN_MASS) = column;
        result(c, 0, COLUMN_PM2_5_MASS) = pm25;
        result(c, 0, SURFACE_MASS) = result(c, 0, MASS_CONCENTRATION);
        result(c, 0, SURFACE_PM2_5_MASS) = result(c, 0, PM2_5_CONCENTRATION);
      });
  copyBack(result, output);
}

void mam4xxCCN(std::size_t cells, std::size_t levels,
               const Mam4xxEnvironment &env, const double *state,
               const double *supersaturations, std::size_t count,
               double *output) {
  if (!Kokkos::is_initialized())
    throw std::logic_error("FATAL ERROR: initialize the EX-aero environment "
                           "before MAM4xx compute");
  const Snapshot data(cells, levels, env, state);
  Kokkos::View<double *> axis("mam4xx_supersaturations", count);
  const Kokkos::View<const double *, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      host_axis(supersaturations, count);
  Kokkos::deep_copy(axis, host_axis);
  Device3 result("mam4xx_ccn", cells, levels, count);
  Kokkos::parallel_for(
      "mam4xx_ccn", Kokkos::RangePolicy<>(0, cells * levels),
      KOKKOS_LAMBDA(int point) {
        const int c = point % cells, k = point / cells;
        // Modal kappa-Kohler integral as in MAM4 ndrop::ccncalc, generalized to
        // the caller's axis.
        const double a = 4 * 0.076 * mam4::Constants::molec_weight_h2o /
                         (mam4::Constants::r_gas * data.temperature(c, k) *
                          mam4::Constants::density_h2o);
        for (int phase = 0; phase < 2; ++phase) {
          const auto p = properties(data, c, k, phase);
          for (int m = 0; m < 4; ++m) {
            if (p.number[m] == 0 || p.kappa[m] == 0)
              continue;
            const double critical =
                Kokkos::sqrt(4 * a * a * a /
                             (27 * p.kappa[m] * Kokkos::pow(p.diameter[m], 3)));
            for (std::size_t s = 0; s < count; ++s)
              if (axis(s) > 0)
                result(c, k, s) +=
                    p.number[m] * 0.5 *
                    Kokkos::erfc(2 * Kokkos::log(critical / axis(s)) /
                                 (3 * Kokkos::sqrt(2.0) *
                                  Kokkos::log(mam4::modes(m).mean_std_dev)));
          }
        }
      });
  copyBack(result, output);
}
} // namespace exaero
