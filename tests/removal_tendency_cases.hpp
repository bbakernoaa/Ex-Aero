#pragma once
/// @file removal_tendency_cases.hpp
/// @brief Shared analytic sink, mapping, and transactional-failure cases.
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exaero/BulkRemoval.hpp>
#include <exaero/BulkRemoval_C.h>
#include <exaero/ModalRemoval.hpp>
#include <exaero/ModalRemoval_C.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace removal_test {
/// @brief Relative comparison that also resolves tiny concentration rates.
inline bool close(double a, double b) {
  return std::abs(a - b) <=
         1e-12 * std::max({1e-300, std::abs(a), std::abs(b)});
}

/// @brief Exercise two populations at distinct cells/levels through C++/C.
inline void tendencies(exaero::IAerosolPackage &package, bool modal) {
  using namespace exaero;
  constexpr std::size_t nc = 2, nl = 2, nm = 2, points = nc * nl;
  const std::size_t nf = modal ? 6 : 2;
  std::vector<double> data(points * nf), mass_rate(points * nm),
      number_rate(points * nm), velocity(points * nm), mass_output(points * nf),
      number_output(points * nm), initial_number(points * nm);
  std::array<double, points> thickness{};
  View3D<double> state(data.data(), nc, nl, nf);
  for (std::size_t m = 0; m < nm; ++m)
    for (std::size_t k = 0; k < nl; ++k)
      for (std::size_t c = 0; c < nc; ++c) {
        const auto point = c + nc * k, i = point + points * m;
        thickness[point] = 10 * (point + 1);
        mass_rate[i] = (m == 0 ? 5e-4 : 2e-4) * (point + 1);
        number_rate[i] = (m == 0 ? 1e-4 : 3e-4) * (point + 1);
        velocity[i] = .01 * (m + 1) * (point + 1);
        if (modal) {
          state(c, k, m * 3) = (2 + point) * 1e-9 * (m + 1);
          state(c, k, m * 3 + 1) = (3 + point) * 1e-9 * (m + 1);
          state(c, k, m * 3 + 2) = (5 + point) * 1e7 * (m + 1);
          initial_number[i] = state(c, k, m * 3 + 2);
        } else {
          state(c, k, m) = (2 + point + m) * 1e-9;
        }
      }
  View3D<const double> input(data.data(), nc, nl, nf),
      mr(mass_rate.data(), nc, nl, nm), nr(number_rate.data(), nc, nl, nm),
      vv(velocity.data(), nc, nl, nm);
  View2D<const double> dz(thickness.data(), nc, nl);
  View3D<double> mo(mass_output.data(), nc, nl, nf),
      no(number_output.data(), nc, nl, nm);
  if (!modal) {
    View3D<double> n(initial_number.data(), nc, nl, nm);
    bulkFixedNumber(package, input, n);
  }
  const auto original = data;
  const auto original_mr = mass_rate, original_nr = number_rate;
  const auto original_v = velocity;
  const auto original_dz = thickness;

  const auto cpp = [&](bool use_velocity) {
    if (modal) {
      if (use_velocity)
        modalRemovalVelocityTendencies(package, dz, input, vv, mo, no);
      else
        modalRemovalRateTendencies(package, input, mr, nr, mo, no);
    } else if (use_velocity) {
      bulkRemovalVelocityTendencies(package, dz, input, vv, mo, no);
    } else {
      bulkRemovalRateTendencies(package, input, mr, mo, no);
    }
  };
  const auto verify = [&](bool use_velocity) {
    for (std::size_t m = 0; m < nm; ++m)
      for (std::size_t k = 0; k < nl; ++k)
        for (std::size_t c = 0; c < nc; ++c) {
          const auto point = c + nc * k, i = point + points * m;
          const double lm =
              use_velocity ? velocity[i] / thickness[point] : mass_rate[i];
          const double ln = use_velocity || !modal ? lm : number_rate[i];
          assert(close(no(c, k, m), -ln * initial_number[i]));
          if (modal) {
            assert(close(mo(c, k, m * 3), -lm * state(c, k, m * 3)));
            assert(close(mo(c, k, m * 3 + 1), -lm * state(c, k, m * 3 + 1)));
            assert(mo(c, k, m * 3 + 2) == 0);
          } else {
            assert(close(mo(c, k, m), -lm * state(c, k, m)));
          }
        }
    assert(data == original);
  };
  char error[256]{};
  int flag = -1;
  const auto abi = [&](bool use_velocity) {
    if (modal) {
      if (use_velocity)
        exaero_modal_removal_velocity_tendencies(
            &package, nc, nl, nf, nm, data.data(), velocity.data(),
            thickness.data(), mass_output.data(), number_output.data(), error,
            &flag);
      else
        exaero_modal_removal_rate_tendencies(
            &package, nc, nl, nf, nm, data.data(), mass_rate.data(),
            number_rate.data(), mass_output.data(), number_output.data(), error,
            &flag);
    } else if (use_velocity) {
      exaero_bulk_removal_velocity_tendencies(
          &package, nc, nl, nf, data.data(), velocity.data(), thickness.data(),
          mass_output.data(), number_output.data(), error, &flag);
    } else {
      exaero_bulk_removal_rate_tendencies(&package, nc, nl, nf, data.data(),
                                          mass_rate.data(), mass_output.data(),
                                          number_output.data(), error, &flag);
    }
  };
  for (bool use_velocity : {false, true}) {
    cpp(use_velocity);
    verify(use_velocity);
    const auto expected_m = mass_output, expected_n = number_output;
    std::fill(mass_output.begin(), mass_output.end(), 77);
    std::fill(number_output.begin(), number_output.end(), 88);
    abi(use_velocity);
    assert(flag == 0 && error[0] == '\0');
    assert(mass_output == expected_m && number_output == expected_n);
  }
  mass_rate[0] = 0;
  cpp(false);
  verify(false);
  assert(mass_output[0] == 0);
  assert(mass_output[points * (modal ? 3 : 1)] < 0);
  mass_rate = original_mr;
  std::fill(velocity.begin(), velocity.end(), 0);
  cpp(true);
  for (double x : mass_output)
    assert(x == 0);
  for (double x : number_output)
    assert(x == 0);
  velocity = original_v;

  const auto reset = [&] {
    std::fill(mass_output.begin(), mass_output.end(), 77);
    std::fill(number_output.begin(), number_output.end(), 88);
  };
  const auto unchanged = [&] {
    for (double x : mass_output)
      assert(x == 77);
    for (double x : number_output)
      assert(x == 88);
  };
  const auto failure = [&](bool use_velocity) {
    reset();
    bool rejected = false;
    try {
      cpp(use_velocity);
    } catch (const std::exception &e) {
      rejected = true;
      assert(std::string(e.what()).starts_with("FATAL ERROR:"));
    }
    assert(rejected);
    unchanged();
  };
  for (double bad : {-1., std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity()}) {
    mass_rate.back() = bad;
    failure(false);
    abi(false);
    assert(flag == 1 && error[0] != '\0');
    unchanged();
    mass_rate = original_mr;
    velocity.back() = bad;
    failure(true);
    velocity = original_v;
    if (modal) {
      number_rate.back() = bad;
      failure(false);
      number_rate = original_nr;
    }
  }
  for (double bad : {0., -1., std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity()}) {
    thickness.back() = bad;
    failure(true);
    abi(true);
    assert(flag == 1);
    unchanged();
    thickness = original_dz;
  }
  mr = View3D<const double>(mass_rate.data(), nc, nl, 1);
  failure(false);
  mr = View3D<const double>(mass_rate.data(), nc, nl, nm);
  for (double bad : {-1., std::numeric_limits<double>::quiet_NaN(),
                     std::numeric_limits<double>::infinity()}) {
    data.back() = bad;
    failure(false);
    failure(true);
    data = original;
  }
  if (modal) {
    data.back() = 0;
    failure(false);
    data = original;
    nr = View3D<const double>(number_rate.data(), nc, nl, 1);
    failure(false);
    nr = View3D<const double>(number_rate.data(), nc, nl, nm);
  }
  mass_rate.back() = std::numeric_limits<double>::max();
  data[points * (modal ? 3 : 1) + points - 1] = 2;
  failure(false);
  mass_rate = original_mr;
  data = original;
  mo = View3D<double>(data.data(), nc, nl, nf);
  bool rejected = false;
  try {
    cpp(false);
  } catch (const std::exception &) {
    rejected = true;
  }
  assert(rejected && data == original);
  mo = View3D<double>(mass_output.data(), nc, nl, nf);
  no = View3D<double>(mass_output.data(), nc, nl, nm);
  reset();
  rejected = false;
  try {
    cpp(true);
  } catch (const std::exception &) {
    rejected = true;
  }
  assert(rejected);
  unchanged();
  no = View3D<double>(number_output.data(), nc, nl, nm);
  no = View3D<double>(mass_rate.data(), nc, nl, nm);
  reset();
  rejected = false;
  try {
    cpp(false);
  } catch (const std::exception &) {
    rejected = true;
  }
  assert(rejected && mass_rate == original_mr);
  unchanged();
  no = View3D<double>(number_output.data(), nc, nl, nm);
  reset();
  if (modal) {
    exaero_modal_removal_rate_tendencies(&package, nc, nl, nf - 1, nm,
                                         data.data(), mass_rate.data(),
                                         number_rate.data(), mass_output.data(),
                                         number_output.data(), error, &flag);
  } else {
    exaero_bulk_removal_rate_tendencies(&package, nc, nl, nf, data.data(),
                                        nullptr, mass_output.data(),
                                        number_output.data(), error, &flag);
  }
  assert(flag == 1);
  unchanged();
  if (modal)
    exaero_modal_removal_velocity_tendencies(
        &package, nc, nl, nf, nm, data.data(), nullptr, thickness.data(),
        mass_output.data(), number_output.data(), error, &flag);
  else
    exaero_bulk_removal_velocity_tendencies(
        nullptr, nc, nl, nf, data.data(), velocity.data(), thickness.data(),
        mass_output.data(), number_output.data(), error, &flag);
  assert(flag == 1);
  unchanged();

  // A zero-mass/zero-number population remains empty; no number is inferred.
  std::fill(data.begin(), data.end(), 0);
  cpp(false);
  for (double x : mass_output)
    assert(x == 0);
  for (double x : number_output)
    assert(x == 0);
  data = original;
  cpp(false);
  const auto ambient_mass = mass_output, ambient_number = number_output;
  for (double &x : data)
    x *= .5;
  cpp(false);
  for (std::size_t i = 0; i < data.size(); ++i)
    assert(data[i] == .5 * original[i] &&
           close(mass_output[i], .5 * ambient_mass[i]));
  for (std::size_t i = 0; i < number_output.size(); ++i)
    assert(close(number_output[i], .5 * ambient_number[i]));
  data = original;
  assert(mass_rate == original_mr && number_rate == original_nr &&
         velocity == original_v && thickness == original_dz);
}
} // namespace removal_test
