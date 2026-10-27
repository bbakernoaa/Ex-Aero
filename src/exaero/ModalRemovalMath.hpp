#pragma once
/// @file ModalRemovalMath.hpp
/// @brief Dry lognormal modal moments used by explicit host removal utilities.
/// @details Pure, dependency-free mathematical operations; no state mutation.
#include <cmath>
#include <limits>
#include <stdexcept>

namespace exaero {
inline double modalNumberVolume(double radius, double sigma) {
  if (!(std::isfinite(radius) && radius > 0 && std::isfinite(sigma) &&
        sigma > 1))
    throw std::invalid_argument(
        "FATAL ERROR: invalid dry modal radius or width");
  const double l = std::log(sigma);
  const double value = (4.0 / 3.0) * std::acos(-1.0) * radius * radius *
                       radius * std::exp(4.5 * l * l);
  if (!(std::isfinite(value) && value > 0))
    throw std::overflow_error(
        "FATAL ERROR: unrepresentable modal particle volume");
  return value;
}

inline double modalRadius(double volume, double number, double sigma) {
  if (!(std::isfinite(volume) && volume > 0 && std::isfinite(number) &&
        number > 0 && std::isfinite(sigma) && sigma > 1))
    throw std::invalid_argument("FATAL ERROR: invalid modal volume or number");
  const double l = std::log(sigma);
  const double r = std::cbrt(volume / number) *
                   std::cbrt(3.0 / (4.0 * std::acos(-1.0))) *
                   std::exp(-1.5 * l * l);
  if (!(std::isfinite(r) && r > 0))
    throw std::overflow_error("FATAL ERROR: unrepresentable modal radius");
  return r;
}

inline double modalNumberAtRadius(double volume, double radius, double sigma) {
  if (!(std::isfinite(volume) && volume > 0))
    throw std::invalid_argument("FATAL ERROR: invalid modal volume");
  const double n = volume / modalNumberVolume(radius, sigma);
  if (!(std::isfinite(n) && n > 0))
    throw std::overflow_error(
        "FATAL ERROR: unrepresentable modal number bound");
  return n;
}

inline double modalMassFractionBelow(double diameter, double radius,
                                     double sigma) {
  if (!(diameter >= 0) || std::isnan(diameter))
    throw std::invalid_argument("FATAL ERROR: invalid diameter edge");
  (void)modalNumberVolume(radius, sigma);
  if (diameter == 0)
    return 0;
  if (std::isinf(diameter))
    return 1;
  const double l = std::log(sigma);
  const double z = (std::log(diameter) - std::log(2.0 * radius) - 3.0 * l * l) /
                   (std::sqrt(2.0) * l);
  return .5 * std::erfc(-z);
}

inline double modalMassFractionAbove(double diameter, double radius,
                                     double sigma) {
  if (!(diameter >= 0) || std::isnan(diameter))
    throw std::invalid_argument("FATAL ERROR: invalid diameter edge");
  (void)modalNumberVolume(radius, sigma);
  if (diameter == 0)
    return 1;
  if (std::isinf(diameter))
    return 0;
  const double l = std::log(sigma);
  const double z = (std::log(diameter) - std::log(2.0 * radius) - 3.0 * l * l) /
                   (std::sqrt(2.0) * l);
  return .5 * std::erfc(z);
}
} // namespace exaero
