#pragma once
/// @file CloudIceMath.hpp
/// @brief Pure number-tail weights for host-selected ice size rules.
#include <algorithm>
#include <cmath>
#include <exaero/CloudIceConfig.hpp>
#include <stdexcept>

namespace exaero {
inline double iceSizeWeight(const IceSizeRule &rule,
                            double dry_number_median_diameter,
                            double geometric_std_dev) {
  validateIceSizeRule(rule);
  if (!(std::isfinite(dry_number_median_diameter) &&
        dry_number_median_diameter > 0 && std::isfinite(geometric_std_dev) &&
        geometric_std_dev >= 1))
    throw std::invalid_argument("FATAL ERROR: invalid ice size distribution");
  switch (rule.kind) {
  case IceSizeRuleKind::None:
    return 0;
  case IceSizeRuleKind::All:
    return 1;
  case IceSizeRuleKind::FixedFraction:
    return rule.value;
  case IceSizeRuleKind::AboveDryDiameter:
    if (geometric_std_dev == 1)
      return dry_number_median_diameter > rule.value ? 1 : 0;
    return std::clamp(
        0.5 * std::erfc(std::log(rule.value / dry_number_median_diameter) /
                        (std::sqrt(2.0) * std::log(geometric_std_dev))),
        0.0, 1.0);
  }
  throw std::invalid_argument("FATAL ERROR: invalid ice size rule kind");
}
} // namespace exaero
