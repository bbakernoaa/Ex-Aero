#pragma once
/// @file CloudIceConfig.hpp
/// @brief Radiation-independent physical roles and ice size rules.
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace exaero {
enum class IceCategory : int {
  Dust = 0,
  Sulfate = 1,
  BlackCarbon = 2,
  StratosphericSulfate = 3
};
enum class IceSizeRuleKind : int {
  None = 0,
  All = 1,
  FixedFraction = 2,
  AboveDryDiameter = 3
};
struct IceSizeRule {
  IceSizeRuleKind kind = IceSizeRuleKind::None;
  double value = 0;
};
struct IcePopulationRule {
  IceCategory category = IceCategory::Dust;
  IceSizeRule normal;
  std::optional<IceSizeRule> preexisting_ice;
  bool heterogeneous_freezing = false;
};

inline bool validIceCategory(IceCategory category) {
  return category >= IceCategory::Dust &&
         category <= IceCategory::StratosphericSulfate;
}
inline void validateIceSizeRule(const IceSizeRule &rule) {
  if (!std::isfinite(rule.value))
    throw std::invalid_argument("FATAL ERROR: nonfinite ice size rule");
  switch (rule.kind) {
  case IceSizeRuleKind::None:
  case IceSizeRuleKind::All:
    if (rule.value == 0)
      return;
    break;
  case IceSizeRuleKind::FixedFraction:
    if (rule.value >= 0 && rule.value <= 1)
      return;
    break;
  case IceSizeRuleKind::AboveDryDiameter:
    if (rule.value > 0)
      return;
    break;
  }
  throw std::invalid_argument("FATAL ERROR: invalid ice size rule");
}
inline void validateIcePopulationRule(const IcePopulationRule &rule) {
  if (!validIceCategory(rule.category))
    throw std::invalid_argument("FATAL ERROR: invalid ice category");
  validateIceSizeRule(rule.normal);
  if (rule.preexisting_ice)
    validateIceSizeRule(*rule.preexisting_ice);
}
inline void
validateIceConfiguration(const std::vector<IceCategory> &categories,
                         const std::vector<IcePopulationRule> &rules) {
  unsigned category_mask = 0;
  unsigned rule_mask = 0;
  for (const auto category : categories) {
    if (!validIceCategory(category))
      throw std::invalid_argument("FATAL ERROR: invalid ice category");
    const unsigned bit = 1u << static_cast<int>(category);
    if (category_mask & bit)
      throw std::invalid_argument("FATAL ERROR: duplicate ice category");
    category_mask |= bit;
  }
  for (const auto &rule : rules) {
    validateIcePopulationRule(rule);
    const unsigned bit = 1u << static_cast<int>(rule.category);
    if (!(category_mask & bit))
      throw std::invalid_argument("FATAL ERROR: ice rule without category");
    if (rule_mask & bit)
      throw std::invalid_argument("FATAL ERROR: duplicate ice rule");
    rule_mask |= bit;
  }
}
inline IceCategory iceCategoryFromName(std::string_view name) {
  if (name == "dust")
    return IceCategory::Dust;
  if (name == "sulfate")
    return IceCategory::Sulfate;
  if (name == "black_carbon")
    return IceCategory::BlackCarbon;
  if (name == "stratospheric_sulfate")
    return IceCategory::StratosphericSulfate;
  throw std::invalid_argument("FATAL ERROR: unknown ice category");
}
inline IceSizeRuleKind iceSizeRuleKindFromName(std::string_view name) {
  if (name == "none")
    return IceSizeRuleKind::None;
  if (name == "all")
    return IceSizeRuleKind::All;
  if (name == "fixed_fraction")
    return IceSizeRuleKind::FixedFraction;
  if (name == "above_dry_diameter")
    return IceSizeRuleKind::AboveDryDiameter;
  throw std::invalid_argument("FATAL ERROR: unknown ice size rule kind");
}
} // namespace exaero
