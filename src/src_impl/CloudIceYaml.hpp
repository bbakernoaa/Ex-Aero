#pragma once
/// @file CloudIceYaml.hpp
/// @brief Shared YAML conversion for bulk and modal physical ice rules.
#include <exaero/CloudIceConfig.hpp>
#include <stdexcept>
#include <string>
#include <yaml-cpp/yaml.h>

namespace exaero {
namespace cloud_ice_yaml {
inline void keys(const YAML::Node &node,
                 std::initializer_list<const char *> allowed) {
  if (!node.IsMap())
    throw std::invalid_argument("FATAL ERROR: ice rule must be a mapping");
  for (const auto &entry : node) {
    const auto key = entry.first.as<std::string>();
    bool found = false;
    for (const auto *name : allowed)
      found = found || key == name;
    if (!found)
      throw std::invalid_argument("FATAL ERROR: unknown ice rule key: " + key);
  }
}
inline IceSizeRule sizeRule(const YAML::Node &node) {
  keys(node, {"kind", "value"});
  if (!node["kind"])
    throw std::invalid_argument("FATAL ERROR: ice size rule needs kind");
  IceSizeRule rule{iceSizeRuleKindFromName(node["kind"].as<std::string>()),
                   node["value"] ? node["value"].as<double>() : 0};
  validateIceSizeRule(rule);
  return rule;
}
inline std::vector<IceCategory> categories(const YAML::Node &node) {
  if (!node)
    return {};
  if (!node.IsSequence())
    throw std::invalid_argument("FATAL ERROR: ice_categories must be a sequence");
  std::vector<IceCategory> result;
  for (const auto &value : node)
    result.push_back(iceCategoryFromName(value.as<std::string>()));
  validateIceConfiguration(result, {});
  return result;
}
inline std::vector<IcePopulationRule> rules(const YAML::Node &node) {
  if (!node)
    return {};
  if (!node.IsSequence())
    throw std::invalid_argument("FATAL ERROR: ice_rules must be a sequence");
  std::vector<IcePopulationRule> result;
  for (const auto &value : node) {
    keys(value, {"category", "kind", "value", "preexisting_ice",
                 "heterogeneous_freezing"});
    if (!value["category"])
      throw std::invalid_argument("FATAL ERROR: ice rule needs category");
    IcePopulationRule rule;
    rule.category =
        iceCategoryFromName(value["category"].as<std::string>());
    if (!value["kind"])
      throw std::invalid_argument("FATAL ERROR: ice size rule needs kind");
    rule.normal = {iceSizeRuleKindFromName(value["kind"].as<std::string>()),
                   value["value"] ? value["value"].as<double>() : 0};
    validateIceSizeRule(rule.normal);
    if (value["preexisting_ice"])
      rule.preexisting_ice = sizeRule(value["preexisting_ice"]);
    if (value["heterogeneous_freezing"])
      rule.heterogeneous_freezing =
          value["heterogeneous_freezing"].as<bool>();
    result.push_back(rule);
  }
  return result;
}
} // namespace cloud_ice_yaml
} // namespace exaero
