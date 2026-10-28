#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <stdexcept>
#include <yaml-cpp/yaml.h>
#include <CloudIceYaml.hpp>

int main() {
  const auto config = YAML::Load(R"(
    ice_categories: [dust]
    ice_rules:
      - category: dust
        kind: above_dry_diameter
        value: 1.0e-7
        preexisting_ice: {kind: all}
        heterogeneous_freezing: true
  )");
  const auto categories =
      exaero::cloud_ice_yaml::categories(config["ice_categories"]);
  const auto rules = exaero::cloud_ice_yaml::rules(config["ice_rules"]);
  assert(categories.size() == 1 &&
         categories[0] == exaero::IceCategory::Dust);
  assert(rules.size() == 1 &&
         rules[0].normal.kind == exaero::IceSizeRuleKind::AboveDryDiameter &&
         rules[0].normal.value == 1e-7 && rules[0].preexisting_ice &&
         rules[0].heterogeneous_freezing);
  exaero::validateIceConfiguration(categories, rules);

  bool rejected = false;
  try {
    const auto invalid = YAML::Load("ice_categories: [unknown]");
    (void)exaero::cloud_ice_yaml::categories(invalid["ice_categories"]);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected);
}
