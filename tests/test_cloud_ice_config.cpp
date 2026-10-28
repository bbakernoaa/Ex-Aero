#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <exaero/CloudIceConfig.hpp>
#include <exaero/CloudIceMath.hpp>
#include <exaero/GocartConfig.hpp>
#include <exaero/ModalConfig.hpp>

int main() {
  using namespace exaero;
  const IceSizeRule all{IceSizeRuleKind::All, 0};
  const IceSizeRule none{IceSizeRuleKind::None, 0};
  const IceSizeRule fixed{IceSizeRuleKind::FixedFraction, 0.25};
  assert(iceSizeWeight(all, 1e-7, 1.5) == 1);
  assert(iceSizeWeight(none, 1e-7, 1.5) == 0);
  assert(iceSizeWeight(fixed, 1e-7, 1.5) == 0.25);

  const IceSizeRule cutoff{IceSizeRuleKind::AboveDryDiameter, 1e-7};
  assert(iceSizeWeight(cutoff, 1e-7, 1) == 0);
  assert(iceSizeWeight(cutoff, 1.01e-7, 1) == 1);
  assert(std::abs(iceSizeWeight(cutoff, 1e-7, std::exp(1.0)) - 0.5) <
         1e-14);
  assert(std::abs(iceSizeWeight({IceSizeRuleKind::AboveDryDiameter,
                                std::exp(1.0) * 1e-7},
                               1e-7, std::exp(1.0)) -
                  0.15865525393145707) <
         1e-14);

  bool rejected = false;
  try {
    (void)iceSizeWeight({IceSizeRuleKind::FixedFraction, 1.1}, 1e-7, 1.5);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected);

  GocartSpeciesConfig bulk;
  bulk.ice_categories = {IceCategory::Dust};
  bulk.ice_rules = {{IceCategory::Dust, all, std::nullopt, true}};
  validateIceConfiguration(bulk.ice_categories, bulk.ice_rules);

  ModalSpeciesConfig material;
  material.ice_categories = {IceCategory::Dust};
  ModalModeConfig mode;
  mode.ice_rules = bulk.ice_rules;
  validateIceConfiguration(material.ice_categories, mode.ice_rules);

  rejected = false;
  try {
    validateIceConfiguration({IceCategory::Dust, IceCategory::Dust},
                             bulk.ice_rules);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected);
  rejected = false;
  try {
    validateIceConfiguration({IceCategory::Sulfate}, bulk.ice_rules);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected);
}
