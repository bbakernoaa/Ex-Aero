// Standalone public-header compile check (contract
//).
//
// PURPOSE: prove the EX-aero public attribute surface is consumable WITHOUT
// Kokkos on the include path. This translation unit is compiled by the
// `exaero_standalone_public_check` target, which deliberately links nothing and
// sets no Kokkos include directories: if any public header pulled in
// <Kokkos_Core.hpp> (directly or transitively) this file would fail to
// configure/compile. That is the gate.
//
// It also exercises the header contents at compile time (types complete,
// constants usable, the POD provenance record is trivially constructible) so
// the check is not vacuous.

#include <exaero/AerosolIndices.hpp> // diagnostic / optical index constants
#include <exaero/AttributeQuery.hpp> // AttributeCategory, AttributeStatus, ProvenanceInfo,
// SpeciesCurveConfig, *_indices namespaces (new surface)
#include <exaero/Environment.hpp> // initialize/finalize_environment declarations
#include <exaero/IAerosolPackage_C.h> // C entry points incl. exaero_get_moment_counts

#include <type_traits>

namespace {

// The status / delivery enums must be int-typed and stable across the C
// boundary.
static_assert(std::is_enum_v<exaero::AttributeStatus>,
              "AttributeStatus is an enum");
static_assert(
    std::is_same_v<std::underlying_type_t<exaero::AttributeStatus>, int>,
    "AttributeStatus underlying type is int (C boundary)");
static_assert(static_cast<int>(exaero::AttributeStatus::AvailableConfig) == 5,
              "AttributeStatus ordinal stable");
static_assert(static_cast<int>(exaero::DeliverySource::Config) == 2,
              "DeliverySource ordinal stable");

// ProvenanceInfo must be a POD that crosses the FFI by value (no Kokkos, no
// std::string).
static_assert(std::is_trivially_copyable_v<exaero::ProvenanceInfo>,
              "ProvenanceInfo must be trivially copyable");

// The polarized-moment element stride is the documented P11,P12,P33,P34,P22,P44
// ordering.
static_assert(exaero::polarized_moment_indices::ELEMENT_STRIDE == 6,
              "element ordering documented (data-model)");

// Curve-mapping config carries dynamic axes: vectors + counts, never fixed-size
// arrays.
static_assert(std::is_same_v<decltype(exaero::SpeciesCurveConfig::radius_nodes),
                             std::vector<double>>,
              "radius axis is dynamic");

// Compile-time use of the category-bit helper (activation mask builder,
//).
constexpr int kMask =
    exaero::attribute_category_bit(exaero::AttributeCategory::PolarizedMoment);
static_assert(kMask == (1 << 2),
              "attribute_category_bit builds (1 << category)");

} // namespace

int main() {
  // Instantiate the public descriptors so their definitions are required (not
  // just declared).
  exaero::ProvenanceInfo prov{};
  prov.num_pol = 6;
  prov.num_moment = 3;
  exaero::SpeciesCurveConfig cfg;
  cfg.species_name = "dust_mode_1";
  cfg.source_label = "DU";
  cfg.radius_nodes = {1.0e-07, 8.0e-07, 3.0e-06};

  // Reference a C entry point symbol (declared, not called: no library linked
  // here). Taking the address proves the extern "C" declaration is visible and
  // well-formed.
  (void)&exaero_get_moment_counts;
  (void)&exaero_query_attribute;

  return (prov.num_pol == 6 && prov.num_moment == 3 &&
          cfg.radius_nodes.size() == 3)
             ? 0
             : 1;
}
