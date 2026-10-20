/// @file test_package_selection.cpp
/// @brief Factory availability and error-path checks in every build
/// configuration.
#include <exaero/IAerosolPackage_C.h>
#include <exaero/PackageFactory.hpp>
#include <stdexcept>

int main() {
  const auto bulk = exaero::createPackage("gocart");
  if (bulk->representation() != exaero::AerosolRepresentation::Bulk ||
      !exaero::packageAvailable("gocart") ||
      exaero::packageAvailable("unknown"))
    return 1;
  char error[256]{};
  int status = 0;
  if (exaero_create_package("unknown", error, &status) || status != 1 ||
      !error[0])
    return 2;
  if (exaero_create_package(nullptr, error, &status) || status != 1)
    return 3;
  auto *legacy = exaero_create_gocart_package();
  if (!legacy || exaero_get_representation(legacy) != 0)
    return 4;
  exaero_free_package(legacy);
  auto *modal = exaero_create_package("mam4xx", error, &status);
  if (exaero::packageAvailable("mam4xx")) {
    if (!modal || status || exaero_get_representation(modal) != 1 ||
        exaero_get_state_field_count(modal) != 64)
      return 5;
  } else if (modal || !status || !error[0]) {
    return 6;
  }
  exaero_free_package(modal);
  return 0;
}
