/// @file test_modal_main.cpp
/// @brief Kokkos-only test lifecycle, isolated from public mdspan headers.
/// @author EX-aero contributors
/// @details History: 2026-10. Usage: CTest. No file input/output.
#define CATCH_CONFIG_RUNNER
#include <Kokkos_Core.hpp>
#include <catch2/catch.hpp>
int main(int argc, char **argv) {
  Kokkos::ScopeGuard guard(argc, argv);
  return Catch::Session().run(argc, argv);
}
