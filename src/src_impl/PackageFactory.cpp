/// @file PackageFactory.cpp
/// @brief Runtime backend selection. No Kokkos headers enter this translation
/// unit.
#include <exaero/PackageFactory.hpp>
#include <gocart/GocartPackage.hpp>
#include <modal/ModalPackage.hpp>
#ifdef EXAERO_ENABLE_MAM4XX
#include <mam4xx/Mam4xxPackage.hpp>
#endif

namespace exaero {
bool packageAvailable(std::string_view backend) noexcept {
  if (backend == "gocart")
    return true;
  if (backend == "modal")
    return true;
#ifdef EXAERO_ENABLE_MAM4XX
  if (backend == "mam4xx")
    return true;
#endif
  return false;
}

std::unique_ptr<IAerosolPackage> createPackage(std::string_view backend) {
  if (backend == "gocart")
    return std::make_unique<GocartPackage>();
  if (backend == "modal")
    return std::make_unique<ModalPackage>();
#ifdef EXAERO_ENABLE_MAM4XX
  if (backend == "mam4xx")
    return std::make_unique<Mam4xxPackage>();
#endif
  throw std::invalid_argument(
      "FATAL ERROR: unknown or disabled aerosol backend: " +
      std::string(backend));
}
} // namespace exaero
