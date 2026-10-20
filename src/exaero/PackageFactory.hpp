#pragma once
/// @file PackageFactory.hpp
/// @brief Dependency-free factory for bulk and modal aerosol packages.
#include <exaero/IAerosolPackage.hpp>
#include <memory>
#include <string_view>

namespace exaero {
/// @brief Construct "gocart" (bulk) or "mam4xx" (modal), without
/// initialization.
/// @throws std::invalid_argument for unknown or disabled backends.
std::unique_ptr<IAerosolPackage> createPackage(std::string_view backend);
/// @brief Whether a backend was compiled into this library.
bool packageAvailable(std::string_view backend) noexcept;
} // namespace exaero
