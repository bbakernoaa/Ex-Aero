#ifndef EXAERO_API_HPP
#define EXAERO_API_HPP
/// @file api.hpp
/// @brief Legacy extern "C" init/solve/finalize lifecycle for the C API.
///
/// The original flat-pointer entry points (pre-IAerosolPackage design):
/// a single global context initialized by name, then one solve call per
/// step over caller-allocated concentration and meteorology arrays.
/// @note The GOCART/GEOSmie production surface is IAerosolPackage_C.h;
/// this header is retained for the chemistry-registry path.
#include <stdint.h>

extern "C" {

/// @brief Create the global run context for a registered mechanism.
/// @param mechanism_name Chemistry mechanism name (must be in the registry).
/// @param total_cells Horizontal grid-cell count (sizes device scratch).
/// @return 0 on success, 1 on unknown mechanism or allocation failure.
int exaero_init(const char *mechanism_name, int32_t total_cells);

/// @brief Advance one solve step over caller-allocated arrays.
/// @param conc_ptr [kg/m^3] concentrations (species, lon, lat, alt).
/// @param met_ptr Meteorology array (10 vars, lon, lat, alt); slot 0 = SZA.
/// @param dt Integration timestep [s].
/// @return 0 on success, 1 when uninitialized or on any exception.
int exaero_solve(double *conc_ptr, double *met_ptr, double dt);

/// @brief Release the global run context (safe to call twice).
void exaero_finalize();
}

#endif // EXAERO_API_HPP
