/// @file api.cpp
/// @brief Legacy extern "C" implementation over the chemistry registry.
///
/// @warning Not part of the active GOCART build (see src/CMakeLists.txt);
/// retained with the api.hpp/sorter.hpp/state.hpp scaffold. The supported
/// host entry point is the IAerosolPackage_C.h CCPP surface. Implements
/// the init/solve/finalize lifecycle: a single global context, an SZA
/// sort, then a registry dispatch of the chemistry mechanism.
///
/// @note Error handling: every entry point is exception-safe (no throw
/// crosses @c extern "C"); failures print a @c "FATAL ERROR:" line and
/// return a nonzero code per the EE2 convention.
#include "exaero/api.hpp"
#include "exaero/registry.hpp"
#include "exaero/sorter.hpp"
#include "exaero/state.hpp"
#include <iostream>
#include <memory>

/// @brief Process-wide run context; created by exaero_init, torn down by
/// exaero_finalize. RAII-owned so a solve-time throw cannot leak device
/// scratch.
static std::unique_ptr<exaero::ExaeroContext> global_ctx = nullptr;

extern "C" {

// @copydoc exaero_init
int exaero_init(const char *mechanism_name, int32_t total_cells) {
  try {
    if (!exaero::Registry::has_mechanism(mechanism_name)) {
      std::cerr << "FATAL ERROR: Mechanism " << mechanism_name
                << " not found in Exaero registry." << std::endl;
      return 1;
    }
    global_ctx =
        std::make_unique<exaero::ExaeroContext>(mechanism_name, total_cells);
    return 0;
  } catch (...) {
    return 1;
  }
}

// @copydoc exaero_solve
int exaero_solve(double *conc_ptr, double *met_ptr, double dt) {
  if (!global_ctx)
    return 1;
  try {
    // Wrap the caller arrays zero-copy, sort by SZA, then dispatch the
    // active mechanism over the sorted cell order.
    exaero::UnmanagedDeviceState state(conc_ptr, met_ptr, 1,
                                       global_ctx->total_grid_cells, 1, 1);
    exaero::SZA_Sorter::sort_workload(*global_ctx, state);
    exaero::Registry::dispatch(*global_ctx, state, dt);
    return 0;
  } catch (...) {
    return 1;
  }
}

// @copydoc exaero_finalize
void exaero_finalize() { global_ctx.reset(); }

} // extern "C"
