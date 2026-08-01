#include "exaero/api.hpp"
#include "exaero/state.hpp"
#include "exaero/sorter.hpp"
#include "exaero/registry.hpp"
#include <iostream>
#include <memory>

static std::unique_ptr<exaero::ExaeroContext> global_ctx = nullptr;

extern "C" {

int exaero_init(const char* mechanism_name, int32_t total_cells) {
    try {
        if (!exaero::Registry::has_mechanism(mechanism_name)) {
            std::cerr << "FATAL ERROR: Mechanism " << mechanism_name << " not found in Exaero registry." << std::endl;
            return 1;
        }
        global_ctx = std::make_unique<exaero::ExaeroContext>(mechanism_name, total_cells);
        return 0;
    } catch (...) {
        return 1;
    }
}

int exaero_solve(double* conc_ptr, double* met_ptr, double dt) {
    if (!global_ctx) return 1;
    try {
        exaero::UnmanagedDeviceState state(conc_ptr, met_ptr, 1, global_ctx->total_grid_cells, 1, 1);
        exaero::SZA_Sorter::sort_workload(*global_ctx, state);
        exaero::Registry::dispatch(*global_ctx, state, dt);
        return 0;
    } catch (...) {
        return 1;
    }
}

void exaero_finalize() {
    global_ctx.reset();
}

} // extern "C"
