#ifndef EXAERO_SORTER_HPP
#define EXAERO_SORTER_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include "exaero/state.hpp"

namespace exaero {

class SZA_Sorter {
public:
    static void sort_workload(ExaeroContext& ctx, UnmanagedDeviceState& state) {
        auto scratch = ctx.sorter_scratch;
        auto indices = ctx.sza_sorted_indices;
        auto meteo = state.meteorology;
        
        Kokkos::parallel_for("PopulateScratch", ctx.total_grid_cells, KOKKOS_LAMBDA(const int i) {
            scratch(i).sza = meteo(0, i, 0, 0);
            scratch(i).index = i;
        });

        // Use Kokkos::sort on the pair struct
        Kokkos::sort(scratch);

        Kokkos::parallel_for("ExtractIndices", ctx.total_grid_cells, KOKKOS_LAMBDA(const int i) {
            indices(i) = scratch(i).index;
        });
    }
};

} // namespace exaero
#endif // EXAERO_SORTER_HPP
