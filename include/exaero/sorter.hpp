#ifndef EXAERO_SORTER_HPP
#define EXAERO_SORTER_HPP
/// @file sorter.hpp
/// @brief Solar-zenith-angle workload sorter for the legacy C API.
///
/// Computes a device-side permutation of grid cells ordered by
/// increasing SZA (meteo slot 0). Photolysis and thermodynamics kernels
/// consume the permutation so sunlit cells are processed together,
/// improving cache locality and load balance across thousands of ranks.
#include "exaero/state.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>

namespace exaero {

/// @brief Builds the SZA-sorted cell permutation in the run context.
class SZA_Sorter {
public:
  /// @brief Populate, sort, and extract the SZA permutation for one step.
  ///
  /// Three device passes: (1) fill scratch with (SZA, index) pairs from
  /// meteo slot 0, (2) Kokkos::sort on the pair struct (ascending SZA),
  /// (3) write the permuted indices back to ctx.sza_sorted_indices.
  /// @param ctx Per-run context (owns the scratch and output views).
  /// @param state Unmanaged device state (meteorology read only).
  static void sort_workload(ExaeroContext &ctx, UnmanagedDeviceState &state) {
    auto scratch = ctx.sorter_scratch;
    auto indices = ctx.sza_sorted_indices;
    auto meteo = state.meteorology;

    Kokkos::parallel_for(
        "PopulateScratch", ctx.total_grid_cells, KOKKOS_LAMBDA(const int i) {
          scratch(i).sza = meteo(0, i, 0, 0);
          scratch(i).index = i;
        });

    // Use Kokkos::sort on the pair struct
    Kokkos::sort(scratch);

    Kokkos::parallel_for(
        "ExtractIndices", ctx.total_grid_cells,
        KOKKOS_LAMBDA(const int i) { indices(i) = scratch(i).index; });
  }
};

} // namespace exaero
#endif // EXAERO_SORTER_HPP
