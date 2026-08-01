#ifndef EXAERO_STATE_HPP
#define EXAERO_STATE_HPP

#include <Kokkos_Core.hpp>
#include <string>

namespace exaero {

struct UnmanagedDeviceState {
    using Layout = Kokkos::LayoutLeft;
    using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
    using View4D = Kokkos::View<double****, Layout, MemorySpace, Kokkos::MemoryUnmanaged>;

    View4D concentrations;
    View4D meteorology;

    UnmanagedDeviceState(double* conc_ptr, double* met_ptr, int n_species, int n_lon, int n_lat, int n_alt)
        : concentrations(conc_ptr, n_species, n_lon, n_lat, n_alt),
          meteorology(met_ptr, 10, n_lon, n_lat, n_alt)
    {}
};

struct SZA_Pair {
    double sza;
    int index;
    KOKKOS_INLINE_FUNCTION bool operator<(const SZA_Pair& other) const {
        return sza < other.sza;
    }
};

struct ExaeroContext {
    std::string active_mechanism;
    int32_t total_grid_cells;
    Kokkos::View<int*> sza_sorted_indices;
    Kokkos::View<SZA_Pair*> sorter_scratch;

    ExaeroContext(const std::string& mech_name, int32_t cells)
        : active_mechanism(mech_name), 
          total_grid_cells(cells),
          sza_sorted_indices("sza_sorted_indices", cells),
          sorter_scratch("sorter_scratch", cells)
    {}
};

} // namespace exaero
#endif // EXAERO_STATE_HPP
