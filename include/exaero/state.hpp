#ifndef EXAERO_STATE_HPP
#define EXAERO_STATE_HPP
/// @file state.hpp
/// @brief Device-resident state and context types for the legacy C API.
///
/// Backing types for the extern "C" surface in api.hpp (pre-IAerosolPackage
/// scaffold): unmanaged zero-copy views over host-allocated concentration
/// and meteorology arrays plus the per-run sorting context.
///
/// @note Layout: all views are Kokkos::LayoutLeft (column-major) so the
/// Fortran host arrays are consumed zero-copy; index order is
/// (species/meteo_var, lon, lat, alt).
#include <Kokkos_Core.hpp>
#include <string>

namespace exaero {

/// @brief Non-owning device views over caller-allocated model arrays.
///
/// Wraps raw pointers as unmanaged Kokkos views (no allocation, no
/// copy) so a Fortran/C++ host owns the memory and the device kernels
/// alias it directly.
struct UnmanagedDeviceState {
  using Layout = Kokkos::LayoutLeft; ///< Column-major (Fortran order).
  using MemorySpace =
      Kokkos::DefaultExecutionSpace::memory_space; ///< Execution memory space
                                                   ///< (device on GPU builds).
  using View4D =
      Kokkos::View<double ****, Layout, MemorySpace,
                   Kokkos::MemoryUnmanaged>; ///< 4D unmanaged view (species,
                                             ///< lon, lat, alt).

  View4D concentrations; ///< Aerosol concentrations (species, lon, lat, alt)
                         ///< [kg/m^3].
  View4D meteorology;    ///< Meteorological fields; slot 0 carries the solar
                         ///< zenith angle (SZA) used for sorting.

  /// @brief Alias caller pointers as unmanaged views (no ownership taken).
  /// @param conc_ptr [kg/m^3] concentration array (n_species, n_lon, n_lat,
  ///        n_alt).
  /// @param met_ptr Meteorology array (10 variables, n_lon, n_lat, n_alt).
  /// @param n_species Species axis length.
  /// @param n_lon Longitude cell count.
  /// @param n_lat Latitude cell count.
  /// @param n_alt Altitude/level count.
  UnmanagedDeviceState(double *conc_ptr, double *met_ptr, int n_species,
                       int n_lon, int n_lat, int n_alt)
      : concentrations(conc_ptr, n_species, n_lon, n_lat, n_alt),
        meteorology(met_ptr, 10, n_lon, n_lat, n_alt) {}
};

/// @brief (SZA, cell index) pair ordered by solar zenith angle.
///
/// Sort key for the SZA-sorted workload ordering: processing cells in
/// increasing SZA (toward local noon) improves photolysis-cache locality
/// and load balance on device.
struct SZA_Pair {
  double sza; ///< Solar zenith angle [deg] (sort key).
  int index;  ///< Original cell index (payload carried through the sort).
  /// @brief Order ascending in SZA.
  KOKKOS_INLINE_FUNCTION bool operator<(const SZA_Pair &other) const {
    return sza < other.sza;
  }
};

/// @brief Per-run context: mechanism selection + reusable device scratch.
///
/// Created once by exaero_init() and held for the run's lifetime so the
/// sorter scratch is not reallocated every timestep.
struct ExaeroContext {
  std::string active_mechanism; ///< Registered chemistry mechanism name.
  int32_t total_grid_cells;     ///< Horizontal cell count (sorter league size).
  Kokkos::View<int *> sza_sorted_indices;  ///< Output permutation: device cell
                                           ///< order after SZA sort.
  Kokkos::View<SZA_Pair *> sorter_scratch; ///< (SZA, index) scratch for the
                                           ///< sort.

  /// @brief Allocate the device scratch views for @p cells grid cells.
  /// @param mech_name Active chemistry mechanism name.
  /// @param cells Total horizontal grid-cell count.
  ExaeroContext(const std::string &mech_name, int32_t cells)
      : active_mechanism(mech_name), total_grid_cells(cells),
        sza_sorted_indices("sza_sorted_indices", cells),
        sorter_scratch("sorter_scratch", cells) {}
};

} // namespace exaero
#endif // EXAERO_STATE_HPP
