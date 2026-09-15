#ifndef EXAERO_AUXILIARY_ENGINES_HPP
#define EXAERO_AUXILIARY_ENGINES_HPP
/// @file auxiliary_engines.hpp
/// @brief Kokkos team-launch drivers for the external photolysis and
/// thermodynamics engines.
///
/// Wraps the CloudJ (photolysis J-rates) and ISORROPIA II (aerosol
/// thermodynamics) Fortran kernels behind Kokkos TeamPolicy launches, one
/// team per SZA-sorted grid cell. Both engines are optional: their calls
/// compile only when EXAERO_WITH_CLOUDJ / EXAERO_WITH_ISORROPIALITE are
/// defined, so the library builds without them.

#include "exaero/state.hpp"
#include <Kokkos_Core.hpp>

namespace exaero {

/// @brief Static facade that drives the optional external physics engines
/// (CloudJ photolysis, ISORROPIA II thermodynamics) over the SZA-sorted
/// cell decomposition owned by ExaeroContext.
class AuxiliaryEngines {
public:
  /// @brief Run CloudJ photolysis over the SZA-sorted cells.
  ///
  /// Validates the meteorology array for NaN/Inf (fail-fast), then
  /// launches one team per cell; each team calls the device J-rate driver
  /// for its sorted cell index.
  /// @param ctx Per-run context (supplies the SZA-sorted index order).
  /// @param state Unmanaged device state (meteorology in, J-rates out).
  /// @param dt Integration timestep [s].
  /// @throws std::runtime_error "FATAL ERROR:" on null or NaN/Inf input.
  static void compute_photolysis(ExaeroContext &ctx,
                                 UnmanagedDeviceState &state, double dt);

  /// @brief Run ISORROPIA II thermodynamic equilibration per cell.
  /// @param ctx Per-run context (supplies the SZA-sorted index order).
  /// @param state Unmanaged device state (concentrations + meteorology).
  /// @param dt Integration timestep [s].
  /// @throws std::runtime_error "FATAL ERROR:" on null or NaN/Inf input.
  static void compute_thermodynamics(ExaeroContext &ctx,
                                     UnmanagedDeviceState &state, double dt);

private:
  /// @brief Team functor: one CloudJ J-rate call per sorted cell.
  struct PhotolysisFunctor {
    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
    using MemberType = typename TeamPolicy::member_type;

    UnmanagedDeviceState state; ///< Non-owning device array views.
    Kokkos::View<int *>
        sorted_indices; ///< SZA permutation (league rank -> original cell).

    /// @brief Bind the state views and the sorted index permutation.
    PhotolysisFunctor(UnmanagedDeviceState s, Kokkos::View<int *> idx)
        : state(s), sorted_indices(idx) {}

    /// @brief Team body: compute this cell's photolysis J-rates.
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &team_member) const;
  };

  /// @brief Team functor: one ISORROPIA equilibration per sorted cell.
  struct ThermoFunctor {
    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
    using MemberType = typename TeamPolicy::member_type;

    UnmanagedDeviceState state; ///< Non-owning device array views.
    Kokkos::View<int *>
        sorted_indices; ///< SZA permutation (league rank -> original cell).

    /// @brief Bind the state views and the sorted index permutation.
    ThermoFunctor(UnmanagedDeviceState s, Kokkos::View<int *> idx)
        : state(s), sorted_indices(idx) {}

    /// @brief Team body: equilibrate this cell's aerosol inorganic mixture.
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &team_member) const;
  };
};

} // namespace exaero

#endif
