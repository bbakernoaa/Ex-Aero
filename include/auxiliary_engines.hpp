#ifndef EXAERO_AUXILIARY_ENGINES_HPP
#define EXAERO_AUXILIARY_ENGINES_HPP
/// @file auxiliary_engines.hpp
/// @brief Kokkos HOST team-launch drivers for the external photolysis and
/// thermodynamics engines.
///
/// Wraps the CloudJ (photolysis J-rates) and ISORROPIA II (aerosol
/// thermodynamics) Fortran kernels behind Kokkos host TeamPolicy launches, one
/// team per SZA-sorted grid cell. Both engines are optional: their calls
/// compile only when EXAERO_WITH_CLOUDJ / EXAERO_WITH_ISORROPIALITE are
/// defined, so the library builds without them.
///
/// @note HOST-CONFINED. The wrapped kernels are host `bind(c)` Fortran that
/// cannot be dereferenced inside a device kernel, so the dispatch runs on
/// @c Kokkos::DefaultHostExecutionSpace and the functors are NOT annotated
/// @c KOKKOS_INLINE_FUNCTION. This is a deliberate honesty constraint (see the
/// device-callability status contract): the chain is `host-confined`, never
/// advertised as device-executable. Callers on a GPU build must keep the state
/// host-accessible.

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
  /// launches one host team per cell; each team calls the host J-rate driver
  /// for its sorted cell index.
  /// @param ctx Per-run context (supplies the SZA-sorted index order).
  /// @param state Unmanaged state (meteorology in, J-rates out).
  /// @param dt Integration timestep [s].
  /// @throws std::runtime_error "FATAL ERROR:" on null/NaN/Inf input, or when
  ///         the CloudJ engine was not compiled in (no silent no-op).
  /// @note host-confined dispatch; see the file-level note.
  static void compute_photolysis(ExaeroContext &ctx,
                                 UnmanagedDeviceState &state, double dt);

  /// @brief Run ISORROPIA II thermodynamic equilibration per cell.
  /// @param ctx Per-run context (supplies the SZA-sorted index order).
  /// @param state Unmanaged state (concentrations + meteorology).
  /// @param dt Integration timestep [s].
  /// @throws std::runtime_error "FATAL ERROR:" on null/NaN/Inf input, or when
  ///         the ISORROPIA engine was not compiled in (no silent no-op).
  /// @note host-confined dispatch; see the file-level note.
  static void compute_thermodynamics(ExaeroContext &ctx,
                                     UnmanagedDeviceState &state, double dt);

private:
  /// @brief Host team functor: one CloudJ J-rate call per sorted cell.
  struct PhotolysisFunctor {
    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultHostExecutionSpace>;
    using MemberType = typename TeamPolicy::member_type;

    UnmanagedDeviceState state; ///< Non-owning host-accessible array views.
    Kokkos::View<int *>
        sorted_indices; ///< SZA permutation (league rank -> original cell).

    /// @brief Bind the state views and the sorted index permutation.
    PhotolysisFunctor(UnmanagedDeviceState s, Kokkos::View<int *> idx)
        : state(s), sorted_indices(idx) {}

    /// @brief Team body: compute this cell's photolysis J-rates.
    /// Host-only: calls a host bind(c) Fortran kernel, so it is intentionally
    /// NOT annotated KOKKOS_INLINE_FUNCTION.
    void operator()(const MemberType &team_member) const;
  };

  /// @brief Host team functor: one ISORROPIA equilibration per sorted cell.
  struct ThermoFunctor {
    using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultHostExecutionSpace>;
    using MemberType = typename TeamPolicy::member_type;

    UnmanagedDeviceState state; ///< Non-owning host-accessible array views.
    Kokkos::View<int *>
        sorted_indices; ///< SZA permutation (league rank -> original cell).

    /// @brief Bind the state views and the sorted index permutation.
    ThermoFunctor(UnmanagedDeviceState s, Kokkos::View<int *> idx)
        : state(s), sorted_indices(idx) {}

    /// @brief Team body: equilibrate this cell's aerosol inorganic mixture.
    /// Host-only bind(c) Fortran call; intentionally NOT device-annotated.
    void operator()(const MemberType &team_member) const;
  };
};

} // namespace exaero

#endif
