/// @file auxiliary_engines.cpp
/// @brief Kokkos HOST team drivers bridging to the CloudJ / ISORROPIA Fortran
/// kernels.
///
/// Each engine call is guarded by its compile-time macro so the library
/// links without the optional dependencies. The team league runs over the
/// SZA-sorted cell order on Kokkos::DefaultHostExecutionSpace: one team per
/// cell, Kokkos::single(PerTeam) so the Fortran kernel (which is not
/// thread-parallel) executes once per cell, and pure .data() pointers plus an
/// explicit extents array are passed so the bind(c) Fortran signatures
/// receive contiguous arrays.
///
/// @note HOST-CONFINED by structural limitation: the wrapped kernels are host
/// bind(c) Fortran and cannot execute inside a device kernel. The functors
/// are therefore NOT annotated KOKKOS_INLINE_FUNCTION, every launch is pinned
/// to the host execution space, and configuring a device backend together
/// with an engine option is a compile error (see the guards below) rather
/// than a silent miscompile. When an engine option is OFF, the entry points
/// throw instead of launching an empty team: no silent no-op compute path.
///
/// @note Defensive validation: both entry points scan the input arrays for
/// NaN/Inf (val != val || val * 0.0 != 0.0) via a host parallel_reduce and
/// throw "FATAL ERROR:" before any kernel launch, so corrupted host data
/// never reaches the physics.
#include "auxiliary_engines.hpp"
#include <iostream>
#include <stdexcept>

// Declare extern "C" bindings strictly matching the Fortran bind(c)
// signatures. These are HOST kernels: no name here carries a _device suffix,
// because no device definition exists (or could exist without porting the
// Fortran itself). The names must match what the optional engine library
// provides.
#ifdef EXAERO_WITH_CLOUDJ
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
#error "EXAERO_WITH_CLOUDJ is host-confined: the CloudJ bind(c) kernel cannot run inside a device launch. Port the Fortran to device (Tier-2) before enabling it in a GPU build."
#endif
extern "C" void cloudj_driver_compute_jrates(const double *meteo_data,
                                             const int *extents,
                                             int cell_idx);
#endif

#ifdef EXAERO_WITH_ISORROPIALITE
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
#error "EXAERO_WITH_ISORROPIALITE is host-confined: the ISORROPIA bind(c) kernel cannot run inside a device launch. Port the Fortran to device (Tier-2) before enabling it in a GPU build."
#endif
extern "C" void isorropia_compute_thermo(double *conc_data,
                                         const double *meteo_data,
                                         const int *extents,
                                         int cell_idx);
#endif

namespace exaero {

// @copydoc AuxiliaryEngines::PhotolysisFunctor::operator()
// Host-only: intentionally no KOKKOS_INLINE_FUNCTION, the body calls a host
// bind(c) Fortran kernel through zero-copy host pointers.
void AuxiliaryEngines::PhotolysisFunctor::operator()(
    const MemberType &team_member) const {
  const int sza_rank = team_member.league_rank();
  const int cell_idx = sorted_indices(sza_rank);

  if (!state.meteorology.data()) {
    Kokkos::abort(
        "FATAL ERROR: state.meteorology.data() is null in PhotolysisFunctor");
  }

  // : Prevent thread starvation with sub-stepping calculations
  Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
    // : Extract pure.data() pointer and explicit dimensions so Fortran can
    // accept it via bind(c)
#ifdef EXAERO_WITH_CLOUDJ
    int extents[4] = {static_cast<int>(state.meteorology.extent(0)),
                      static_cast<int>(state.meteorology.extent(1)),
                      static_cast<int>(state.meteorology.extent(2)),
                      static_cast<int>(state.meteorology.extent(3))};
    cloudj_driver_compute_jrates(state.meteorology.data(), extents, cell_idx);
#endif
  });
}

// @copydoc AuxiliaryEngines::ThermoFunctor::operator()
// Host-only: intentionally no KOKKOS_INLINE_FUNCTION, the body calls a host
// bind(c) Fortran kernel through zero-copy host pointers.
void AuxiliaryEngines::ThermoFunctor::operator()(
    const MemberType &team_member) const {
  const int cell_idx = sorted_indices(team_member.league_rank());

  if (!state.concentrations.data() || !state.meteorology.data()) {
    Kokkos::abort("FATAL ERROR: state.concentrations.data() or "
                  "state.meteorology.data() is null in ThermoFunctor");
  }

  Kokkos::single(Kokkos::PerTeam(team_member), [&]() {
#ifdef EXAERO_WITH_ISORROPIALITE
    int extents[4] = {static_cast<int>(state.concentrations.extent(0)),
                      static_cast<int>(state.concentrations.extent(1)),
                      static_cast<int>(state.concentrations.extent(2)),
                      static_cast<int>(state.concentrations.extent(3))};
    isorropia_compute_thermo(state.concentrations.data(),
                             state.meteorology.data(), extents, cell_idx);
#endif
  });
}

// @copydoc AuxiliaryEngines::compute_photolysis
void AuxiliaryEngines::compute_photolysis(ExaeroContext &ctx,
                                          UnmanagedDeviceState &state,
                                          double dt) {
#ifndef EXAERO_WITH_CLOUDJ
  // Fail loud instead of launching an empty team: with the option OFF the
  // functor body compiles to nothing, and a silent no-op "photolysis" would
  // misreport physics that never ran. Photolysis J-rates are an MKPP *input*
  // (mkpp_set_photolysis_ptrs); this driver exists only for builds that link
  // a real CloudJ engine.
  (void)ctx;
  (void)state;
  (void)dt;
  throw std::runtime_error(
      "FATAL ERROR: CloudJ photolysis driver not compiled in "
      "(EXAERO_WITH_CLOUDJ=OFF); compute_photolysis would be a silent no-op. "
      "Supply J-rates through the MKPP photolysis interface instead.");
#else
  // Host-confined dispatch (see file header): the Fortran kernel and the
  // zero-copy host pointers behind the unmanaged views cannot be touched from
  // a device launch. Validation runs on the host for the same reason.
  using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

  if (!state.meteorology.data()) {
    throw std::runtime_error("FATAL ERROR: state.meteorology.data() is null");
  }

  int meteo_errors = 0;
  Kokkos::parallel_reduce(
      "Validate_Meteo",
      Kokkos::RangePolicy<HostExecSpace>(0, state.meteorology.size()),
      [=](const size_t i, int &lsum) {
        double val = state.meteorology.data()[i];
        if (val != val || val * 0.0 != 0.0) {
          lsum += 1;
        }
      },
      meteo_errors);

  if (meteo_errors > 0) {
    throw std::runtime_error(
        "FATAL ERROR: NaN/Inf detected in meteorology array");
  }

  Kokkos::parallel_for(
      "CloudJ_Team_Dispatch",
      Kokkos::TeamPolicy<HostExecSpace>(ctx.total_grid_cells, Kokkos::AUTO),
      PhotolysisFunctor(state, ctx.sza_sorted_indices));
  (void)dt;
#endif
}

// @copydoc AuxiliaryEngines::compute_thermodynamics
void AuxiliaryEngines::compute_thermodynamics(ExaeroContext &ctx,
                                              UnmanagedDeviceState &state,
                                              double dt) {
#ifndef EXAERO_WITH_ISORROPIALITE
  // Fail loud instead of launching an empty team (same rationale as
  // compute_photolysis): the ISORROPIA coupling stays host-confined until the
  // Tier-2 port lands; with the option OFF nothing here could execute.
  (void)ctx;
  (void)state;
  (void)dt;
  throw std::runtime_error(
      "FATAL ERROR: ISORROPIA thermodynamic driver not compiled in "
      "(EXAERO_WITH_ISORROPIALITE=OFF); compute_thermodynamics would be a "
      "silent no-op.");
#else
  using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

  if (!state.concentrations.data() || !state.meteorology.data()) {
    throw std::runtime_error("FATAL ERROR: state.concentrations.data() or "
                             "state.meteorology.data() is null");
  }

  int thermo_errors = 0;
  Kokkos::parallel_reduce(
      "Validate_Thermo",
      Kokkos::RangePolicy<HostExecSpace>(0, state.concentrations.size()),
      [=](const size_t i, int &lsum) {
        double val = state.concentrations.data()[i];
        if (val != val || val * 0.0 != 0.0) {
          lsum += 1;
        }
      },
      thermo_errors);

  if (thermo_errors > 0) {
    throw std::runtime_error(
        "FATAL ERROR: NaN/Inf detected in concentrations array");
  }

  Kokkos::parallel_for(
      "ISORROPIA_Team_Dispatch",
      Kokkos::TeamPolicy<HostExecSpace>(ctx.total_grid_cells, Kokkos::AUTO),
      ThermoFunctor(state, ctx.sza_sorted_indices));
  (void)dt;
#endif
}

} // namespace exaero
