#include "auxiliary_engines.hpp"
#include <iostream>
#include <stdexcept>

// T008: Declare extern "C" bindings strictly matching Fortran bind(c) signatures
#ifdef EXAERO_WITH_CLOUDJ
extern "C" void cloudj_driver_compute_jrates_device(const double* meteo_data, const int* extents, int cell_idx);
#endif

#ifdef EXAERO_WITH_ISORROPIALITE
extern "C" void isorropia_compute_thermo_device(double* conc_data, const double* meteo_data, const int* extents, int cell_idx);
#endif

namespace exaero {

KOKKOS_INLINE_FUNCTION
void AuxiliaryEngines::PhotolysisFunctor::operator()(const MemberType& team_member) const {
    const int sza_rank = team_member.league_rank();
    const int cell_idx = sorted_indices(sza_rank);

    if (!state.meteorology.data()) {
        Kokkos::abort("FATAL ERROR: state.meteorology.data() is null in PhotolysisFunctor");
    }

    // T012: Prevent thread starvation with sub-stepping calculations
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, 1), [&](const int& s) {
        // T010: Preprocessor conditional
#ifdef EXAERO_WITH_CLOUDJ
        // T009: Extract pure .data() pointer and explicit dimensions so Fortran can accept it via bind(c)
        int extents[4] = {
            static_cast<int>(state.meteorology.extent(0)),
            static_cast<int>(state.meteorology.extent(1)),
            static_cast<int>(state.meteorology.extent(2)),
            static_cast<int>(state.meteorology.extent(3))
        };
        cloudj_driver_compute_jrates_device(state.meteorology.data(), extents, cell_idx);
#endif
    });
}

KOKKOS_INLINE_FUNCTION
void AuxiliaryEngines::ThermoFunctor::operator()(const MemberType& team_member) const {
    const int cell_idx = sorted_indices(team_member.league_rank());

    if (!state.concentrations.data() || !state.meteorology.data()) {
        Kokkos::abort("FATAL ERROR: state.concentrations.data() or state.meteorology.data() is null in ThermoFunctor");
    }

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, 1), [&](const int& s) {
#ifdef EXAERO_WITH_ISORROPIALITE
        int extents[4] = {
            static_cast<int>(state.concentrations.extent(0)),
            static_cast<int>(state.concentrations.extent(1)),
            static_cast<int>(state.concentrations.extent(2)),
            static_cast<int>(state.concentrations.extent(3))
        };
        isorropia_compute_thermo_device(state.concentrations.data(), state.meteorology.data(), extents, cell_idx);
#endif
    });
}

void AuxiliaryEngines::compute_photolysis(ExaeroContext& ctx, UnmanagedDeviceState& state, double dt) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;

    if (!state.meteorology.data()) {
        throw std::runtime_error("FATAL ERROR: state.meteorology.data() is null");
    }

    int meteo_errors = 0;
    Kokkos::parallel_reduce("Validate_Meteo", state.meteorology.size(), KOKKOS_LAMBDA(const size_t i, int& lsum) {
        double val = state.meteorology.data()[i];
        if (val != val || val * 0.0 != 0.0) {
            lsum += 1;
        }
    }, meteo_errors);

    if (meteo_errors > 0) {
        throw std::runtime_error("FATAL ERROR: NaN/Inf detected in meteorology array");
    }

    Kokkos::parallel_for("CloudJ_Team_Dispatch",
        Kokkos::TeamPolicy<ExecSpace>(ctx.total_grid_cells, Kokkos::AUTO),
        PhotolysisFunctor(state, ctx.sza_sorted_indices)
    );
}

void AuxiliaryEngines::compute_thermodynamics(ExaeroContext& ctx, UnmanagedDeviceState& state, double dt) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;

    if (!state.concentrations.data() || !state.meteorology.data()) {
        throw std::runtime_error("FATAL ERROR: state.concentrations.data() or state.meteorology.data() is null");
    }

    int thermo_errors = 0;
    Kokkos::parallel_reduce("Validate_Thermo", state.concentrations.size(), KOKKOS_LAMBDA(const size_t i, int& lsum) {
        double val = state.concentrations.data()[i];
        if (val != val || val * 0.0 != 0.0) {
            lsum += 1;
        }
    }, thermo_errors);

    if (thermo_errors > 0) {
        throw std::runtime_error("FATAL ERROR: NaN/Inf detected in concentrations array");
    }

    Kokkos::parallel_for("ISORROPIA_Team_Dispatch",
        Kokkos::TeamPolicy<ExecSpace>(ctx.total_grid_cells, Kokkos::AUTO),
        ThermoFunctor(state, ctx.sza_sorted_indices)
    );
}

} // namespace exaero
