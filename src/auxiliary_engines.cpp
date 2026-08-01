#include "auxiliary_engines.hpp"
#include <iostream>

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
    Kokkos::parallel_for("CloudJ_Team_Dispatch",
        Kokkos::TeamPolicy<ExecSpace>(ctx.total_grid_cells, Kokkos::AUTO),
        PhotolysisFunctor(state, ctx.sza_sorted_indices)
    );
}

void AuxiliaryEngines::compute_thermodynamics(ExaeroContext& ctx, UnmanagedDeviceState& state, double dt) {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    Kokkos::parallel_for("ISORROPIA_Team_Dispatch",
        Kokkos::TeamPolicy<ExecSpace>(ctx.total_grid_cells, Kokkos::AUTO),
        ThermoFunctor(state, ctx.sza_sorted_indices)
    );
}

} // namespace exaero
