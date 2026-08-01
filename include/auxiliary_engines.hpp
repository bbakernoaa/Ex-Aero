#ifndef EXAERO_AUXILIARY_ENGINES_HPP
#define EXAERO_AUXILIARY_ENGINES_HPP

#include <Kokkos_Core.hpp>
#include "exaero/state.hpp"

namespace exaero {

class AuxiliaryEngines {
public:
    static void compute_photolysis(ExaeroContext& ctx, UnmanagedDeviceState& state, double dt);
    static void compute_thermodynamics(ExaeroContext& ctx, UnmanagedDeviceState& state, double dt);
    
private:
    struct PhotolysisFunctor {
        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using MemberType = typename TeamPolicy::member_type;
        
        UnmanagedDeviceState state;
        Kokkos::View<int*> sorted_indices;
        
        PhotolysisFunctor(UnmanagedDeviceState s, Kokkos::View<int*> idx) : state(s), sorted_indices(idx) {}

        KOKKOS_INLINE_FUNCTION
        void operator()(const MemberType& team_member) const;
    };
    
    struct ThermoFunctor {
        using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
        using MemberType = typename TeamPolicy::member_type;
        
        UnmanagedDeviceState state;
        Kokkos::View<int*> sorted_indices;
        
        ThermoFunctor(UnmanagedDeviceState s, Kokkos::View<int*> idx) : state(s), sorted_indices(idx) {}

        KOKKOS_INLINE_FUNCTION
        void operator()(const MemberType& team_member) const;
    };
};

} // namespace exaero

#endif
