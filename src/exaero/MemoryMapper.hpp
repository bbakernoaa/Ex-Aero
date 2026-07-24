#pragma once
#include <Kokkos_Core.hpp>
#include <type_traits>

// Forward declare custom isolated layout types to avoid including any mdspan headers here
namespace exaero_mdspan {
    struct layout_right;
    struct layout_left;
}

namespace exaero {

    // Type trait to map exaero_mdspan layouts to Kokkos layouts
    template <typename Layout>
    struct mdspan_to_kokkos_layout;

    template <>
    struct mdspan_to_kokkos_layout<exaero_mdspan::layout_right> {
        using type = Kokkos::LayoutRight;
    };

    template <>
    struct mdspan_to_kokkos_layout<exaero_mdspan::layout_left> {
        using type = Kokkos::LayoutLeft;
    };

    // Generic unmanaged Kokkos View converter template (works for any conforming mdspan View)
    template <typename SpanType>
    auto make_unmanaged_kokkos_view(SpanType span) {
        using LayoutType = typename SpanType::layout_type;
        using KokkosLayout = typename mdspan_to_kokkos_layout<LayoutType>::type;
        
        using T = typename SpanType::element_type;
        using NonConstT = std::remove_const_t<T>;
        
        if constexpr (SpanType::rank() == 3) {
            using ViewType = std::conditional_t<std::is_const_v<T>, const NonConstT***, NonConstT***>;
            return Kokkos::View<ViewType, KokkosLayout, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                const_cast<NonConstT*>(span.data_handle()),
                span.extent(0),
                span.extent(1),
                span.extent(2)
            );
        } else {
            using ViewType = std::conditional_t<std::is_const_v<T>, const NonConstT**, NonConstT**>;
            return Kokkos::View<ViewType, KokkosLayout, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                const_cast<NonConstT*>(span.data_handle()),
                span.extent(0),
                span.extent(1)
            );
        }
    }

} // namespace exaero
