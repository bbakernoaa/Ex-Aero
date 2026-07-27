#pragma once
#include <Kokkos_Core.hpp>
#include <exaero/IAerosolPackage.hpp>
#include <type_traits>

namespace exaero {

// Type trait to map standard exaero_mdspan layouts to Kokkos layouts
template <typename Layout> struct mdspan_to_kokkos_layout;

template <> struct mdspan_to_kokkos_layout<exaero_mdspan::layout_right> {
  using type = Kokkos::LayoutRight;
};

template <> struct mdspan_to_kokkos_layout<exaero_mdspan::layout_left> {
  using type = Kokkos::LayoutLeft;
};

// Helper for 3D dynamic Views: mdspan -> unmanaged Kokkos::View
template <typename T, typename Layout, typename Accessor>
auto make_unmanaged_kokkos_view(
    exaero_mdspan::mdspan<
        T,
        exaero_mdspan::extents<size_t, std::dynamic_extent, std::dynamic_extent,
                               std::dynamic_extent>,
        Layout, Accessor>
        span) {
  using KokkosLayout = typename mdspan_to_kokkos_layout<Layout>::type;

  using NonConstT = std::remove_const_t<T>;
  using ViewType = std::conditional_t<std::is_const_v<T>, const NonConstT ***,
                                      NonConstT ***>;

  return Kokkos::View<ViewType, KokkosLayout, Kokkos::HostSpace,
                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      const_cast<NonConstT *>(span.data_handle()), span.extent(0),
      span.extent(1), span.extent(2));
}

// Helper for 2D dynamic Views: mdspan -> unmanaged Kokkos::View
template <typename T, typename Layout, typename Accessor>
auto make_unmanaged_kokkos_view(
    exaero_mdspan::mdspan<T,
                          exaero_mdspan::extents<size_t, std::dynamic_extent,
                                                 std::dynamic_extent>,
                          Layout, Accessor>
        span) {
  using KokkosLayout = typename mdspan_to_kokkos_layout<Layout>::type;

  using NonConstT = std::remove_const_t<T>;
  using ViewType =
      std::conditional_t<std::is_const_v<T>, const NonConstT **, NonConstT **>;

  return Kokkos::View<ViewType, KokkosLayout, Kokkos::HostSpace,
                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      const_cast<NonConstT *>(span.data_handle()), span.extent(0),
      span.extent(1));
}

} // namespace exaero
