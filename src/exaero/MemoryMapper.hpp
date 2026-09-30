#pragma once
/// @file MemoryMapper.hpp
/// @brief Zero-copy bridge from public mdspan views to Kokkos::View.
///
/// Internal (Kokkos-bearing) helper: converts the non-owning
/// @c exaero_mdspan views that cross the public boundary into unmanaged
/// @c Kokkos::View objects for device kernels — no allocation, no
/// transpose. The @c layout_left / @c LayoutLeft pairing preserves the
/// Fortran column-major contract so a host array is read with stride-1
/// locality in the inner (cell) dimension.
///
/// @note This header includes <Kokkos_Core.hpp> and is therefore for
/// private/implementation use only; it must never appear in a public
/// include chain (ADR-001).
#include <Kokkos_Core.hpp>
#include <exaero/IAerosolPackage.hpp>
#include <type_traits>

namespace exaero {

/// @brief Type trait mapping an @c exaero_mdspan layout tag to its
/// equivalent Kokkos layout policy.
/// @tparam Layout The mdspan layout tag to translate.
template <typename Layout> struct mdspan_to_kokkos_layout;

/// @brief Row-major (C-order) mdspan layout maps to Kokkos::LayoutRight.
template <> struct mdspan_to_kokkos_layout<exaero_mdspan::layout_right> {
  using type = Kokkos::LayoutRight; ///< The equivalent Kokkos layout policy.
};

/// @brief Column-major (Fortran-order) mdspan layout maps to
/// Kokkos::LayoutLeft — the default for grid state arrays.
template <> struct mdspan_to_kokkos_layout<exaero_mdspan::layout_left> {
  using type = Kokkos::LayoutLeft; ///< The equivalent Kokkos layout policy.
};

/// @brief Wrap a 3D dynamic mdspan as an unmanaged Kokkos::View (zero-copy).
///
/// The result aliases the mdspan's data handle directly (no ownership,
/// HostSpace); const @c T yields a const-qualified view type. Used to
/// hand (cell, level, species) state to device kernels without copying.
/// @tparam T Element type (may be const).
/// @tparam Layout mdspan layout tag (resolved to the Kokkos counterpart).
/// @param span The 3D non-owning view to alias.
/// @return An unmanaged 3D Kokkos::View over the same memory.
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

/// @brief Wrap a 2D dynamic mdspan as an unmanaged Kokkos::View (zero-copy).
///
/// The 2D overload of make_unmanaged_kokkos_view() for environmental
/// fields laid out as (cell, level).
/// @tparam T Element type (may be const).
/// @tparam Layout mdspan layout tag (resolved to the Kokkos counterpart).
/// @param span The 2D non-owning view to alias.
/// @return An unmanaged 2D Kokkos::View over the same memory.
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
