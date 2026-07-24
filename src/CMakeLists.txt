# Public Interface Target
add_library(exaero SHARED dummy.cpp)

target_include_directories(exaero PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(exaero PUBLIC Kokkos::kokkos)

add_subdirectory(src_impl)
