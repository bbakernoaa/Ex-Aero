# Public Interface Target (builds C++, C, and Fortran interfaces)
add_library(exaero SHARED
    dummy.cpp
    exaero/IAerosolPackage_C.cpp
    exaero/exaero_interface.f90
)

target_include_directories(exaero PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(exaero PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src_impl
)

target_link_libraries(exaero PUBLIC Kokkos::kokkos)

add_subdirectory(src_impl)
