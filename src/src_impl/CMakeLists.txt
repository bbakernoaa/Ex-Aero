find_package(yaml-cpp REQUIRED)

# Private Header Target for mdspan backport
add_library(exaero_mdspan INTERFACE)
target_include_directories(exaero_mdspan INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/mdspan/include>
)

# 1. The GPU Solver (Kokkos only, NO private mdspan headers)
add_library(exaero_solver STATIC
    gocart/GocartSolver.cpp
)
target_include_directories(exaero_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(exaero_solver PUBLIC Kokkos::kokkos)

# 2. The Package Interface Implementation (mdspan, yaml-cpp, NO Kokkos)
add_library(exaero_impl OBJECT
    gocart/GocartPackage.cpp
)
target_include_directories(exaero_impl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(exaero_impl PUBLIC exaero_mdspan yaml-cpp::yaml-cpp)

# Link exaero public library to private implementations and all external library dependencies
target_link_libraries(exaero PRIVATE exaero_impl exaero_solver yaml-cpp::yaml-cpp Kokkos::kokkos)
target_link_libraries(exaero PUBLIC exaero_mdspan)
