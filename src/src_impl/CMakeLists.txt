find_package(yaml-cpp REQUIRED)

# --- Build-time baked-in MIE table generation ---
# The pinned GEOSmie RRTMG snapshot is converted into a deterministic float64
# header. Listing the generated header as a source makes every consumer depend
# on the rule.
set(MIE_GENERATED_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/data/generated/MieTableData.hpp)
set(MIE_GENERATOR ${CMAKE_SOURCE_DIR}/tools/generate_mie_tables.py)
set(MIE_SNAPSHOT_DIR ${CMAKE_SOURCE_DIR}/tools/geosmie_snapshot)

add_custom_command(
  OUTPUT ${MIE_GENERATED_HEADER}
  COMMAND Python3::Interpreter ${MIE_GENERATOR} ${MIE_GENERATED_HEADER}
          ${MIE_SNAPSHOT_DIR}
  DEPENDS ${MIE_GENERATOR}
          ${MIE_SNAPSHOT_DIR}/opticsBands_DU.v15_3.RRTMG.nc
          ${MIE_SNAPSHOT_DIR}/opticsBands_SS.v3_3.RRTMG.nc
          ${MIE_SNAPSHOT_DIR}/opticsBands_SU.v1_3.RRTMG.nc
          ${MIE_SNAPSHOT_DIR}/opticsBands_BC.v1_3.RRTMG.nc
          ${MIE_SNAPSHOT_DIR}/opticsBands_OC.v1_3.RRTMG.nc
          ${MIE_SNAPSHOT_DIR}/opticsBands_NI.v2_5.RRTMG.nc
  COMMENT "Generating baked-in MieTableData.hpp from pinned GEOSmie snapshot"
  VERBATIM)

add_custom_target(exaero_mie_data DEPENDS ${MIE_GENERATED_HEADER})

# Private Header Target for mdspan backport
add_library(exaero_mdspan INTERFACE)
target_include_directories(
  exaero_mdspan
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/mdspan/include>)

# The vendored mdspan auto-selects operator() only when the bracket operator is
# unavailable (pre-C++23). The project standard is now C++23, where
# __cpp_multidimensional_subscript is defined and the auto-selection flips to
# operator[] only, breaking every EX-aero call site that uses view(i, j, k).
# Force the paren operator ON (the header's #ifndef guard leaves [] available
# too) so the public mdspan vocabulary keeps its documented () convention across
# all targets and standards.
target_compile_definitions(exaero_mdspan INTERFACE MDSPAN_USE_PAREN_OPERATOR=1)

# 1. The GPU Solver (Kokkos only, NO private mdspan headers)
add_library(exaero_solver STATIC gocart/GocartSolver.cpp)
target_include_directories(exaero_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                                 ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(exaero_solver PUBLIC Kokkos::kokkos)

# 1. The Package Interface Implementation (mdspan, yaml-cpp, NO Kokkos)
add_library(
  exaero_impl OBJECT
  gocart/GocartPackage.cpp utils/LutGenerator.cpp
  utils/DubovikSpheroidDatabase.cpp loader/MieTableStore.cpp
  loader/MieFileLoader.cpp loader/Provenance.cpp)
add_dependencies(exaero_impl exaero_mie_data)
target_include_directories(exaero_impl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                               ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(exaero_impl PUBLIC exaero_mdspan yaml-cpp::yaml-cpp)

# Link exaero public library to private implementations and all external library
# dependencies
target_link_libraries(exaero PRIVATE exaero_impl exaero_solver
                                     yaml-cpp::yaml-cpp Kokkos::kokkos)
target_link_libraries(exaero PUBLIC exaero_mdspan)
