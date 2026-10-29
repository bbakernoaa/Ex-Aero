# Adopt a yaml-cpp the parent build already provides instead of re-running
# find_package: config-mode discovery would re-create the imported target and
# hard-fail ("already exists") when the parent vendored it via FetchContent.
# Same guard pattern as Kokkos in the top-level CMakeLists.
if(NOT TARGET yaml-cpp::yaml-cpp AND NOT TARGET yaml-cpp)
  find_package(yaml-cpp REQUIRED)
endif()

# --- Build-time baked-in MIE table generation ---
# The pinned GEOSmie RRTMG snapshot is converted into a deterministic float64
# header. Listing the generated header as a source makes every consumer depend
# on the rule.
set(MIE_GENERATED_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/data/generated/MieTableData.hpp)
set(MIE_GENERATOR ${EXAERO_REPO_ROOT}/tools/generate_mie_tables.py)
set(MIE_SNAPSHOT_DIR ${EXAERO_REPO_ROOT}/tools/geosmie_snapshot)

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

# Modal curves are generated offline with the pinned GEOSmie Mie engine. Normal
# builds only serialize the committed snapshot; no GEOSmie or network required.
set(MODAL_MIE_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/data/generated/ModalOpticsData.hpp)
add_custom_command(
  OUTPUT ${MODAL_MIE_HEADER}
  COMMAND
    Python3::Interpreter ${EXAERO_REPO_ROOT}/tools/generate_modal_optics.py
    --snapshot ${EXAERO_REPO_ROOT}/tools/modal_optics/modal_mie_v1.nc --header
    ${MODAL_MIE_HEADER}
  DEPENDS ${EXAERO_REPO_ROOT}/tools/generate_modal_optics.py
          ${EXAERO_REPO_ROOT}/tools/modal_optics/modal_mie_v1.nc
  COMMENT "Baking committed GEOSmie modal optical curves"
  VERBATIM)
add_custom_target(exaero_modal_mie_data DEPENDS ${MODAL_MIE_HEADER})

set(MAM4XX_MIE_HEADER
    ${CMAKE_CURRENT_SOURCE_DIR}/data/generated/Mam4xxOpticsData.hpp)
add_custom_command(
  OUTPUT ${MAM4XX_MIE_HEADER}
  COMMAND
    Python3::Interpreter ${EXAERO_REPO_ROOT}/tools/generate_modal_optics.py
    --snapshot ${EXAERO_REPO_ROOT}/tools/modal_optics/mam4xx_mie_v1.nc --header
    ${MAM4XX_MIE_HEADER} --namespace mie_mam4xx
  DEPENDS ${EXAERO_REPO_ROOT}/tools/generate_modal_optics.py
          ${EXAERO_REPO_ROOT}/tools/modal_optics/mam4xx_mie_v1.nc
  COMMENT "Baking committed GEOSmie MAM4xx component curves"
  VERBATIM)
add_custom_target(exaero_mam4xx_mie_data DEPENDS ${MAM4XX_MIE_HEADER})

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
add_library(exaero_solver STATIC gocart/GocartSolver.cpp modal/ModalSolver.cpp
                                 modal/ModalChemistrySolver.cpp)
target_include_directories(exaero_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                                 ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(exaero_solver PUBLIC Kokkos::kokkos)

# 1. The Package Interface Implementation (mdspan, yaml-cpp, NO Kokkos)
add_library(
  exaero_impl OBJECT
  PackageFactory.cpp
  HostCoupling.cpp
  Removal.cpp
  PopulationProperties.cpp
  PopulationProperties_C.cpp
  Removal_C.cpp
  HostCoupling_C.cpp
  CloudIceActivation.cpp
  CloudIceCoupling.cpp
  CommonCloudIce_C.cpp
  CloudIceCoupling_C.cpp
  modal/ModalPackage.cpp
  modal/ModalChemistry.cpp
  modal/ModalChemistry_C.cpp
  modal/ModalRemoval.cpp
  modal/ModalCloudIce.cpp
  modal/ModalRemoval_C.cpp
  modal/ModalPresets.cpp
  modal/ModalCoupling_C.cpp
  gocart/GocartPackage.cpp
  gocart/BulkRemoval.cpp
  gocart/BulkCloudIce.cpp
  gocart/BulkRemoval_C.cpp
  utils/LutGenerator.cpp
  utils/DubovikSpheroidDatabase.cpp
  loader/MieTableStore.cpp
  loader/MieFileLoader.cpp
  loader/Provenance.cpp)
add_dependencies(exaero_impl exaero_mie_data)
add_dependencies(exaero_impl exaero_modal_mie_data exaero_mam4xx_mie_data)
target_include_directories(exaero_impl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                               ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(exaero_impl PUBLIC exaero_mdspan yaml-cpp::yaml-cpp)

# Link exaero public library to private implementations and all external library
# dependencies
target_link_libraries(exaero PRIVATE exaero_impl exaero_solver
                                     yaml-cpp::yaml-cpp Kokkos::kokkos)
target_link_libraries(exaero PUBLIC exaero_mdspan)

if(EXAERO_ENABLE_MAM4XX)
  include(${EXAERO_REPO_ROOT}/cmake/Mam4xx.cmake)
  target_sources(exaero_impl PRIVATE mam4xx/Mam4xxPackage.cpp)
  target_link_libraries(exaero_solver PRIVATE exaero_mam4xx_dependency)
  target_compile_definitions(exaero_solver PRIVATE EXAERO_ENABLE_MAM4XX=1)
  target_compile_definitions(exaero_impl PRIVATE EXAERO_ENABLE_MAM4XX=1)
  add_library(exaero_mam4xx_solver STATIC mam4xx/Mam4xxSolver.cpp)
  target_include_directories(
    exaero_mam4xx_solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
                                 ${EXAERO_REPO_ROOT}/src)
  target_link_libraries(exaero_mam4xx_solver PRIVATE ${EXAERO_MAM4XX_TARGET})
  target_link_libraries(exaero PRIVATE exaero_mam4xx_solver)
endif()
