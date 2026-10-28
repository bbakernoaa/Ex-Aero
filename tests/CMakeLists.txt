# GOCART Diagnostics and YAML Parsing Test Runner
find_package(Catch2 2 CONFIG REQUIRED)
add_executable(exaero_cloud_ice_config_test test_cloud_ice_config.cpp)
target_link_libraries(exaero_cloud_ice_config_test PRIVATE exaero_mdspan)
target_include_directories(exaero_cloud_ice_config_test PRIVATE ${EXAERO_REPO_ROOT}/src)
target_compile_features(exaero_cloud_ice_config_test PRIVATE cxx_std_23)
add_test(NAME ExaeroCloudIceConfig COMMAND exaero_cloud_ice_config_test)
add_executable(exaero_modal_test test_modal.cpp test_modal_mam4.cpp
                                 test_modal_chemistry.cpp test_modal_main.cpp)
add_executable(exaero_modal_removal_math_test test_modal_removal_math.cpp)
add_test(NAME ExaeroModalRemovalMath COMMAND exaero_modal_removal_math_test)
target_include_directories(exaero_modal_removal_math_test PRIVATE ${EXAERO_REPO_ROOT}/src)
target_compile_features(exaero_modal_removal_math_test PRIVATE cxx_std_23)
add_executable(exaero_modal_removal_test test_modal_removal.cpp)
target_link_libraries(exaero_modal_removal_test PRIVATE exaero)
target_compile_features(exaero_modal_removal_test PRIVATE cxx_std_23)
add_test(NAME ExaeroModalRemoval COMMAND exaero_modal_removal_test)
add_executable(exaero_bulk_removal_test test_bulk_removal.cpp)
target_link_libraries(exaero_bulk_removal_test PRIVATE exaero)
target_compile_features(exaero_bulk_removal_test PRIVATE cxx_std_23)
add_test(NAME ExaeroBulkRemoval COMMAND exaero_bulk_removal_test)
add_executable(exaero_bulk_cloud_ice_test test_bulk_cloud_ice.cpp)
target_link_libraries(exaero_bulk_cloud_ice_test PRIVATE exaero)
target_compile_features(exaero_bulk_cloud_ice_test PRIVATE cxx_std_23)
add_test(NAME ExaeroBulkCloudIce COMMAND exaero_bulk_cloud_ice_test)
target_link_libraries(exaero_modal_test PRIVATE exaero Catch2::Catch2)
target_compile_definitions(
  exaero_modal_test
  PRIVATE EXAERO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
          EXAERO_TEST_SOURCE_DIR="${EXAERO_REPO_ROOT}")
add_test(NAME ExaeroModal COMMAND exaero_modal_test)
set_tests_properties(
  ExaeroModal PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false")
add_test(
  NAME ExaeroModalMieDeterminism
  COMMAND
    ${Python3_EXECUTABLE} ${EXAERO_REPO_ROOT}/tools/generate_modal_optics.py
    --snapshot ${EXAERO_REPO_ROOT}/tools/modal_optics/modal_mie_v1.nc --header
    ${EXAERO_REPO_ROOT}/src/src_impl/data/generated/ModalOpticsData.hpp --check)
add_test(
  NAME ExaeroMam4xxMieDeterminism
  COMMAND
    ${Python3_EXECUTABLE} ${EXAERO_REPO_ROOT}/tools/generate_modal_optics.py
    --snapshot ${EXAERO_REPO_ROOT}/tools/modal_optics/mam4xx_mie_v1.nc --header
    ${EXAERO_REPO_ROOT}/src/src_impl/data/generated/Mam4xxOpticsData.hpp
    --namespace mie_mam4xx --check)
add_executable(exaero_modal_fortran_test test_modal_fortran.f90)
target_link_libraries(exaero_modal_fortran_test PRIVATE exaero)
add_test(NAME ExaeroModalFortran COMMAND exaero_modal_fortran_test)
set_tests_properties(
  ExaeroModalFortran PROPERTIES ENVIRONMENT
                                "OMP_NUM_THREADS=2;OMP_PROC_BIND=false")

add_executable(exaero_package_selection_test test_package_selection.cpp)
target_link_libraries(exaero_package_selection_test PRIVATE exaero)
add_test(NAME ExaeroPackageSelection COMMAND exaero_package_selection_test)
add_executable(exaero_backend_fortran_test test_backend_fortran.f90)
target_link_libraries(exaero_backend_fortran_test PRIVATE exaero)
add_test(NAME ExaeroBackendFortran COMMAND exaero_backend_fortran_test)

if(EXAERO_ENABLE_MAM4XX)
  find_package(Catch2 CONFIG REQUIRED)
  add_executable(exaero_mam4xx_test test_mam4xx.cpp)
  target_link_libraries(exaero_mam4xx_test
                        PRIVATE exaero exaero_mam4xx_dependency Catch2::Catch2)
  add_test(NAME ExaeroMam4xx COMMAND exaero_mam4xx_test)
  set_tests_properties(
    ExaeroMam4xx PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false")

  add_executable(exaero_modal_mam4_comparison modal_mam4_comparison.cpp
                                              modal_mam4_reference.cpp)
  target_link_libraries(exaero_modal_mam4_comparison
                        PRIVATE exaero exaero_mam4xx_dependency)
  target_compile_features(exaero_modal_mam4_comparison PRIVATE cxx_std_23)
  add_test(NAME ExaeroModalMam4Comparison COMMAND exaero_modal_mam4_comparison
                                                  --check)
  set_tests_properties(
    ExaeroModalMam4Comparison
    PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2;OMP_PROC_BIND=false" TIMEOUT 60)

  option(EXAERO_COMPARISON_MKPP
         "Compile and execute the comparison YAML with the installed MKPP" OFF)
  if(EXAERO_COMPARISON_MKPP)
    set(comparison_generated "${CMAKE_CURRENT_BINARY_DIR}/comparison_generated")
    set(comparison_rates
        "${comparison_generated}/modal_mam4_comparison/rates.cpp")
    file(MAKE_DIRECTORY "${comparison_generated}")
    configure_file(${EXAERO_REPO_ROOT}/tools/modal_comparison/reactions.yaml
                   "${comparison_generated}/reactions.yaml" COPYONLY)
    configure_file(${EXAERO_REPO_ROOT}/tools/modal_comparison/environment.yaml
                   "${comparison_generated}/environment.yaml" COPYONLY)
    add_custom_command(
      OUTPUT
        "${comparison_rates}"
        "${comparison_generated}/modal_mam4_comparison/host_declarations.hpp"
      COMMAND ${Python3_EXECUTABLE} -m mkpp.cli compile reactions.yaml --out .
              --test-env environment.yaml
      DEPENDS ${EXAERO_REPO_ROOT}/tools/modal_comparison/reactions.yaml
              ${EXAERO_REPO_ROOT}/tools/modal_comparison/environment.yaml
      COMMENT "Compile modal MAM4 comparison reaction YAML with MKPP"
      WORKING_DIRECTORY "${comparison_generated}"
      VERBATIM)
    target_sources(exaero_modal_mam4_comparison PRIVATE "${comparison_rates}")
    target_include_directories(exaero_modal_mam4_comparison
                               PRIVATE "${comparison_generated}")
    target_compile_definitions(
      exaero_modal_mam4_comparison
      PRIVATE EXAERO_COMPARISON_MKPP=1
              MKPP_SEPARATE_HOST_KERNELS_MODAL_MAM4_COMPARISON=1)
  endif()
endif()

add_executable(exaero_test_runner test_main.cpp)
target_include_directories(exaero_test_runner
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl)
# Absolute fixture directory injected at configure time so runtime-file tests
# resolve tests/data/ hermetically regardless of the CTest working directory.
target_compile_definitions(
  exaero_test_runner
  PRIVATE EXAERO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(exaero_test_runner PRIVATE exaero Kokkos::kokkos)
add_test(NAME ExaeroUnitTests COMMAND exaero_test_runner)

# Standalone public-header compile gate (invariant ): compiles the new public
# attribute surface with NO Kokkos include path and links NOTHING (OBJECT
# library - the gate is compilation itself). If any public header pulled
# <Kokkos_Core.hpp> in transitively, this target fails to build.
add_library(exaero_standalone_public_check OBJECT standalone_public_check.cpp)
target_compile_features(exaero_standalone_public_check PRIVATE cxx_std_23)
target_include_directories(exaero_standalone_public_check
                           PRIVATE ${EXAERO_REPO_ROOT}/src)
target_link_libraries(exaero_standalone_public_check PRIVATE exaero_mdspan)
# Deliberately NO target_link_libraries against exaero/Kokkos: keep the include
# path clean.
add_dependencies(exaero_standalone_public_check exaero_impl)

# Property-Based and Fuzz Test Runner
add_executable(exaero_fuzz_test test_fuzzer.cpp)
target_include_directories(exaero_fuzz_test
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl)
target_link_libraries(exaero_fuzz_test PRIVATE exaero Kokkos::kokkos)
add_test(NAME ExaeroFuzzTests COMMAND exaero_fuzz_test)

# Fortran C-Interoperability Integration Test
add_executable(exaero_fortran_test test_fortran_interface.f90)
target_include_directories(exaero_fortran_test
                           PRIVATE ${EXAERO_FORTRAN_MODULE_DIR})
target_link_libraries(exaero_fortran_test PRIVATE exaero)
add_test(NAME ExaeroFortranTest COMMAND exaero_fortran_test)
# The moments round-trip loads the runtime-file fixture by absolute path; pass
# the data directory through the environment so the test stays hermetic under
# CTest.
set_tests_properties(
  ExaeroFortranTest
  PROPERTIES ENVIRONMENT
             "EXAERO_TEST_DATA_DIR=${CMAKE_CURRENT_SOURCE_DIR}/data")

# Generator determinism gate re-run the baked-in table generator in --check mode
# and fail on ANY diff to the committed header. This is the CI-side twin of the
# build-time custom command (which rewrites on drift); here we assert the
# committed artifact is exactly reproducible from the pinned snapshot, so a
# stale/hand-edited MieTableData.hpp blocks the merge rather than shipping
# non-deterministic physics.
add_test(
  NAME ExaeroMieDeterminism
  COMMAND
    ${Python3_EXECUTABLE} ${EXAERO_REPO_ROOT}/tools/generate_mie_tables.py
    --check ${EXAERO_REPO_ROOT}/src/src_impl/data/generated/MieTableData.hpp
    ${EXAERO_REPO_ROOT}/tools/geosmie_snapshot)
set_tests_properties(
  ExaeroMieDeterminism PROPERTIES PASS_REGULAR_EXPRESSION
                                  "deterministic \\(check mode")

# MKPP chemistry gtest suite (added by 963d356). Gated with the chemistry
# driver: test_registry/test_dispatcher include exaero/registry.hpp and
# exaero/dispatcher.hpp, which 37f0ce5 removed, so they cannot build until that
# workstream's sources are restored. Nothing here was deleted.
if(EXAERO_BUILD_CHEMISTRY)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https:
    //github.com/google/googletest.git
    GIT_TAG release-1.12.1)
  set(gtest_force_shared_crt
      ON
      CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  add_executable(test_registry test_registry.cpp)
  target_link_libraries(test_registry PRIVATE exaero_chem gtest_main)
  add_test(NAME RegistryTest COMMAND test_registry)

  add_executable(test_sorter test_sorter.cpp)
  target_link_libraries(test_sorter PRIVATE exaero_chem gtest_main)
  add_test(NAME SorterTest COMMAND test_sorter)

  add_executable(test_dispatcher test_dispatcher.cpp)
  target_link_libraries(test_dispatcher PRIVATE exaero_chem gtest_main)
  add_test(NAME DispatcherTest COMMAND test_dispatcher)

  add_executable(test_state test_state.cpp)
  target_link_libraries(test_state PRIVATE exaero_chem gtest_main)
  add_test(NAME StateTest COMMAND test_state)
  add_executable(test_auxiliary_engines test_auxiliary_engines.cpp)
  target_link_libraries(test_auxiliary_engines PRIVATE exaero_chem gtest_main)
  add_test(NAME AuxiliaryEnginesTest COMMAND test_auxiliary_engines)
endif()
