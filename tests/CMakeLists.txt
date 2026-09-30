# GOCART Diagnostics and YAML Parsing Test Runner
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
