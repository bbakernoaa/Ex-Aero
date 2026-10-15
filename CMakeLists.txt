cmake_minimum_required(VERSION 3.24)
project(ExaeroDriver CXX Fortran)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Resolve this repository's root from the current list file, not from the
# top-level project.  When EX-aero is consumed via add_subdirectory from a
# parent build (e.g. CATChem), CMAKE_SOURCE_DIR points at the PARENT, so the
# build-time MIE generator and the test fixtures would resolve against the wrong
# tree.  CMAKE_CURRENT_SOURCE_DIR here is always this repo's root because this
# list file lives at that root, in standalone and subproject mode alike.  Every
# in-repo absolute path below is anchored to this variable.
set(EXAERO_REPO_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# Single documented home for the Fortran module files (exaero_interface.mod),
# defined at this level so both the library scope (src/) and the test scope
# (tests/) resolve it, and consuming hosts get it through the exaero target's
# PUBLIC include interface.
set(EXAERO_FORTRAN_MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/src/fortran_modules)

# Adopt a Kokkos the parent already located instead of re-running find_package
# against a possibly different on-disk install (which hard-fails in the
# exported-targets file with "Some (but not all) targets already defined"). Same
# guard MKPP uses.
if(NOT TARGET Kokkos::kokkos)
  find_package(Kokkos REQUIRED)
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# ---------------------------------------------------------------------------
# Build-graph note (2026-09-09)
#
# Two independent target groups in this repository both used the name `exaero`:
# * the aerosol interface library (src/, src/src_impl/, include) built via
# add_subdirectory(src) below -- SHARED `exaero` plus the targets `exaero_impl`
# (mdspan, no Kokkos) and `exaero_solver` (Kokkos, no mdspan); * the MKPP
# chemistry driver (top-level add_library), generated from generate_registry.py
# + include/exaero/registry.hpp.
#
# Commit d3f976d replaced the aerosol-only top-level with the chemistry driver
# and dropped add_subdirectory(src), so the aerosol sources targeted by were no
# longer compiled or tested. The chemistry driver is itself currently
# unbuildable: generate_registry.py, include/exaero/registry.hpp,
# include/exaero/dispatcher.hpp, src/dispatcher.cpp and src/registry.cpp were
# removed by 37f0ce5, src/sorter.cpp is empty, and the bare
# add_custom_command(OUTPUT ${REGISTRY_SRC}) expanded to zero arguments and
# aborted configuration outright. It is therefore gated OFF by default and its
# target renamed `exaero_chem` to free the `exaero` name for the interface
# library -- restoring those inputs is a one-line re-enable. No chemistry source
# file was deleted by this change.
# ---------------------------------------------------------------------------

option(
  EXAERO_BUILD_CHEMISTRY
  "Build the MKPP chemistry driver (requires generate_registry.py + registry/dispatcher headers)"
  OFF)

option(
  EXAERO_BUILD_TESTS
  "Build the EX-aero test runners (unit, fuzz, Fortran interop, MIE determinism)"
  ON)

# --- Aerosol interface library: public C++/C/Fortran boundary + impl/solver ---
add_subdirectory(src)

if(EXAERO_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(EXAERO_BUILD_CHEMISTRY)
  set(EXAERO_MECHANISMS
      "ChapmanCycle"
      CACHE STRING "Semicolon separated list of mechanisms to compile")

  option(EXAERO_WITH_CLOUDJ "Enable Cloud-J photolysis forcing support" OFF)
  option(EXAERO_WITH_ISORROPIALITE
         "Enable ISORROPIA-Lite thermodynamic support" OFF)

  set(REGISTRY_SRC "${CMAKE_CURRENT_BINARY_DIR}/registry_generated.cpp")
  add_custom_command(
    OUTPUT ${REGISTRY_SRC}
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate_registry.py
      ${REGISTRY_SRC} "${EXAERO_MECHANISMS}"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generate_registry.py
    COMMENT "Generating mechanism registry for: ${EXAERO_MECHANISMS}")

  add_library(exaero_chem STATIC src/api.cpp ${REGISTRY_SRC} src/sorter.cpp
                                 src/auxiliary_engines.cpp)

  set(MKPP_INCLUDE_DIR
      "${CMAKE_SOURCE_DIR}/../../build/mkpp-generated"
      CACHE PATH "Path to MKPP generated headers")

  target_include_directories(
    exaero_chem
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
           $<BUILD_INTERFACE:${MKPP_INCLUDE_DIR}> $<INSTALL_INTERFACE:include>)

  target_link_libraries(exaero_chem PUBLIC Kokkos::kokkos)

  if(EXAERO_WITH_CLOUDJ)
    # : Instead of Fortran linkage we just use the native C++ port headers
    target_include_directories(exaero_chem
                               PRIVATE ${CMAKE_SOURCE_DIR}/../Cloud-J/include)
    target_compile_definitions(exaero_chem PUBLIC EXAERO_WITH_CLOUDJ)
  endif()

  if(EXAERO_WITH_ISORROPIALITE)
    # The ISORROPIA-Lite refactor is vectorized C++ according to the master plan
    # (V8.md section 2.6.2)
    add_subdirectory(${CMAKE_SOURCE_DIR}/../isorropialite
                     ${CMAKE_BINARY_DIR}/isorropialite)
    target_link_libraries(exaero_chem PRIVATE isorropia)
    target_compile_definitions(exaero_chem PUBLIC EXAERO_WITH_ISORROPIALITE)
  endif()

  if(EXISTS ${CMAKE_SOURCE_DIR}/../../tests)
    add_subdirectory(${CMAKE_SOURCE_DIR}/../../tests
                     ${CMAKE_BINARY_DIR}/e2e_tests)
  endif()
endif()
