cmake_minimum_required(VERSION 3.24)
project(ExaeroDriver CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Kokkos REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(EXAERO_MECHANISMS "ChapmanCycle" CACHE STRING "Semicolon separated list of mechanisms to compile")

option(EXAERO_WITH_CLOUDJ "Enable Cloud-J photolysis forcing support" OFF)
option(EXAERO_WITH_ISORROPIALITE "Enable ISORROPIA-Lite thermodynamic support" OFF)

set(REGISTRY_SRC "${CMAKE_CURRENT_BINARY_DIR}/registry_generated.cpp")
add_custom_command(
    OUTPUT ${REGISTRY_SRC}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate_registry.py ${REGISTRY_SRC} "${EXAERO_MECHANISMS}"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/generate_registry.py
    COMMENT "Generating mechanism registry for: ${EXAERO_MECHANISMS}"
)

add_library(exaero STATIC
    src/api.cpp
    ${REGISTRY_SRC}
    src/sorter.cpp
    src/dispatcher.cpp
    src/auxiliary_engines.cpp
)

set(MKPP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../build/mkpp-generated" CACHE PATH "Path to MKPP generated headers")

target_include_directories(exaero PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${MKPP_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(exaero PUBLIC Kokkos::kokkos)

if(EXAERO_WITH_CLOUDJ)
    # T004: Instead of Fortran linkage we just use the native C++ port headers
    target_include_directories(exaero PRIVATE ${CMAKE_SOURCE_DIR}/../Cloud-J/include)
    target_compile_definitions(exaero PUBLIC EXAERO_WITH_CLOUDJ)
endif()

if(EXAERO_WITH_ISORROPIALITE)
    # The ISORROPIA-Lite refactor is vectorized C++ according to the master plan (V8.md section 2.6.2)
    # Let's verify what the src/isorropialite directory contains
    add_subdirectory(${CMAKE_SOURCE_DIR}/../isorropialite ${CMAKE_BINARY_DIR}/isorropialite)
    target_link_libraries(exaero PRIVATE isorropia)
    target_compile_definitions(exaero PUBLIC EXAERO_WITH_ISORROPIALITE)
endif()

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
if(BUILD_TESTING)
    add_subdirectory(${CMAKE_SOURCE_DIR}/../../tests ${CMAKE_BINARY_DIR}/e2e_tests)
endif()

# Add dummy implementations for tests to link against if needed
