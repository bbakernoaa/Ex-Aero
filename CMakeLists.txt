cmake_minimum_required(VERSION 3.18)
project(EX-aero LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Kokkos dependency
find_package(Kokkos REQUIRED)

add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
