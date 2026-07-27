# NWS Office of Modeling and Development - CMake Guidelines

## 1. Scope & Target
* **Applies to:** `**/CMakeLists.txt`, `**/*.cmake`
* **Purpose:** Ensure stable, cross-platform, HPC-ready builds that bridge C++, Fortran, and external libraries (e.g., `mam4xx`) while adhering to NOAA UFS operational standards.

## 2. Core CMake Standards
* **Minimum Version:** Always specify a modern CMake minimum version (e.g., `cmake_minimum_required(VERSION 3.24)`) capable of handling modern Fortran/C++ interoperability.
* **Target-Based CMake:** Never use directory-level includes or link directories (`include_directories()`, `link_directories()`). Always use target-specific commands (`target_include_directories()`, `target_link_libraries()`, `target_compile_options()`).
* **Scope Visibility:** Explicitly define the scope (`PUBLIC`, `PRIVATE`, `INTERFACE`) for all target properties.

## 3. Language & Compiler Standards
* **C++ Standard:** Enforce C++23 for modern interoperability (e.g., `set_target_properties(target PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON)`).
* **Fortran Support:** When enabling Fortran alongside C++, ensure modules are generated in the correct directories (`set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)`).
* **Interoperability:** Ensure the CMake configuration supports seamless linking of Fortran objects and C++ objects, adhering to the EE2 and HPC requirements.

## 4. High-Performance Computing (HPC) & Dependencies
* **MPI Integration:** Always use `find_package(MPI REQUIRED)` and link against the imported targets (`MPI::MPI_CXX`, `MPI::MPI_Fortran`) rather than assuming compiler wrappers (like `mpicxx`) handle it automatically.
* **External Modules:** Manage local dependencies in `extern/` (like `mam4xx`) carefully. Avoid polluting the global CMake cache with downstream variables.
* **OpenMP:** If utilizing shared-memory parallelism, find and link OpenMP properly (`find_package(OpenMP REQUIRED)` and target link `OpenMP::OpenMP_CXX`).

## 5. Security & Build Environments
* **No Hardcoded Paths:** Never hardcode absolute paths for libraries, includes, or compilers. Always use relative paths based on `${CMAKE_CURRENT_SOURCE_DIR}` or rely on standard package locators.
* **Out-of-Source Builds:** Build systems should strictly enforce out-of-source builds to avoid cluttering the repository.
