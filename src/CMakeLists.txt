# Public Interface Target (builds C++, C, and Fortran interfaces)
add_library(exaero SHARED dummy.cpp exaero/IAerosolPackage_C.cpp
                          exaero/exaero_interface.f90)

# Explicit, documented home for the generated Fortran module files
# (exaero_interface.mod).  Relying on the compiler's implicit default makes the
# location mode-dependent (it shifts under add_subdirectory) and gives consuming
# hosts nothing stable to add to their Fortran include path.  The directory
# variable itself is defined at the repository top level.
set_target_properties(exaero PROPERTIES Fortran_MODULE_DIRECTORY
                                        ${EXAERO_FORTRAN_MODULE_DIR})
target_include_directories(
  exaero PUBLIC $<BUILD_INTERFACE:${EXAERO_FORTRAN_MODULE_DIR}>
                $<INSTALL_INTERFACE:include>)

target_include_directories(
  exaero PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:include>)

target_include_directories(exaero PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src_impl)

target_link_libraries(exaero PUBLIC Kokkos::kokkos)

add_subdirectory(src_impl)
