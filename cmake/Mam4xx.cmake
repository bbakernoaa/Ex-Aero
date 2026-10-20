# Build only the pinned MAM4xx library, using the host's Kokkos and EKAT.
# Upstream's standalone driver fetches its own EKAT/Kokkos and validation data;
# those are deliberately outside this library's dependency graph.
if(TARGET mam4xx)
  set(EXAERO_MAM4XX_TARGET mam4xx)
else()
  if(NOT TARGET ekat::KokkosUtils)
    find_package(Ekat CONFIG REQUIRED)
  endif()
  if(NOT TARGET ekat::KokkosUtils)
    message(
      FATAL_ERROR
        "FATAL ERROR: MAM4xx requires EKAT with Kokkos utilities enabled")
  endif()
  set(EXAERO_MAM4XX_SOURCE_DIR
      "${EXAERO_REPO_ROOT}/extern/mam4xx"
      CACHE PATH "Pinned MAM4xx source checkout")
  if(NOT EXISTS "${EXAERO_MAM4XX_SOURCE_DIR}/src/mam4xx/mam4_config.hpp.in")
    message(
      FATAL_ERROR
        "FATAL ERROR: initialize extern/mam4xx with git submodule update --init --recursive"
    )
  endif()
  set(MAM4XX_REAL_TYPE double)
  set(MAM4XX_DOUBLE_PRECISION 1)
  set(MAM4XX_NUM_VERTICAL_LEVELS 72)
  set(MAM4XX_ENABLE_GPU OFF)
  if(Kokkos_ENABLE_CUDA
     OR Kokkos_ENABLE_HIP
     OR Kokkos_ENABLE_SYCL)
    set(MAM4XX_ENABLE_GPU ON)
  endif()
  set(mam4xx_generated_dir "${CMAKE_CURRENT_BINARY_DIR}/mam4xx_generated")
  configure_file("${EXAERO_MAM4XX_SOURCE_DIR}/src/mam4xx/mam4_config.hpp.in"
                 "${mam4xx_generated_dir}/mam4xx/mam4_config.hpp" @ONLY)
  add_library(
    exaero_mam4xx STATIC
    "${EXAERO_MAM4XX_SOURCE_DIR}/src/mam4xx/aero_modes.cpp"
    "${EXAERO_MAM4XX_SOURCE_DIR}/src/mam4xx/aero_species.cpp"
    "${EXAERO_MAM4XX_SOURCE_DIR}/src/mam4xx/gas_species.cpp")
  target_include_directories(
    exaero_mam4xx SYSTEM PUBLIC "${EXAERO_MAM4XX_SOURCE_DIR}/src"
                                "${mam4xx_generated_dir}")
  target_link_libraries(exaero_mam4xx PUBLIC ekat::KokkosUtils Kokkos::kokkos)
  set(EXAERO_MAM4XX_TARGET exaero_mam4xx)
endif()
add_library(exaero_mam4xx_dependency INTERFACE)
target_link_libraries(exaero_mam4xx_dependency
                      INTERFACE ${EXAERO_MAM4XX_TARGET})
