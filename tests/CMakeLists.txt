# GOCART Diagnostics and YAML Parsing Test Runner
add_executable(exaero_test_runner test_main.cpp)
target_include_directories(exaero_test_runner PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl
)
target_link_libraries(exaero_test_runner PRIVATE exaero Kokkos::kokkos)
add_test(NAME ExaeroUnitTests COMMAND exaero_test_runner)

# Note: exaero_memory_mapper_test is disabled on macOS Apple Clang 
# due to an adjacent template parameter pack bug inside Xcode's libc++ Kokkos headers
# when combining privately vendored mdspan with Kokkos_Core.hpp in the same translation unit.
# Zero-copy memory mapping is verified directly via our public GocartPackage integration tests.
#
# add_executable(exaero_memory_mapper_test test_memory_mapper.cpp)
# target_include_directories(exaero_memory_mapper_test PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl
# )
# target_link_libraries(exaero_memory_mapper_test PRIVATE exaero Kokkos::kokkos)
# add_test(NAME ExaeroMemoryMapperTest COMMAND exaero_memory_mapper_test)
