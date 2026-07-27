# GOCART Diagnostics and YAML Parsing Test Runner
add_executable(exaero_test_runner test_main.cpp)
target_include_directories(exaero_test_runner PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl
)
target_link_libraries(exaero_test_runner PRIVATE exaero Kokkos::kokkos)
add_test(NAME ExaeroUnitTests COMMAND exaero_test_runner)

# Property-Based and Fuzz Test Runner
add_executable(exaero_fuzz_test test_fuzzer.cpp)
target_include_directories(exaero_fuzz_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl
)
target_link_libraries(exaero_fuzz_test PRIVATE exaero Kokkos::kokkos)
add_test(NAME ExaeroFuzzTests COMMAND exaero_fuzz_test)

# Fortran C-Interoperability Integration Test
add_executable(exaero_fortran_test test_fortran_interface.f90)
target_include_directories(exaero_fortran_test PRIVATE
    ${CMAKE_BINARY_DIR}/src
)
target_link_libraries(exaero_fortran_test PRIVATE exaero)
add_test(NAME ExaeroFortranTest COMMAND exaero_fortran_test)
