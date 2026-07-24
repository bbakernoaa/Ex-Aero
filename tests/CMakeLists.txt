# GOCART Diagnostics and YAML Parsing Test Runner
add_executable(exaero_test_runner test_main.cpp)
target_include_directories(exaero_test_runner PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/src_impl
)
target_link_libraries(exaero_test_runner PRIVATE exaero Kokkos::kokkos)
add_test(NAME ExaeroUnitTests COMMAND exaero_test_runner)
