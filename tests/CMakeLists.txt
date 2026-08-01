include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.12.1
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test_registry test_registry.cpp)
target_link_libraries(test_registry PRIVATE exaero gtest_main)
add_test(NAME RegistryTest COMMAND test_registry)

add_executable(test_sorter test_sorter.cpp)
target_link_libraries(test_sorter PRIVATE exaero gtest_main)
add_test(NAME SorterTest COMMAND test_sorter)

add_executable(test_dispatcher test_dispatcher.cpp)
target_link_libraries(test_dispatcher PRIVATE exaero gtest_main)
add_test(NAME DispatcherTest COMMAND test_dispatcher)

add_executable(test_state test_state.cpp)
target_link_libraries(test_state PRIVATE exaero gtest_main)
add_test(NAME StateTest COMMAND test_state)
add_executable(test_auxiliary_engines test_auxiliary_engines.cpp)
target_link_libraries(test_auxiliary_engines PRIVATE exaero gtest_main)
add_test(NAME AuxiliaryEnginesTest COMMAND test_auxiliary_engines)
