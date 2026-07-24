# Private Header Target for mdspan backport
add_library(exaero_mdspan INTERFACE)
target_include_directories(exaero_mdspan INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party/mdspan/include>
)

# Implementation Library
add_library(exaero_impl STATIC dummy.cpp)
target_include_directories(exaero_impl PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)
target_link_libraries(exaero_impl PUBLIC exaero_mdspan)

# Link exaero to private implementation
target_link_libraries(exaero PRIVATE exaero_impl)
