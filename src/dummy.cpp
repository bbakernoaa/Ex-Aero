/// @file dummy.cpp
/// @brief Placeholder translation unit so the @c exaero shared-library
/// target has at least one object before the interface is populated.
///
/// Historically needed to satisfy CMake until IAerosolPackage was
/// implemented; kept as a stable link anchor for the public library.
void exaero_dummy_function() {}
